nextflow.enable.dsl=2

/*
========================================================================================
    3D Genome & Metagenomic Hi-C Pipeline
========================================================================================
*/

// 定義輸入參數與預設值
params.fastq_dir    = "/home/dhllove/NGS_DATA/HiC"
params.reads        = "${params.fastq_dir}/4DN/*_R{1,2}_001.fastq"
params.fasta        = "/home/dhllove/NGS_DATA/reference_FASTA/Homo_sapiens_assembly38.fasta"
params.outdir       = "/home/dhllove/Work/HiC"
params.chrom_sizes  = "/home/dhllove/NGS_DATA/reference_FASTA/chrom.sizes.primary"
params.base_resolution = 10000
params.cooler_resolutions = "10000,20000,50000,100000,250000,500000,1000000"
params.sif1         = "/home/dhllove/DockerImage/bin3c_2.sif" // 指向包含工具鏈的 SIF 檔案
params.sif2         = "/home/dhllove/DockerImage/pairtools.sif" // 指向包含工具鏈的 SIF 檔案
params.sif3         = "/home/dhllove/DockerImage/juicer.sif" // 指向包含工具鏈的 SIF 檔案
params.sif          = params.sif1

/*
========================================================================================
    PROCESS DEFINITIONS (全數指定容器執行環境)
========================================================================================
*/

process ALIGN_BWA {
    tag "Sample: ${sample_id}"
    container "${params.sif1}"
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "${sample_id}/aligned_bam/${filename}" }

    input:
    tuple val(sample_id), path(reads)
    val fasta

    output:
    tuple val(sample_id), path("${sample_id}_aligned.bam"), emit: bam

    script:
    def samtools_sort_threads = Math.max(1, Math.floor(task.cpus * 0.25).toInteger())
    def bwa_threads = Math.max(1, task.cpus - samtools_sort_threads - 2)    

    """
	echo "=============================================="
	echo "BWA Alignment" 
	echo "Sample ID : ${sample_id}"
	echo "Read 1 : ${reads[0]}"
	echo "Read 2 : ${reads[1]}"
	echo "=============================================="

    bwa-mem2 mem -5SP -t ${bwa_threads} ${fasta} ${reads[0]} ${reads[1]} | \\
    samtools view -S -b -@ 4 - | \\
    samtools sort -n -@ ${samtools_sort_threads} -m 2G -o ${sample_id}_aligned.bam -

    """
}


process RUN_PAIRTOOLS {
    tag "Sample: ${sample_id}"
    container "${params.sif2}"   
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "${sample_id}/pairtools/${filename}" }

    input:
    tuple val(sample_id), path(bam)
    path chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.valid.pairs.gz"), emit: valid_pairs

    script:
    """
	echo "======================================================"
	echo "Pairtools" 
	echo "Sample ID : ${sample_id}"
	echo "parses alignments into Hi-C pairs  --> Sorting --> "
	echo "PCR-duplicate removal --> pair-type filtering"
	echo "======================================================"
	export PT_TMPDIR=\$(mktemp -d /tmp/pairtools_${sample_id}.XXXXXX)
    trap 'rm -rf "\$PT_TMPDIR"' EXIT
	
    pairtools parse \\
        --chroms-path ${chrom_sizes} \\
        --walks-policy 5unique \\
        --add-columns mapq \\
        --output ${sample_id}.parsed.pairsam.gz \\
        ${bam}

    pairtools sort --nproc ${task.cpus} \\
        --tmpdir \$PT_TMPDIR \\
        --output ${sample_id}.sorted.pairsam.gz \\
        ${sample_id}.parsed.pairsam.gz

    pairtools dedup --mark-dups \\
        --output ${sample_id}.dedup.pairsam.gz \\
        --output-stats ${sample_id}.dedup.stats \\
        ${sample_id}.sorted.pairsam.gz

    pairtools select \\
        '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \\
        --output ${sample_id}.valid.pairs.gz \\
        ${sample_id}.dedup.pairsam.gz
    """
}
process BUILD_COOLER {
    tag "Sample: ${sample_id}"
    container "${params.sif2}"
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "${sample_id}/cooler/${filename}" }

    input:
    tuple val(sample_id), path(pairs)
    path chrom_sizes
    val base_resolution   // e.g. 10000

    output:
    tuple val(sample_id), path("${sample_id}.mcool"), emit: mcool

    script:
    """
    cooler cload pairs \\
        -c1 2 -p1 3 -c2 4 -p2 5 \\
        ${chrom_sizes}:${base_resolution} \\
        ${pairs} \\
        ${sample_id}.${base_resolution}.cool

    cooler zoomify \\
        --nproc ${task.cpus} \\
        --resolutions ${params.cooler_resolutions} \\
        --balance \\
        -o ${sample_id}.mcool \\
        ${sample_id}.${base_resolution}.cool
    """
}

process RUN_COMPARTMENTS {
    tag "Sample: ${sample_id}"
    container "${params.sif2}"
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "${sample_id}/compartments/${filename}" }

    input:
    tuple val(sample_id), path(mcool)
    path fasta
    val resolution      // e.g. 100000
    val target_chroms   // comma-separated list, e.g. "chr1,chr2,...,chr22" — which chromosomes to plot in the track figure

    output:
    path "${sample_id}.compartments.bedgraph", emit: compartments
    path "${sample_id}_track_*.png", optional: true, emit: track_plots
    path "${sample_id}_demo_track.png", optional: true, emit: track_demo_plot
    path "${sample_id}_saddle.png", optional: true, emit: saddle_plot
    path "${sample_id}_saddle.npz", optional: true, emit: saddle_data

    script:
    """
    python3 <<'PYEOF'
import bioframe
import cooler
import cooltools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# 1. 計算 E1 Eigenvector (Compartment Score)
# -------------------------------------------------------------------
clr = cooler.Cooler("${mcool}::resolutions/${resolution}")
bins = clr.bins()[:]
view_df = bioframe.make_viewframe(clr.chromsizes)
gc_track = bioframe.frac_gc(bins[["chrom", "start", "end"]], bioframe.load_fasta("${fasta}"))

# 導出 E1 數值
_, eigvecs = cooltools.eigs_cis(clr, gc_track, view_df=view_df, n_eigs=3)
track_df = eigvecs[["chrom", "start", "end", "E1"]].copy()

# 儲存 BedGraph
track_df.to_csv(
    "${sample_id}.compartments.bedgraph", sep="\\t", header=False, index=False
)

# -------------------------------------------------------------------
# 2. 視覺化：Graphical Visualization of Continuous Score Track
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# 2a. 視覺化：每條染色體各自輸出一張全解析度 Track 圖
# -------------------------------------------------------------------
def draw_track(ax, chrom, sub_df, show_xlabel):
    pos_mb = sub_df['start'] / 1e6
    scores = sub_df['E1']
    ax.plot(pos_mb, scores, color='black', linewidth=0.5)
    ax.axhline(0, color='grey', linestyle='--', linewidth=0.8)
    ax.fill_between(pos_mb, scores, 0, where=(scores >= 0), color='crimson', alpha=0.6, label='Comp A')
    ax.fill_between(pos_mb, scores, 0, where=(scores < 0), color='steelblue', alpha=0.6, label='Comp B')
    ax.set_ylabel(f"{chrom}\\nE1 Score", fontsize=10)
    ax.set_ylim(-1.1, 1.1)
    if show_xlabel:
        ax.set_xlabel("Genomic Position (Mb)", fontsize=10)

try:
    # 可選染色體清單，由 Nextflow 參數傳入 (預設 chr1-chr22)，
    # 並與 cooler 中實際存在的染色體名稱取交集，避免命名不符造成錯誤
    requested_chroms = [c.strip() for c in "${target_chroms}".split(",") if c.strip()]
    target_chroms = [c for c in requested_chroms if c in clr.chromnames]

    missing = sorted(set(requested_chroms) - set(target_chroms))
    if missing:
        print(f"Warning: requested chromosomes not found in cooler, skipping: {missing}")

    if not target_chroms:
        # 找不到任何相符染色體時，退回前 3 條作為保底
        target_chroms = clr.chromnames[:3]

    for chrom in target_chroms:
        sub_df = track_df[track_df['chrom'] == chrom].dropna(subset=['E1'])
        fig, ax = plt.subplots(figsize=(14, 3))
        draw_track(ax, chrom, sub_df, show_xlabel=True)
        ax.set_title(f"Compartment Score Track - {chrom} (Sample: ${sample_id})", fontsize=11)
        plt.tight_layout()
        plt.savefig(f"${sample_id}_track_{chrom}.png", dpi=300)
        plt.close(fig)

    print(f"Saved {len(target_chroms)} individual per-chromosome track plot(s).")
except Exception as e:
    print(f"Failed to generate per-chromosome track plots: {e}")

# -------------------------------------------------------------------
# 2b. 視覺化：固定 3 條染色體的總覽圖 (Demo, 固定 3-row layout)
# -------------------------------------------------------------------
try:
    demo_chroms = [c for c in ['chr1', 'chr2', 'chr3', '1', '2', '3'] if c in clr.chromnames][:3]
    if not demo_chroms:
        demo_chroms = clr.chromnames[:3]

    fig, axes = plt.subplots(3, 1, figsize=(12, 6), sharex=False)
    fig.suptitle(f"Compartment Score Overview - Demo (Sample: ${sample_id})", fontsize=14, fontweight='bold')

    for idx, chrom in enumerate(demo_chroms):
        sub_df = track_df[track_df['chrom'] == chrom].dropna(subset=['E1'])
        draw_track(axes[idx], chrom, sub_df, show_xlabel=(idx == 2))

    plt.tight_layout()
    plt.savefig("${sample_id}_demo_track.png", dpi=600)
    plt.close(fig)
    print("Demo overview plot saved successfully.")
except Exception as e:
    print(f"Failed to generate demo track plot: {e}")

# -------------------------------------------------------------------
# 3. 繪製 Cooltools Saddle Plot (Compartment Strength Analysis)
# -------------------------------------------------------------------
try:
    Q_LO = 0.025   # 忽略 E1 最低 2.5% 的 bin
    Q_HI = 0.975   # 忽略 E1 最高 2.5% 的 bin
    q_bins = 20    # 其餘 95% 分成 20 個等量分位區間

    # 計算 cis expected (P(s)) — saddle() 需要用它做 O/E 正規化
    cvd = cooltools.expected_cis(clr, view_df=view_df)

    # saddle() 直接吃未離散化的 track_df (chrom,start,end,E1)，
    # 依 n_bins/qrange 在內部自動離散化 — 不需要另外呼叫 digitize()
    interaction_sum, interaction_count = cooltools.saddle(
        clr,
        cvd,
        track_df,
        'cis',
        n_bins=q_bins,
        qrange=(Q_LO, Q_HI),
        view_df=view_df,
    )

    with np.errstate(divide='ignore', invalid='ignore'):
        matrix = np.log10(interaction_sum / interaction_count)

    np.savez(
        "${sample_id}_saddle.npz",
        interaction_sum=interaction_sum,
        interaction_count=interaction_count,
    )

    # 繪製 Saddle Heatmap
    fig = plt.figure(figsize=(12, 10))
    plt.imshow(
        matrix,
        cmap='coolwarm',
        vmax=0.5, vmin=-0.5,
        extent=[0, q_bins, 0, q_bins],
        origin='lower'
    )
    cbar = plt.colorbar(label='log10 (Observed / Expected)')
    plt.title("Saddle Plot - Compartment Signal (${sample_id})", fontsize=12)
    plt.xlabel("B --> A (Quantiles)", fontsize=10)
    plt.ylabel("B --> A (Quantiles)", fontsize=10)
    
    plt.tight_layout()
    plt.savefig("${sample_id}_saddle.png", dpi=600)
    plt.close()
    print("Saddle plot saved successfully.")
except Exception as e:
    print(f"Failed to generate saddle plot: {e}")

PYEOF
    """
}

// 3. Juicer Tools Processing (修復：在 input 明確拆出 sample_id 變數)
process RUN_JUICER_PRE {
    tag "Sample: ${sample_id}"
    container "${params.sif3}"
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "${sample_id}/juicer/${filename}" }

    input:
    tuple val(sample_id), path(pairs)
    path chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.hic"), emit: hic_matrix
    path "loops/*", optional: true, emit: loops
    path "tads/*", optional: true, emit: tads

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 8
    """
    echo "=============================================="
    echo "Juicer Tools (.hic generation)"
    echo "Sample ID : ${sample_id}"
    echo "=============================================="


    # Step 4.2: 使用 Juicer Tools pre 建立 .hic 矩陣
    java -Xmx${avail_mem}g -jar /opt/juicer/juicer_tools.jar pre \\
        ${pairs} \\
        ${sample_id}.hic \\
        ${chrom_sizes}
    """
}

// JUICER平行任務 A: TADs Calling (Arrowhead)
process RUN_ARROWHEAD {
    tag "Sample: ${sample_id}"
    container "${params.sif3}"
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "${sample_id}/juicer/tads/${filename}" }

    input:
    tuple val(sample_id), path(hic)

    output:
    tuple val(sample_id), path("*", optional: true), emit: tads

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 8
    """
    echo "=============================================="
    echo "Juicer Arrowhead (TAD calling)"
    echo "Sample ID : ${sample_id}"
    echo "=============================================="

    java -Xmx${avail_mem}g -jar /opt/juicer/juicer_tools.jar arrowhead \\
        --ignore-sparsity \\
        -m 2000 \\
        ${hic} \\
        . || echo "Arrowhead failed for ${sample_id} -- see .command.log" > ${sample_id}_arrowhead_status.txt
    """
}
// JUICER平行任務 B: Loops Calling (HiCCUPS)
process RUN_HICCUPS {
    tag "Sample: ${sample_id}"
    container "${params.sif3}"
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "${sample_id}/juicer/loops/${filename}" }

    input:
    tuple val(sample_id), path(hic)

    output:
    tuple val(sample_id), path("*", optional: true), emit: loops

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 8
    """
    echo "=============================================="
    echo "Juicer HiCCUPS (Loop calling)"
    echo "Sample ID : ${sample_id}"
    echo "=============================================="
    echo "--- GPU state at job start ---"
    nvidia-smi || echo "nvidia-smi not available in container"
    echo "------------------------------"

    java -Xmx${avail_mem}g -jar /opt/juicer/juicer_tools.jar hiccups \\
        --ignore-sparsity \\
        ${hic} \\
        . || echo "HiCCUPS failed for ${sample_id} -- see .command.log" > ${sample_id}_hiccups_status.txt
    """
}
/*
========================================================================================
    WORKFLOW DEFINITION
========================================================================================
*/

workflow {
    log.info """\\
      ================================================================
       3 D   G E N O M E   P I P E L I N E (DSL2)
      ================================================================
       Fastq Directory    : ${params.fastq_dir}
       Reads Pattern      : ${params.reads}
       Reference FASTA    : ${params.fasta}
       Chrom Sizes File   : ${params.chrom_sizes}
       Output Dir         : ${params.outdir}
      ================================================================
   """
    ch_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true).first()
    ch_chrom_sizes = Channel.fromPath(params.chrom_sizes, checkIfExists: true).first()

    // 預設繪製 chr1-chr22 (體染色體)，可用 --compartment_chroms 覆寫，
    // 例如 --compartment_chroms "chr1,chrX" 只畫特定染色體
    def compartment_chroms = params.compartment_chroms ?: (1..22).collect { "chr${it}" }.join(',')

    // 1. Alignment
    ALIGN_BWA(
        ch_reads,
        params.fasta
    )
    RUN_PAIRTOOLS(ALIGN_BWA.out.bam, ch_chrom_sizes)
    RUN_JUICER_PRE(RUN_PAIRTOOLS.out.valid_pairs, ch_chrom_sizes)
    BUILD_COOLER(RUN_PAIRTOOLS.out.valid_pairs, ch_chrom_sizes, params.base_resolution)

    RUN_COMPARTMENTS(BUILD_COOLER.out.mcool, ch_fasta, 100000, compartment_chroms)

    // 4. 並列執行：將 .hic 同時餵給 ARROWHEAD 與 HICCUPS找TADs & Loops
    RUN_ARROWHEAD(RUN_JUICER_PRE.out.hic_matrix)
    RUN_HICCUPS(RUN_JUICER_PRE.out.hic_matrix)

}