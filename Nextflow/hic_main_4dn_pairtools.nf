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
params.chrom_sizes  = "/home/dhllove/NGS_DATA/reference_FASTA/chrom.sizes"
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
	
    pairtools parse \\
        --chroms-path ${chrom_sizes} \\
        --walks-policy 5unique \\
        --add-columns mapq \\
        --output ${sample_id}.parsed.pairsam.gz \\
        ${bam}

    pairtools sort --nproc ${task.cpus} \\
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
    val resolution   // e.g. 100000

    output:
    path "${sample_id}.compartments.bedgraph", emit: compartments

    script:
    """
    python3 <<'PYEOF'
import bioframe, cooler, cooltools

clr = cooler.Cooler("${mcool}::resolutions/${resolution}")
bins = clr.bins()[:]
view_df = bioframe.make_viewframe(clr.chromsizes)
gc_track = bioframe.frac_gc(bins[["chrom", "start", "end"]], bioframe.load_fasta("${fasta}"))

_, eigvecs = cooltools.eigs_cis(clr, gc_track, view_df=view_df, n_eigs=3)
eigvecs[["chrom", "start", "end", "E1"]].to_csv(
    "${sample_id}.compartments.bedgraph", sep="\\t", header=False, index=False
)
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
    tuple val(sample_id), path("*"), emit: tads

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 8
    """
    echo "=============================================="
    echo "Juicer Arrowhead (TAD calling)"
    echo "Sample ID : ${sample_id}"
    echo "=============================================="

    java -Xmx${avail_mem}g -jar /opt/juicer/juicer_tools.jar arrowhead \\
        -m 2000 \\
        ${hic} \\
        .
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
    tuple val(sample_id), path("*"), emit: loops

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 8
    """
    echo "=============================================="
    echo "Juicer HiCCUPS (Loop calling)"
    echo "Sample ID : ${sample_id}"
    echo "=============================================="

    java -Xmx${avail_mem}g -jar /opt/juicer/juicer_tools.jar hiccups \\
        ${hic} \\
        . || echo "HiCCUPS completed with status/skipped (Requires GPU or specific resolution)."
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
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
    ch_chrom_sizes = Channel.fromPath(params.chrom_sizes, checkIfExists: true)

    // 1. Alignment
    ALIGN_BWA(
        ch_reads,
        params.fasta
    )
    RUN_PAIRTOOLS(ALIGN_BWA.out.bam, ch_chrom_sizes)
    RUN_JUICER_PRE(RUN_PAIRTOOLS.out.valid_pairs, ch_chrom_sizes)
    BUILD_COOLER(RUN_PAIRTOOLS.out.valid_pairs, ch_chrom_sizes, params.base_resolution)

    RUN_COMPARTMENTS(BUILD_COOLER.out.mcool, ch_fasta, 100000)

    // 4. 並列執行：將 .hic 同時餵給 ARROWHEAD 與 HICCUPS找TADs & Loops
    RUN_ARROWHEAD(RUN_JUICER_PRE.out.hic_matrix)
    RUN_HICCUPS(RUN_JUICER_PRE.out.hic_matrix)

}