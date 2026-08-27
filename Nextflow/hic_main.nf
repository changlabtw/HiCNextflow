nextflow.enable.dsl=2

/*
========================================================================================
    3D Genome & Metagenomic Hi-C Pipeline
========================================================================================
*/

// 定義輸入參數與預設值
params.fastq_dir    = "/home/dhllove/NGS_DATA/HiC"
params.reads        = "${params.fastq_dir}/*_{1,2}.fastq"
params.fasta        = "/home/dhllove/NGS_DATA/reference_FASTA/Homo_sapiens_assembly38.fasta"
params.outdir       = "/home/dhllove/Work/HiC"
params.chrom_sizes  = "/home/dhllove/NGS_DATA/reference_FASTA/chrom.sizes"
params.restriction_bed = "/home/dhllove/NGS_DATA/reference_FASTA/Homo_sapiens_assembly38_hindiii.bed"
params.sif1         = "/home/dhllove/DockerImage/bin3c_2.sif" // 指向包含工具鏈的 SIF 檔案
params.sif2         = "/home/dhllove/DockerImage/hicpro.sif" // 指向包含工具鏈的 SIF 檔案
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
    tuple val(sample_id), path("${sample_id}_aligned.bam"), path("${sample_id}_aligned.bam.bai"), emit: bam

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
    samtools sort -@ ${samtools_sort_threads} -m 2G -o ${sample_id}_aligned.bam -
    
    samtools index ${sample_id}_aligned.bam
    """
}

// 2. bin3C Processing (修復：在 input 明確拆出 sample_id 變數)
process RUN_BIN3C {
    tag "Sample: ${sample_id}"
    container "${params.sif1}"
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "${sample_id}/bin3C/${filename}" }

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    path "bin3C_out/*", emit: bin3c_results

    script:
    """
	echo "=============================================="
	echo "bin3C"
	echo "Sample ID : ${sample_id}"
	echo "BAM : ${bam}"
	echo "=============================================="
	
    mkdir -p bin3C_out
    bin3C.py mkmap -b ${bam} -o bin3C_out/map.h5
    bin3C.py cluster -m bin3C_out/map.h5 -o bin3C_out
    """
}

process RUN_HICPRO {
    tag "Sample: ${sample_id}"
    container "${params.sif2}" // 假設 sif2 含有 HiC-Pro
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "${sample_id}/hicpro/${filename}" }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path restriction_bed
    path chrom_sizes

    output:
    tuple val(sample_id), path("${sample_id}.validPairs"), emit: valid_pairs
    tuple val(sample_id), path("${sample_id}.merged_nodups.txt.gz"), emit: merged_nodups

    script:
    """
    echo "=============================================="
    echo "HiC-Pro (BAM -> validPairs -> Juicebox format)"
    echo "Sample ID : ${sample_id}"
    echo "=============================================="

    # Step 3.1: 使用 mapped_2hic_fragments.py 計算片段並篩選 validPairs
    python /opt/HiC-Pro/bin/utils/mapped_2hic_fragments.py \\
        -f ${restriction_bed} \\
        -r ${bam} \\
        -o . \\
        -v

    # 確保產生的檔名符合 Sample ID
    if [ -f "out.validPairs" ]; then
        mv out.validPairs ${sample_id}.validPairs
    fi

    # Step 3.2: 使用 hicpro2juicebox.sh 將 .validPairs 轉為 Juicer 格式
    /opt/HiC-Pro/bin/utils/hicpro2juicebox.sh \\
        -i ${sample_id}.validPairs \\
        -g ${chrom_sizes} \\
        -o .

    if [ -f "merged_nodups.txt" ]; then
        gzip -c merged_nodups.txt > ${sample_id}.merged_nodups.txt.gz
    elif [ -f "merged_nodups.txt.gz" ]; then
        mv merged_nodups.txt.gz ${sample_id}.merged_nodups.txt.gz
    fi
    """
}

// 3. Juicer Tools Processing (修復：在 input 明確拆出 sample_id 變數)
process RUN_JUICER_PRE {
    tag "Sample: ${sample_id}"
    container "${params.sif3}"
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "${sample_id}/juicer/${filename}" }

    input:
    tuple val(sample_id), path(merged_nodups)
    path chrom_sizes

    output:
    path "${sample_id}.hic", emit: hic_matrix
    path "loops/*", optional: true, emit: loops
    path "tads/*", optional: true, emit: tads

    script:
    def avail_mem = task.memory ? task.memory.toGiga() : 8
    """
    echo "=============================================="
    echo "Juicer Tools (.hic generation & Loop/TAD calling)"
    echo "Sample ID : ${sample_id}"
    echo "=============================================="

    # Step 4.1: 排序轉譯後的短格式數據
    zcat ${merged_nodups} | \
    sort -k2,2d -k6,6d -S ${avail_mem}G --parallel=${task.cpus} > merged_nodups_sorted.txt

    # Step 4.2: 使用 Juicer Tools pre 建立 .hic 矩陣
    java -Xmx${avail_mem}g -jar /opt/juicer/juicer_tools.jar pre \\
        merged_nodups_sorted.txt \\
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
       3 D   G E N O M E   P I P E L I N E (DSL2) - Local Node Mode
      ================================================================
       Fastq Directory  : ${params.fastq_dir}
       Reads Pattern    : ${params.reads}
       Reference FASTA  : ${params.fasta}
       Chrom Sizes File : ${params.chrom_sizes}
       Restriction BED  : ${params.restriction_bed}
       Output Dir       : ${params.outdir}
      ================================================================
   """
    ch_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    //ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
    ch_chrom_sizes = Channel.fromPath(params.chrom_sizes, checkIfExists: true)
    ch_restriction_bed = Channel.fromPath(params.restriction_bed, checkIfExists: true)

    // 1. Alignment
    ALIGN_BWA(
        ch_reads,
        params.fasta
    )

    // 2. 下游分流 (bin3C 與 HiC-Pro 同時平行處理)
    RUN_BIN3C(ALIGN_BWA.out.bam)
    RUN_HICPRO(ALIGN_BWA.out.bam, ch_restriction_bed, ch_chrom_sizes)

    // 3. Juicer Pre (生成 .hic 矩陣)
    RUN_JUICER_PRE(RUN_HICPRO.out.merged_nodups, ch_chrom_sizes)

    // 4. 並列執行：將 .hic 同時餵給 ARROWHEAD 與 HICCUPS找TADs & Loops
    RUN_ARROWHEAD(RUN_JUICER_PRE.out.hic_matrix)
    RUN_HICCUPS(RUN_JUICER_PRE.out.hic_matrix)

}