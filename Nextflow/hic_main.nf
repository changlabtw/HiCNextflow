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
params.sif1         = "/home/dhllove/DockerImage/bin3c.sif" // 指向包含工具鏈的 SIF 檔案
params.sif2         = "/home/dhllove/DockerImage/juicer.sif" // 指向包含工具鏈的 SIF 檔案
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
    def samtools_threads = Math.max(1, task.cpus - 2)

    """
	echo "=============================================="
	echo "BWA Alignment" 
	echo "Sample ID : ${sample_id}"
	echo "Read 1 : ${reads[0]}"
	echo "Read 2 : ${reads[1]}"
	echo "=============================================="

    bwa mem -5SP -t ${task.cpus} ${fasta} ${reads[0]} ${reads[1]} | \
    samtools view -S -b -@ ${samtools_threads} - | \
    samtools sort -@ ${samtools_threads} -o ${sample_id}_aligned.bam -
    
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
    bin3C mkmap -b ${bam} -o bin3C_out/map.h5
    bin3C cluster -m bin3C_out/map.h5 -o bin3C_out
    """
}

// 3. Juicer Tools Processing (修復：在 input 明確拆出 sample_id 變數)
process RUN_JUICER_TOOLS {
    tag "Sample: ${sample_id}"
    container "${params.sif1}"
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename -> "${sample_id}/juicer/${filename}" }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path chrom_sizes

    output:
    path "${sample_id}.hic", emit: hic_matrix

    script:
    """
	echo "=============================================="
	echo "Juicer Tools"
	echo "Sample ID : ${sample_id}"
	echo "BAM : ${bam}"
	echo "=============================================="
	
    samtools view -h -@ ${task.cpus} ${bam} | \
    awk -f /opt/HiC-Pro/bin/utils/hicpro2juicer.awk | \
    sort -k2,2d -k6,6d -S 4G --parallel=${task.cpus} > merged_nodups.txt

    java -Xmx${task.memory.toGiga()}g -jar /opt/juicer/juicer_tools.jar pre \
        merged_nodups.txt \
        ${sample_id}.hic \
        ${chrom_sizes}
    """
}
/*
========================================================================================
    WORKFLOW DEFINITION
========================================================================================
*/

workflow {
   log.info """\
   ================================================================
    3 D   G E N O M E   P I P E L I N E (DSL2)
   ================================================================
    Fastq Directory  : ${params.fastq_dir}
    Reads Pattern    : ${params.reads}
    Reference FASTA  : ${params.fasta}
    Output Dir       : ${params.outdir}
    Singularity SIF  : ${params.sif}
   ================================================================
   """

    ch_reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
    //ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
    ch_chrom_sizes = Channel.fromPath(params.chrom_sizes, checkIfExists: true)

    // 比對 (有多組 FASTQ 時，Nextflow 會自動多工/平行處理每一組 sample)
    ALIGN_BWA(
        ch_reads,
        params.fasta
    )

    // 下游分析 (各自傳入對應樣本的 BAM)
    RUN_BIN3C(ALIGN_BWA.out.bam)

    RUN_JUICER_TOOLS(
        ALIGN_BWA.out.bam,
		ch_chrom_sizes
    )
}