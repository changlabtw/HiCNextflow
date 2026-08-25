#!/bin/sh
#SBATCH -J HIC_job                          # Job name
#SBATCH -p cputest                          # Partition Name 
#SBATCH -N 1                                # Maximum number of "nodes" to be allocated 
#SBATCH -n 1                                # Maximum number of 'tasks' to be allocated 
#SBATCH --cpus-per-task=96
#SBATCH --mem=180g                          # memory size
#SBATCH -t 72:00:00                         #runing time
#SBATCH -o hic_workflow_out.%j.log          # Path to the output file 
#SBATCH -e hic_workflow_err.%j.log          # Path to the error log



SOFT_DIR="/home/dhllove/Soft/GenomeAnalysis"
SOFT_IMG="/home/dhllove/DockerImage"

DATA_DIR="/home/dhllove/NGS_DATA"
WORK_DIR="/home/dhllove/Work"
REF_DIR="reference_FASTA"
REFG_FILE19="Homo_sapiens_assembly19.fasta"
REFG_FILE38="Homo_sapiens_assembly38.fasta"
VCF_GRCH38="Homo_sapiens_assembly38.known_indels.vcf"
GRCH_PATH=${DATA_DIR}/${REF_DIR}
VCF_PATH=${DATA_DIR}/reference_VCF


module load singularity/v3.8.7

export PATH="/home/dhllove/Soft/JAVA17/bin:$PATH"
export PATH="/home/dhllove/Soft/Nextflow:$PATH"

export NXF_FILE_LOCKS=false
export SINGULARITYENV_NXF_FILE_LOCKS=false
export APPTAINERENV_NXF_FILE_LOCKS=false

export NXF_HOME="/home/dhllove/work_dir/.nextflow"
export SINGULARITYENV_NXF_HOME="/home/dhllove/work_dir/.nextflow"
export APPTAINERENV_NXF_HOME="/home/dhllove/work_dir/.nextflow"

srun nextflow run hic_main.nf \
        -profile local ;


