#!/usr/bin/env bash 

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1000
#SBATCH --time=00:20:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_couple_sample
#SBATCH --array=1-2
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/output/output_fastqc_%A_SRR7821949_%a.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/error/error_fastqc_%A_SRR7821949_%a.err

DIR="/data/courses/rnaseq_course/toxoplasma_de/reads_Blood"

FILE="${DIR}/SRR7821950_${SLURM_ARRAY_TASK_ID}.fastq.gz"

OUTPUT_DIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/data/external_data/fastqc_output/"

# running the fastQC using a container 
apptainer exec --bind /data /containers/apptainer/fastqc-0.12.1.sif fastqc --threads $SLURM_CPUS_PER_TASK -o "$OUTPUT_DIR" "$FILE"

