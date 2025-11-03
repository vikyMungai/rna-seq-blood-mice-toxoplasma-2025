#!/usr/bin/env bash 

# I know that I have 30 files so I need 30 CPU in total, as I need one CPU for each thread and I am creating on ethread for each file. 

#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=1000
#SBATCH --time=00:30:00
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastqc_single_sample
#SBATCH --array=0-29
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/output/output_fastqc_%A_SRR7821949_%a.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/error/error_fastqc_%A_SRR7821949_%a.err


DIR="/data/courses/rnaseq_course/toxoplasma_de/reads_Blood"
FILES=`ls $DIR/*.fastq.gz`

OUTPUT_DIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/data/external_data/fastqc_output/"

#echo "$FILES[$SLURM_ARRAY_TASK_ID]"
# using SLURM job arrays 
apptainer exec --bind /data /containers/apptainer/fastqc-0.12.1.sif fastqc --threads $SLURM_CPUS_PER_TASK -o "$OUTPUT_DIR" "$FILES"
