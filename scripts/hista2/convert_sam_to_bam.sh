#!/bin/bash

#SBATCH --array=1-15
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=convert_sam_to_bam
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/hista2/convert_sam_to_bam/output_convert_sam_to_bam_%A_%a.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/hista2/convert_sam_to_bam/error_convert_sam_to_bam_%A_%a.err

# This file uses the file "sam_bam_files_list.tsv" generated with the script "generate_sam_bam_lists.sh"


# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# file where both path of .sam and future .bam files are stored 
SAM_BAM_FILES_LIST="$WORKDIR/results/hista2/convert_sam_to_bam/intermediate_results/sam_bam_files_list.tsv"

# read the path of the sam file from 1st column the SAM_BAM_FILES_LIST
SAM_FILE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAM_BAM_FILES_LIST`
# read the path of the future .bam file from the 2nd column SAM_BAM_FILES_LIST
BAM_FILE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAM_BAM_FILES_LIST`

# with the slurm array it convert each .sam file into a .bam file 
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools view -hbS "$SAM_FILE" "$BAM_FILE"