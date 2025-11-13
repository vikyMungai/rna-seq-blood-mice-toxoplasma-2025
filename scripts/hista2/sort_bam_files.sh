#!/bin/bash

#SBATCH --array=1-15
#SBATCH --time=03:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sort_bam_files
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/hista2/sort_bam_files/output_sort_bam_files_%A_%a.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/hista2/sort_bam_files/error_sort_bam_files_%A_%a.err

# This file uses the file "sam_bam_files_list.tsv" generated with the script "generate_sam_bam_lists.sh"


# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# file where both path of .sam and future .bam files are stored 
SAM_BAM_FILES_LIST="$WORKDIR/results/hista2/convert_sam_to_bam/intermediate_results/sam_bam_files_list.tsv"
# output directory for the sorted bam files 
OUTPUT_DIR="$WORKDIR/results/hista2/sort_bam_files"

# read the path of the future .bam file from the 2nd column SAM_BAM_FILES_LIST
BAM_FILE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAM_BAM_FILES_LIST`

# I create the absolute path of the sorted bam file
BASE_NAME=`basename $BAM_FILE`
SORTED_BAM_FILE="$OUTPUT_DIR/${BASE_NAME%*.bam}_sorted.bam"



apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools sort -m $SLURM_MEM_PER_CPU -@ $SLURM_CPUS_PER_TASK -o "$SORTED_BAM_FILE" -T temp "$BAM_FILE"