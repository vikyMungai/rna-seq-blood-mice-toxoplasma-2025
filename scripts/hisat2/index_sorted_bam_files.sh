#!/bin/bash

#SBATCH --array=1-15
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=index_sorted_bam_files
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/hisat2/index_sorted_bam_files/output_index_sorted_bam_files_%A_%a.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/hisat2/index_sorted_bam_files/error_index_sorted_bam_files_%A_%a.err


# This file uses the file "sorted_bam_list.tsv" generated with the script "generate_sorted_bam_list.sh"

# container path for hisat samtools 
CONTAINER_HISAT_SAMTOOLS="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# file where sorted .bam files are stored 
SORTED_BAM_FILES_LIST="$WORKDIR/results/hisat2/index_sorted_bam_files/intermediate_results/sorted_bam_list.tsv"

# output directory for the sorted bam files 
OUTPUT_DIR="$WORKDIR/results/hisat2/index_sorted_bam_files"



# read the path of the .bam file from the 1st column SORTED_BAM_FILES_LIST
SORTED_BAM_FILE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SORTED_BAM_FILES_LIST`

# I create the absolute path of the indexed sorted bam file
BASE_NAME=`basename $SORTED_BAM_FILE`
INDEXED_S_BAM_FILE="$OUTPUT_DIR/${BASE_NAME%*.bam}.index"

# for debugging 
echo "DB: SORTED_BAM_FILE: $SORTED_BAM_FILE, INDEXED_S_BAM_FILE: $INDEXED_S_BAM_FILE, BASE_NAME: $BASE_NAME"

apptainer exec --bind /data/ $CONTAINER_HISAT_SAMTOOLS samtools index "$SORTED_BAM_FILE" -o "$INDEXED_S_BAM_FILE"