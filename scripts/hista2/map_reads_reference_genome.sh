#!/bin/bash


#SBATCH --array=1-15
#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=hista2_mapping
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/hista2/mapping/output_hista2_map_%A_%a.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/hista2/mapping/error_hista2_map_%A_%a.err


# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# index basename to name all the index file 
DIRECTORY_INDEX_FILES="$WORKDIR/results/hista2/index_files"
INDEX_BASENAME="$DIRECTORY_INDEX_FILES/index_file"


# file with the list of the fragments. See script `generate_sampleslist.sh` for further details
SAMPLELIST="$WORKDIR/results/fastqc/intermediate_results/sampleslist.tsv"

# SLURM_ARRAY_TASK_ID we use the ID of the array to access the right line in the file samplelist.tsv

# take the first col with awk 
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
# take the second col with awk 
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
# take the third col with awk 
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# File for SAM output
OUTPUT_SAM_FILE="$WORKDIR/results/hista2/map_reads/${SAMPLE}_mapped_read.sam"

# mapping the reads mate1 and mate2 with the slurm array's indexes 
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2 -x "$INDEX_BASENAME" -1 $READ1 -2 $READ2 -S $OUTPUT_SAM_FILE -p $SLURM_CPUS_PER_TASK