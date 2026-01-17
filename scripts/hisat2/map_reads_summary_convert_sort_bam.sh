#!/bin/bash


#SBATCH --array=1-15
#SBATCH --time=08:00:00
#SBATCH --mem=55G
#SBATCH --cpus-per-task=5
#SBATCH --job-name=hisat2_mapping_summary
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/hisat2/map_reads_summary_convert_sort_bam/output_map_%A_%a.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/hisat2/map_reads_summary_convert_sort_bam/error_map_%A_%a.err

# this files uses the file SAMPLELIST to read the path of the fastq files. So you have to check that the file exists.
# In case it is not present, run the script '/scripts/shared/generate_sampleslist.sh' to generate it 

# container path for hisat samtools 
CONTAINER_HISAT_SAMTOOLS="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# index basename to name all the index file 
DIRECTORY_INDEX_FILES="$WORKDIR/results/hisat2/index_files"
INDEX_BASENAME="$DIRECTORY_INDEX_FILES/index_file"

# only if the directory does not exist it will be created 
if [ ! -d $DIRECTORY_INDEX_FILES ]; then 
    # option -p create the parents' folders if they do not exist
    mkdir -p $DIRECTORY_INDEX_FILES
fi 


# file with the list of the fragments. See script `generate_sampleslist.sh` for further details
SAMPLELIST="$WORKDIR/results/fastqc_after_trimming/intermediate_results/sampleslist.tsv"

# SLURM_ARRAY_TASK_ID we use the ID of the array to access the right line in the file samplelist.tsv
# each slurm job will run all the command for each sample 

# take the first col with awk 
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
# take the second col with awk 
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
# take the third col with awk 
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`


# File for SAM output
# we use the SCRATCH directory to not storage the sam files on the cluster given their size 
SAM_FILE="$SCRATCH/${SAMPLE}_mapped_read.sam"
#directory for the sorted bam files 
BAM_FILES_DIR="$WORKDIR/results/hisat2/map_reads_summary_convert_sort_bam"
# directory for the alignment summary from the hisat2 mapping command 
SUMMARY_DIR="$WORKDIR/results/hisat2/map_reads_summary_convert_sort_bam/alignment_summary"
SUMMARY_FILE="$SUMMARY_DIR/${SAMPLE}_alignment.summary"

SORTED_BAM_FILE="$BAM_FILES_DIR/${SAMPLE}_sorted.bam"

# only if the directory does not exist it will be created 
if [ ! -d $BAM_FILES_DIR ]; then 
    # option -p create the parents' folders if they do not exist
    mkdir -p $BAM_FILES_DIR
fi 

# only if the directory does not exist it will be created 
if [ ! -d $SUMMARY_DIR ]; then 
    # option -p create the parents' folders if they do not exist
    mkdir -p $SUMMARY_DIR
fi 

# mapping the reads mate1 and mate2 of the sample with hisat2 
apptainer exec --bind /data/ \
     $CONTAINER_HISAT_SAMTOOLS hisat2 --new-summary --summary-file "$SUMMARY_FILE" -x "$INDEX_BASENAME" -1 $READ1 -2 $READ2 -S $SAM_FILE -p $SLURM_CPUS_PER_TASK
# the options --new-summary --summary-file are used to print alignment summary in a new style, which is more machine-friendly.

# pipeline: convert sam file to bam file and sort the bam file 
apptainer exec --bind /data/ \
    $CONTAINER_HISAT_SAMTOOLS samtools view -hbS "$SAM_FILE" | 
    apptainer exec --bind /data/ \
    $CONTAINER_HISAT_SAMTOOLS samtools sort -m 30G -@ $SLURM_CPUS_PER_TASK \
    -o "$SORTED_BAM_FILE" -T temp
