#!/bin/bash


#SBATCH --time=03:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=feature_counts
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/feature_counts/output_feature_counts_%J.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/feature_counts/error_feature_counts_%J.err

# container path for hisat samtools 
CONTAINER="/containers/apptainer/subread_2.0.6.sif"

# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
ANNOTATION_FILE="$WORKDIR/data/processed_data/Mus_musculus.GRCm39.115.gtf.gz"
# all the sorted and indexed bam files 
BAM_FILES=`ls $WORKDIR/results/hisat2/map_reads_convert_sort_bam/*_sorted.bam`
# output directory for featureCounts 
OUTPUT_DIR="$WORKDIR/results/feature_counts"

OUTPUT_FILE="$OUTPUT_DIR/feature_count.txt"


# Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2]

# count the number of read per gene 
apptainer exec --bind /data/ $CONTAINER featureCounts -p -a $ANNOTATION_FILE -T $SLURM_CPUS_PER_TASK -s 2 -o $OUTPUT_FILE $BAM_FILES

# - The annotation file has the format .gtf, which is the default one for featureCounts. 
#   So,it is not specified with the -F option 
# - The option `-s 2` is set in order to perform a reversely stranded read counting. 
# - The option `-p` specify to featureCounts that we have paired-end strand. 


