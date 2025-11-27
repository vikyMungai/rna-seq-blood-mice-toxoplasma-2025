#!/bin/bash


#SBATCH --time=03:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=formatted_feature_counts
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/feature_counts/formatted_feature_counts/output_%J.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/feature_counts/formatted_feature_counts/error_%J.err


# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"

# output directory for featureCounts 
FEATURE_COUNT_DIR="$WORKDIR/results/feature_counts"
FEATURE_COUNT_FILE="$FEATURE_COUNT_DIR/feature_count.txt"

FORMATTED_FEATURE_COUNT_DIR="$WORKDIR/results/feature_counts/formatted_files"
FORMATTED_FEATURE_COUNT_FILE="$FORMATTED_FEATURE_COUNT_DIR/formatted_feature_count.txt"

# only if the directory does not exist it will be created 
if [ ! -d $FORMATTED_FEATURE_COUNT_DIR ]; then 
    # option -p create the parents' folders if they do not exist
    mkdir -p $FORMATTED_FEATURE_COUNT_DIR
fi 

echo $FEATURE_COUNT_FILE
ls -l $FEATURE_COUNT_FILE

# remove the first line and thw columns: Chr, Start, End, Strand and Length
tail -n +2 $FEATURE_COUNT_FILE |  cut -f1,7-  > $FORMATTED_FEATURE_COUNT_FILE