#!/bin/bash


#SBATCH --time=03:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=multiqc_feature_counts
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/feature_counts/multiqc/output_feature_counts_%J.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/feature_counts/multiqc/error_feature_counts_%J.err


# path of the container multiqc
CONTAINER_MULTIQC="/containers/apptainer/multiqc-1.19.sif"

# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"

# directory for featureCounts' summary
FEATURE_COUNTS_DIR="$WORKDIR/results/feature_counts"
# output directory for multiqc report 
OUTPUT_DIR_MULTIQC="$WORKDIR/results/feature_counts/multiqc"

SUMMARY_FILE=`ls $FEATURE_COUNTS_DIR/*.summary`

# same for multiqc
if [ ! -d $OUTPUT_DIR_MULTIQC ]; then 
    mkdir -p $OUTPUT_DIR_MULTIQC
fi 

apptainer exec --bind /data/ $CONTAINER_MULTIQC multiqc -f -o $OUTPUT_DIR_MULTIQC $SUMMARY_FILE
