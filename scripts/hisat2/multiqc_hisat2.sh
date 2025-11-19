#!/bin/bash


#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=multiqc_hisat2
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/hisat2/multiqc_hisat2/output_multiqc_%J.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/hisat2/multiqc_hisat2/error_multiqc_%J.err


# path of the container multiqc
CONTAINER="/containers/apptainer/multiqc-1.19.sif"

# directory of the project (this is an absolute path, so depending on which machine it is run it should be changed) 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# directory of the mapping summary 
SUMMARY_MAP_DIR=`ls $WORKDIR/results/hisat2/map_reads_summary_convert_sort_bam/alignment_summary/*`


# output directory for the sorted bam files 
OUTPUT_DIR="$WORKDIR/results/hisat2/multiqc"


# only if the directory does not exist it will be created 
if [ ! -d $OUTPUT_DIR ]; then 
    # option -p create the parents' folders if they do not exist
    mkdir -p $OUTPUT_DIR
fi 

# run multiqc quality check using the container 
apptainer exec \
    --bind /data \
    $CONTAINER multiqc -f -o "$OUTPUT_DIR" $SUMMARY_MAP_DIR

# the option -f was added to force the overwrting of the report to simplify the re-execution of the command 
# as in case the report already exists the multiqc command will give an error 