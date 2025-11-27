#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=run_multiqc
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/multiqc/output_multiqc_%J.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/multiqc/error_multiqc_%J.err

#   Parameters: 
#           ABSOLUTE_PATH_INPUT_FILES (str): the relative path (from the project) of the files that multiqc has to process. 
#                                       As it is a relative do not put '/' at the beginning   
#           NB: it can be one or more paths  

# To run multiqc for fastqc before trimming, fastqc after trimming, hisat and featureCounts I passed each directory 
# because the results were divided in structured folders into 'results'


# sbatch run_multiqc_final_report.sh results/fastqc results/fastp results/fastqc_after_trimming results/hisat2/map_reads_summary_convert_sort_bam/alignment_summary results/feature_counts



# path of the container multiqc
CONTAINER="/containers/apptainer/multiqc-1.19.sif"

# directory of the project (this is an absolute path, so depending on which machine it is run it should be changed) 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# my result directory 
OUTDIR="$WORKDIR/results/multiqc_final_report/"


CONFIG_FILE="$WORKDIR/scripts/multiqc/config_multiqc_final_report.yaml"

# only if the directory does not exist it will be created 
if [ ! -d $OUTDIR ]; then 
    # option -p create the parents' folders if they do not exist
    mkdir -p $OUTDIR
fi 

# convert the relative path to absolute for all the arguments passed 
DIR_LIST=""

for REL_PATH in $@
do
        DIR_LIST="$DIR_LIST $WORKDIR/$REL_PATH"
done



# run multiqc quality check using the container 
apptainer exec \
    --bind /data \
    $CONTAINER multiqc -f --config "$CONFIG_FILE" -o "$OUTDIR" $DIR_LIST

# the option -f was added to force the overwrting of the report to simplify the re-execution of the command 
# as in case the report already exists the multiqc command will give an error 
 
