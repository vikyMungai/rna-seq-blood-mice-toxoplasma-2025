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
#           RELATIVE_PATH_OUTPUT_FILE (str): the relative path (from the project) of the outputfile. 
#                                       As it is a relative do not put '/' at the beginning            

# For the fastqc with raw blood samples
# sbatch run_multiqc.sh results/fastqc results/fastqc results/multiqc

# For the fastqc with trimmed samples
# sbatch run_multiqc.sh results/fastqc_after_trimming results/multiqc_after_trimming


# path of the container multiqc
CONTAINER="/containers/apptainer/multiqc-1.19.sif"

# directory of the project (this is an absolute path, so depending on which machine it is run it should be changed) 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# my result directory 
OUTDIR="$WORKDIR/$2/"
# path of the blood samples directory 
FASTQ_FOLDER="/data/courses/rnaseq_course/toxoplasma_de/reads_Blood"
# store the list of fastqc zips into a variable 
FASTQC_LIST=`ls $WORKDIR/$1/*fastqc.zip`


# only if the directory does not exist it will be created 
if [ ! -d $OUTDIR ]; then 
    # option -p create the parents' folders if they do not exist
    mkdir -p $OUTDIR
fi 

# run multiqc quality check using the container 
apptainer exec \
    --bind /data \
    $CONTAINER multiqc -f -o "$OUTDIR" $FASTQC_LIST

# the option -f was added to force the overwrting of the report to simplify the re-execution of the command 
# as in case the report already exists the multiqc command will give an error 
