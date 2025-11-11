#!/usr/bin/env bash 

# This is the directory where the blood samples are 
DIR="/data/courses/rnaseq_course/toxoplasma_de/reads_Blood"
# I list all the fastq.gz files that are into teh directory and save them as a list 
# in order to go through them with a loop
FILES=`ls $DIR/*.fastq.gz`

# The script that will be executed needs to be in the same directory 
FASTQC_SINGLE_SAMPLE_SCRIPT="run_fastqc_single_sample.sh"

# for each loop I execute the script FASTQC_SINGLE_SAMPLE_SCRIPT passing the filename 
# of one of the fastq files 
for FILE in $FILES; do 
    sbatch --export=FILE=$FILE "$FASTQC_SINGLE_SAMPLE_SCRIPT"
done 