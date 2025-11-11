#!/bin/bash

#   Writes into a file a sequence of lines with this format 
#   col1            col2                            col3
#   "[FRAGMENT_ID]  [MATE_1_FILE_ABSOLUTE_PATH]     [MATE_2_FILE_ABSOLUTE_PATH]"
#   It needs to be run before running run_fastqc_with_slum_array.sh in order to have a formatted file 

#   The name of the output name is fixed as it was needed also for the script run_fastqc_with_slum_array.sh 
#   ./[scriptname] [FASTQ_FOLDER] 
#   For the blood sample 
#   ./generate_sampleslist.sh 


#   Parameters: 
#        FASTQ_FOLDER(string): path of the directory where all the fastq files are stored. Don't put the '/' in the end. 
#                               for the blood samples of this project: 
#                               /data/courses/rnaseq_course/toxoplasma_de/reads_Blood

#    Returns: 
#       Generate the file sampleslist.tsv
        

# save the path of the directory into the variable FASTQ_FOLDER
FASTQ_FOLDER="$1"
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# path of the output file that will contain all the samples and fastq files' paths 
OUTPUT_FILE="$WORKDIR/results/fastqc/intermediate_results/sampleslist.tsv"

# It goes through all the files that has the extention ".fastq.gz" and the name ends with "_1"
# So it searches all the fastq fukes with the reads of the forward strand of the fragment. 
for FILE in $FASTQ_FOLDER/*_*1.fastq.gz
do 
    # for each FILE it removes its end, which correspond to the pattern "_1.fastq.gz" and
    #it saves the remaning part of the path into the variable PREFIX
    PREFIX="${FILE%_*.fastq.gz}"
    # with the command basename it extract the ID of the fragment and save into the variable SAMPLE 
    SAMPLE=`basename $PREFIX`
    # prints to stdout the string with this format: "[FRAGMENT_ID]  [MATE_1_FILE_ABSOLUTE_PATH]  [MATE_2_FILE_ABSOLUTE_PATH]"
    # the fields are separated by '\t' 
    echo -e "${SAMPLE}\t$FILE\t${FILE%?.fastq.gz}2.fastq.gz" 
done > "$OUTPUT_FILE"
