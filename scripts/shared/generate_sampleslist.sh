#!/bin/bash

#   Writes into a file a sequence of lines with this format 
#   col1            col2                            col3
#   "[FRAGMENT_ID]  [MATE_1_FILE_ABSOLUTE_PATH]     [MATE_2_FILE_ABSOLUTE_PATH]"
#   It needs to be run before running run_fastqc_with_slum_array.sh in order to have a formatted file 


#   Parameters: 
#           FASTQ_FOLDER (string): path of the directory where all the fastq files are stored. Don't put the '/' in the end. 
#                               for the trimmed fastq file 
#                               /data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/results/fastp
#           SUFFIX1 (str): the ending part of the name and the extension
#                               for the trimmed fastq file of this project 
#                               _1_trimmed.fastq.gz
#           
#           RELATIVE_PATH_OUTPUT_FILE (str): the relative path (from the project) of the outputfile. 
#                                       As it is a relative do not put '/' at the beginning
#    Returns: 
#       Generate the output file describe at the RELATIVE_PATH_OUTPUT_FILE

#   The name of the output name is fixed as it was needed also for the script run_fastqc_with_slum_array.sh 
#   ./[scriptname] [FASTQ_FOLDER] [SUFFIX1] [RELATIVE_PATH_OUTPUT_FILE]
#   
#   For the trimmed fastq files  
#   ./generate_sampleslist.sh /data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/results/fastp _1_trimmed.fastq.gz _2_trimmed.fastq.gz results/fastqc_after_trimming/intermediate_results/sampleslist.tsv
#   For the blood raw samples 
#   ./generate_sampleslist.sh /data/courses/rnaseq_course/toxoplasma_de/reads_Blood _1.fastq.gz _2.fastq.gz results/fastqc/intermediate_results/sampleslist.tsv
 

# save the path of the directory into the variable FASTQ_FOLDER
FASTQ_FOLDER="$1"
# the ending part of the mate1's name 
SUFFIX1="$2"
# teh ending part of the mate2's name 
SUFFIX2="$3"
# the file to store the path of the fastq files
RELATIVE_PATH_OUTPUT_FILE="$4"
# 


WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# path of the output file that will contain all the samples and fastq files' paths 
OUTPUT_FILE="$WORKDIR/$RELATIVE_PATH_OUTPUT_FILE"


# It goes through all the files that has the extention ".fastq.gz" and the name ends with "_1"
# So it searches all the fastq fukes with the reads of the forward strand of the fragment. 
for FILE in $FASTQ_FOLDER/*$SUFFIX1
do 
    # for each FILE it removes its end, which correspond to the pattern "_1.fastq.gz" and
    #it saves the remaning part of the path into the variable PREFIX
    PREFIX="${FILE%${SUFFIX1}}"
    # with the command basename it extract the ID of the fragment and save into the variable SAMPLE 
    SAMPLE=`basename $PREFIX`
    # prints to stdout the string with this format: "[FRAGMENT_ID]  [MATE_1_FILE_ABSOLUTE_PATH]  [MATE_2_FILE_ABSOLUTE_PATH]"
    # the fields are separated by '\t' 
    echo -e "${SAMPLE}\t$FILE\t${PREFIX}${SUFFIX2}" 
done > "$OUTPUT_FILE"
