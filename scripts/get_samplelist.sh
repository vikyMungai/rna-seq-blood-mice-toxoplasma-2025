#!/bin/bash

#   Writes to standard output a sequence of lines with this format 
#   col1            col2                            col3
#   "[FRAGMENT_ID]  [MATE_1_FILE_ABSOLUTE_PATH]     [MATE_2_FILE_ABSOLUTE_PATH]"
#   It needs to be run before running run_fastqc_with_slum_array.sh in order to have a formatted file 

#   It is suggested to redirect the output to a file in order to make the list of sample names and files
#   available to the new script to run fastqc (the name of the output name is fixed): 
#   ./[scriptname] [FASTQ_FOLDER] > samplelist.tsv


#   Parameters: 
#        FASTQ_FOLDER(string): path of the directory where all the fastq files are stored

#    Returns: 
#       None: write to stdout the formatted lines 
        

# save the path of the directory into the variable FASTQ_FOLDER
FASTQ_FOLDER=$1

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
done
