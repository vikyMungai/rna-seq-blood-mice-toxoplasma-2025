#
# It creates a file with a list of all the sam file in the specified directory 
# and respective bam file associated 
# [SAM_FILE_PATH][BAM_FILE_PATH] 

# To use the sam file generated with the script 
# ./generate_sam_bam_lists.sh results/hista2/convert_sam_to_bam/intermediate_results/sam_bam_files_list.tsv 

#   Parameters: 
#           
#           RELATIVE_PATH_OUTPUT_FILE (str): the relative path (from the project) of the outputfile. 
#                                       As it is a relative do not put '/' at the beginning
#    Returns: 
#       Generate the output file describe at the RELATIVE_PATH_OUTPUT_FILE


# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# directory of SAMs files
DIR_SAM_FILE="$WORKDIR/results/hista2/map_reads"

# output file's path 
OUTPUT_FILE="$WORKDIR/$1"

# list all the .sam file 
SAM_FILES=`ls $DIR_SAM_FILE/*_mapped_read.sam`
# directory where the .bam file will be stored
BAM_FILES_DIR="$WORKDIR/results/hista2/convert_sam_to_bam"



for FILE in $SAM_FILES; do 
    BASE_NAME=`basename $FILE`
    
    echo "$FILE" "$BAM_FILES_DIR/${BASE_NAME%*.sam}.bam" 

done > "$OUTPUT_FILE"