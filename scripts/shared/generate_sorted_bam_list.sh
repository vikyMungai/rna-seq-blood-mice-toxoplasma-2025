#
# It creates a file with a list of all the sorted file in the specified directory 
# in each line you have [SORTED_BAM_PATH] 

# To use the sam file generated with the script 
# ./generate_sam_bam_lists.sh results/hista2/index_sorted_bam_files/intermediate_results/sorted_bam_list.tsv 

#   Parameters: 
#           
#           RELATIVE_PATH_OUTPUT_FILE (str): the relative path (from the project) of the outputfile. 
#                                       As it is a relative do not put '/' at the beginning
#    Returns: 
#       Generate the output file describe at the RELATIVE_PATH_OUTPUT_FILE


# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# directory of sorted BAM files
DIR_BAM_FILE="$WORKDIR/results/hista2/sort_bam_files"


# output file's path 
OUTPUT_FILE="$WORKDIR/$1"

# list all the .bam file 
BAM_FILES=`ls $DIR_BAM_FILE/*_sorted.bam`



for FILE in $BAM_FILES; do 
    
    echo "$FILE" 

done > "$OUTPUT_FILE"