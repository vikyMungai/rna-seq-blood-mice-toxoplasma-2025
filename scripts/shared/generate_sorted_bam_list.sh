#
# It creates a file with a list of all the sorted file in the specified directory 
# in each line you have [SORTED_BAM_PATH] 

# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# directory of sorted BAM files
DIR_BAM_FILE="$WORKDIR/results/hisat2/map_reads_summary_convert_sort_bam"
# dirtectory output file 
OUTPUT_DIR="$WORKDIR/results/hisat2/index_sorted_bam_files/intermediate_results"

# output file's path 
OUTPUT_FILE="$OUTPUT_DIR/sorted_bam_list.tsv"

# only if the directory does not exist it will be created 
if [ ! -d $OUTPUT_DIR ]; then 
    # option -p create the parents' folders if they do not exist
    mkdir -p $OUTPUT_DIR
fi 

# list all the .bam file 
BAM_FILES=`ls $DIR_BAM_FILE/*_sorted.bam`


for FILE in $BAM_FILES; do 
    
    echo "$FILE" 

done > "$OUTPUT_FILE"