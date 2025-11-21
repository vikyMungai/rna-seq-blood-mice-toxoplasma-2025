#!/bin/bash


#SBATCH --time=03:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=feature_counts
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/feature_counts/output_feature_counts_%J.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/feature_counts/error_feature_counts_%J.err

# path of the container subread for featureCounts
CONTAINER_HISAT_SAMTOOLS="/containers/apptainer/subread_2.0.6.sif"


# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
ANNOTATION_FILE="$WORKDIR/data/processed_data/Mus_musculus.GRCm39.115.gtf.gz"
# directory where all the sorted bam files 
cd "$WORKDIR/results/hisat2/map_reads_convert_sort_bam"

# all the sorted and indexed bam files 
BAM_FILES=`ls ./*_sorted.bam`
# output directory for featureCounts 
OUTPUT_DIR="$WORKDIR/results/feature_counts"



OUTPUT_FILE="$OUTPUT_DIR/feature_count.txt"

# only if the directory does not exist it will be created 
if [ ! -d $OUTPUT_DIR ]; then 
    # option -p create the parents' folders if they do not exist
    mkdir -p $OUTPUT_DIR
fi 


# count the number of read per gene 
apptainer exec --bind /data/ $CONTAINER_HISAT_SAMTOOLS featureCounts -p -a $ANNOTATION_FILE -T $SLURM_CPUS_PER_TASK -s 2 -o $OUTPUT_FILE $BAM_FILES

# Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2]
# - The annotation file has the format .gtf, which is the default one for featureCounts. 
#   So,it is not specified with the -F option 
# - The option `-s 2` is set in order to perform a reversely stranded read counting. 
# - The option `-p` specify to featureCounts that we have paired-end strand. 

# -t exon 
# -g (by default )

