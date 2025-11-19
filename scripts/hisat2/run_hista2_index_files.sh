#!/bin/bash


#SBATCH --time=03:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=hisat2_index_files
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/hisat2/index_files/output_index_files_%J.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/hisat2/index_files/error_index_files_%J.err

# container path for hisat samtools 
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# path where the reference sequence fasta is located
REFERENCE_SEQUENCE_FASTA="$WORKDIR/data/processed_data/Mus_musculus.GRCm39.dna.primary_assembly.fa"
# directory where all the index files are stored 
DIRECTORY_INDEX_FILES="$WORKDIR/results/hisat2/index_files"
# index basename to name all the index file 
INDEX_BASENAME="$DIRECTORY_INDEX_FILES/index_file"


# only if the directory does not exist it will be created 
if [ ! -d $DIRECTORY_INDEX_FILES ]; then 
    # option -p create the parents' folders if they do not exist
    mkdir -p $DIRECTORY_INDEX_FILES
fi 

apptainer exec --bind /data/ \
         $CONTAINER hisat2-build $REFERENCE_SEQUENCE_FASTA $INDEX_BASENAME