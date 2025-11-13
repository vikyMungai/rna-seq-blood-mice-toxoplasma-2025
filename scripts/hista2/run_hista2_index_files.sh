#!/bin/bash


#SBATCH --time=03:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=hista2_index_files
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/hista2/output_hista2_index_files_%J.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/hista2/error_hista2_index_files_%J.err


# working directory of the project 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# path where the reference sequence fasta is 
REFERENCE_SEQUENCE_FASTA="$WORKDIR/data/processed_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
# index basename to name all the index file 
DIRECTORY_INDEX_FILES="$WORKDIR/results/hista2/index_files"
INDEX_BASENAME="$DIRECTORY_INDEX_FILES/index_file"

apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2-build $REFERENCE_SEQUENCE_FASTA $INDEX_BASENAME