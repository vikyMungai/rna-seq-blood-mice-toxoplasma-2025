#!/bin/bash

# the number of arrays matches the number of fragments 
#SBATCH --array=1-15
#SBATCH --time=02:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=run_fastqc_with_slum_array
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/fastqc/output_fastqc_%A_%a.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/fastqc/error_fastqc_%A_%a.err

# path of the container fastqc
CONTAINER="/containers/apptainer/fastqc-0.12.1.sif"

# directory of the project (this is an absolute path, so depending on which machine it is run it should be changed) 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# my result directory 
OUTDIR="$WORKDIR/results/fastqc/"
# file with the list of the fragments. See script `generate_sampleslist.sh` for further details
SAMPLELIST="$WORKDIR/results/fastqc/intermediate_results/sampleslist.tsv"


# SLURM_ARRAY_TASK_ID we use the ID of the array to access the right line in the file 

# take the first col with awk 
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
# take the second col with awk 
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
# take the third col with awk 
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# for debugging 
# echo "DB: read1: $READ1, read2: $READ2, sampleslist: $SAMPLELIST"


# only if the directory does not exist it will be created 
if [ ! -d $OUTDIR ]; then 
    # option -p create the parents' folders if they do not exist
    mkdir -p $OUTDIR
fi 

# exectuing the quality check with fastqc using the container, two fastq files are checked at the same time 
apptainer exec --bind /data $CONTAINER fastqc --threads $SLURM_CPUS_PER_TASK -o "$OUTDIR" "$READ1" "$READ2"


sleep 30



