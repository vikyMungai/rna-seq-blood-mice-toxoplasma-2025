#!/bin/bash

# range of the array we want to use (they are used to have the right input)
#SBATCH --array=1-12
#SBATCH --time=00:10:00
#SBATCH --mem=1g
#SBATCH --cpus-per-task=1
#SBATCH --job-name=slurm_array
#SBATCH --output=array_%J.out # I would put also the %j 
#SBATCH --error=array_%J.err # I would put also the %j 
#SBATCH --partition=pibu_el8

CONTAINER="/containers/apptainer/fastqc-0.12.1.sif"

# define variables
# my data user directory 
WORKDIR="/data/users/<username>/rnaseq_course"
# my result directory 
OUTDIR="$WORKDIR/results"
# I create a txt file that creates all the sample data 
SAMPLELIST="$WORKDIR/metadata/samplelist.tsv"


# SLURM_ARRAY_TASK_ID we use the ID of the array to access the right line in the file 

# take the first col with awk 
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
# take the second col with awk 
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
# take the third col with awk 
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

OUTFILE="$OUTDIR/${SAMPLE}.txt"

############################


mkdir -p $OUTDIR

apptainer exec --bind /data $CONTAINER fastqc --threads $SLURM_CPUS_PER_TASK -o "$OUTPUT_DIR" "$READ1" "$READ2"


sleep 30



