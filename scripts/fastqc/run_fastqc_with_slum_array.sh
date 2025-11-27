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

#   Parameters: 
#           SAMPLELIST (str): the file .tsv with the list of all the fastq files  
#                                       As it is a relative do not put '/' at the beginning 
#           RELATIVE_PATH_OUTPUT_FILE (str): the relative path (from the project) of the outputfile. 
#                                       As it is a relative do not put '/' at the beginning            

# this files uses the file SAMPLELIST to read the path of the fastq files. So you have to check that the file exists.
# In case it is not present, run the script '/scripts/shared/generate_sampleslist.sh' to generate it 

# For the fastqc with raw blood samples
# sbatch run_fastqc_with_slum_array.sh results/fastqc/intermediate_results/sampleslist.tsv results/fastqc/

# For the fastqc with trimmed samples
# sbatch run_fastqc_with_slum_array.sh results/fastqc_after_trimming/intermediate_results/sampleslist.tsv results/fastqc_after_trimming/


# path of the container fastqc
CONTAINER="/containers/apptainer/fastqc-0.12.1.sif"

# directory of the project (this is an absolute path, so depending on which machine it is run it should be changed) 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# file with the list of the fragments. See script `generate_sampleslist.sh` for further details
SAMPLELIST="$WORKDIR/$1"
# my result directory 
OUTDIR="$WORKDIR/$2"



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



