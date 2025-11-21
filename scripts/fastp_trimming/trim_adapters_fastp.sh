#!/bin/bash

# the number of arrays matches the number of fragments 
#SBATCH --array=1-15
#SBATCH --time=02:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=trim_adapters_fastp
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=vittoria.mungai@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --output=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/outputs/fastp/output_fastp_%A_%a.out
#SBATCH --error=/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/logs/errors/fastp/error_fastp_%A_%a.err


# this files uses the file SAMPLELIST to read the path of the fastq files. So you have to check that the file exists.
# In case it is not present, run the script '/scripts/fastqc/generate_sampleslist.sh' to generate it 

# path of the container fastp
CONTAINER="/containers/apptainer/fastp_0.24.1.sif"

# directory of the project (this is an absolute path, so depending on which machine it is run it should be changed) 
WORKDIR="/data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025"
# my result directory 
OUTDIR="$WORKDIR/results/fastp"
# file with the list of the fragments. See script `generate_sampleslist.sh` for further details
SAMPLELIST="$WORKDIR/results/fastqc/intermediate_results/sampleslist.tsv"


# SLURM_ARRAY_TASK_ID we use the ID of the array to access the right line in the file 

# take the first col with awk 
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
# take the second col with awk 
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
# take the third col with awk 
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# only if the directory does not exist it will be created 
if [ ! -d $OUTDIR ]; then 
    # option -p create the parents' folders if they do not exist
    mkdir -p $OUTDIR
fi 

# trimming the adapters of the samples 
apptainer exec --bind /data \
                 $CONTAINER fastp \
                 --detect_adapter_for_pe \
                 --in1 "$READ1" --in2 "$READ2" \
                 --out1 "$OUTDIR/${SAMPLE}_1_trimmed.fastq.gz" --out2 "$OUTDIR/${SAMPLE}_2_trimmed.fastq.gz" 
# --in1 : fastq file for the Mate1 of the sample 
# --in2 : fastq file for the Mate2 of the sample 
# --out1 : file where you have the output of the fastp for the trimmed Mate1 
# --out2 : file where you have the output of the fastp for the trimmed Mate2


# add --detect_adapter_for_pe