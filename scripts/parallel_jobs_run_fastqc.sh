#!/usr/bin/env bash 

DIR="/data/courses/rnaseq_course/toxoplasma_de/reads_Blood"
FILES=`ls $DIR/*.fastq.gz`


for FILE in $FILES; do 
    sbatch --export=FILE=$FILE /data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/scripts/run_fastqc_single_sample.sh
done 