= Quality checks 

For the quality check it has been decided to use the tool fastqc. 
It has been decided to use a container. 

To have more efficiency I should have used parallel job with slurm job arrays, but I didn't mangaed to do so. 
Nevertheless, it this case it would have been usefult to count the files I have in order to know the length of my array. 

This is the algorith that I could have used. 
```sh
DIR="/data/courses/rnaseq_course/toxoplasma_de/reads_Blood"
FILES=`ls $DIR/*.fastq.gz`

COUNT=0
for FILE in $FILES; do 
echo $FILE; 
((COUNT++)) 
done

echo $COUNT
```

=== scripts for this step 
I decided 

When launching the fastqc command I specify where to write as the directory where I am reading the fastq files is read-only. I am using this argument:  
```sh 
    -o --outdir     Create all output files in the specified output directory.
                    Please note that this directory must exist as the program
                    will not create it.  If this option is not set then the 
                    output file for each sequence file is created in the same
                    directory as the sequence file which was processed.
``` 

=== Note for me 
This is an example of a fastqc ended in the right way. 
In the error file: 
```
Started analysis of SRR7821949_1.fastq.gz
Approx 5% complete for SRR7821949_1.fastq.gz
Approx 10% complete for SRR7821949_1.fastq.gz
Approx 15% complete for SRR7821949_1.fastq.gz
Approx 20% complete for SRR7821949_1.fastq.gz
Approx 25% complete for SRR7821949_1.fastq.gz
Approx 30% complete for SRR7821949_1.fastq.gz
Approx 35% complete for SRR7821949_1.fastq.gz
Approx 40% complete for SRR7821949_1.fastq.gz
Approx 45% complete for SRR7821949_1.fastq.gz
Approx 50% complete for SRR7821949_1.fastq.gz
Approx 55% complete for SRR7821949_1.fastq.gz
Approx 60% complete for SRR7821949_1.fastq.gz
Approx 65% complete for SRR7821949_1.fastq.gz
Approx 70% complete for SRR7821949_1.fastq.gz
Approx 75% complete for SRR7821949_1.fastq.gz
Approx 80% complete for SRR7821949_1.fastq.gz
Approx 85% complete for SRR7821949_1.fastq.gz
Approx 90% complete for SRR7821949_1.fastq.gz
Approx 95% complete for SRR7821949_1.fastq.gz
```
in the output file 
```
application/gzip
Analysis complete for SRR7821949_2.fastq.gz
```