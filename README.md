# RNA sequencing pipeline 
project name : rna-seq-blood-mice-toxoplasma-2025

## Author 
Vittoria Mungai 

## Introduction  
This repository contains a bioinformatics pipeline designed for analysing fastq files of blood samples of mice infected with Toxoplasma gondii or not obtained from Gene Expression Omnibus (GEO), accession GSE119855, and they are a subset of Singhania et al. (2019). 
The pipeline is capable of performing quality control, mapping to the reference genome, counting the reads per gene, differetial analysis and overrepresentation analysis. 

## Structure of the project 
```
rna-seq-blood-mice-toxoplasma-2025/
├── scripts/
│   ├── data_analysis/
│   │   ├── 2_factor_data_analysis.R
│   │   └── des_eq_custom_functions.R
│   │
│   ├── fastp_trimming/
│   │   └── trim_adapters_fastp.sh
│   │
│   ├── fastqc/
│   │   └── run_fastqc_with_slurm_array.sh
│   │
│   ├── feature_counts/
│   │   ├── feature_counts.sh
│   │   ├── format_feature_counts_file.sh
│   │   └── multiqc_feature_counts.sh
│   │
│   ├── hisat2/
│   │   ├── index_sorted_bam_files.sh
│   │   ├── map_reads_summary_convert_sam.sh
│   │   ├── multiqc_hisat2.sh
│   │   └── run_hisat2_index_files.sh
│   │
│   ├── multiqc/
│   │   ├── config_multiqc_final_report.yaml
│   │   ├── run_multiqc.sh
│   │   └── run_multiqc_final_report.sh
│   │
│   └── shared/
│       ├── generate_sampleslist.sh
│       └── generate_sorted_bam_list.sh
├── .gitignore
└── README.md
```

## Sample data 
The fastq files from the blood samples are stored in the directory: `/data/courses/rnaseq_course/toxoplasma_de/reads_Blood/`. It's possible to access the folder directly when needed or create a logical link with the command `ln -s [name_link] [name_original_file]`

## Workflow 
Starting from the fastq files of the samples the workflow is divided in 7 main steps: 
1. Accessing the raw data from the cluster 
2. Quality checks
3. Map reads to the reference genome
4. Count the number of reads per gene
5. Exploratory data analysis
6. Differential expression analysis
7. Overrepresentation analysis

Usage notes: 
- Change the path of the working directory `WORK_DIR` in the `.sh`
- In the `R` scripts change the variabile `r_scripts_dir` with the path where the script `des_eq_custom_functions.R`is downloaded and the variable `work_dir` with the directory where the results should be saved. 
- In scripts where parameters are requested there is the documentation that explain what is needed. 
Both in this README (in the [pipeline steps](#pipeline-steps)) and in the documentation inside the script there is the specific command with the argument for the specific cases. 
- If a scripts access a file which needs to be generated, a explicit comment in the script will tell you. Nevertheless, in the step of the workflow it would be clearly state when and how to generate these files. 
- The output of the script is generated in the folder ./results. Although, further details are given below. 

## Pipeline steps 
Below is given the workflow in detail and in which order the scripts has to be run, they are all in the directory `./scripts`. If not other information is given, to run the script use `sbatch [filename]`. 

NB: The scripts that don't use sbatch must be run on an interactive node and not in the login node. 

Before submitting a batch script create the logs folders.  
```
logs/
├── errors/
│   ├── fastp/
│   ├── fastqc/
│   ├── feature_counts/
│   ├── feature_counts_not_zipped_gtf/
│   ├── hisat2/
│   └── multiqc/
│
└── outputs/
    ├── fastp/
    ├── fastqc/
    ├── feature_counts/
    ├── feature_counts_not_zipped_gtf/
    ├── hisat2/
    └── multiqc/
```

### 1. Accessing the raw data from the cluster 
The fastq files for the Blood sample are stored in the directory `/data/courses/rnaseq_course/toxoplasma_de/reads_Blood/`. 

### 2. Quality checks 
1. Check the quality of the raw reads with fastqc 
    - `shared/generate_sampleslist.sh` to generate the file `results/fastqc/intermediate_results/sampleslist.tsv` with the list of all the samples and the corrispective fastq files. Before create the folder `results/fastqc/intermediate_results/sampleslist.tsv`
        ``` bash 
        ./generate_sampleslist.sh /data/courses/rnaseq_course/toxoplasma_de/reads_Blood _1.fastq.gz _2.fastq.gz results/fastqc/intermediate_results/sampleslist.tsv
        ```
    - `fastqc/run_fastqc_with_slum_array.sh` to generate the fastqc files for each fragment. 
        ``` bash 
        sbatch run_fastqc_with_slum_array.sh results/fastqc/intermediate_results/sampleslist.tsv results/fastqc/
        ```
2. Check the quality of the raw reads with multiqc 
    - `multiqc/run_multiqc.sh` to generate the multiqc report and have a overview of all the fragments in one single report. 
        ``` bash 
        sbatch run_multiqc.sh results/fastqc results/multiqc
        ```
3. Trimming the raw reads to have a better quality and remove the adapters 
    - `fastp_trimming/trim_adapters_fastp.sh` to trim the fastq files. 
    This script uses the file `results/fastqc/intermediate_results/sampleslist.tsv` generated by the script `/scripts/fastqc/generate_sampleslist.sh` in the previous steps. 
    - `shared/generate_sampleslist.sh` with the arguments specified in the script to generate the list of trimmed files used to run fastqc with the trimmed reads. Create the folder `results/fastqc_after_trimming/intermediate_results` before running the script.
        ```
        ./generate_sampleslist.sh /data/users/vmungai/rna_seq_projects/rna-seq-blood-mice-toxoplasma-2025/results/fastp _1_trimmed.fastq.gz _2_trimmed.fastq.gz results/fastqc_after_trimming/intermediate_results/sampleslist.tsv
        ```
    It will generate the file `results/fastqc_after_trimming/intermediate_results/sampleslist.tsv`

4. Check the quality of the trimmed reads with fastqc 
    - `fastqc/run_fastqc_with_slum_array.sh` to generate the fastqc files for each trimmed sample 
        ```
        sbatch run_fastqc_with_slum_array.sh results/fastqc_after_trimming/intermediate_results/sampleslist.tsv results/fastqc_after_trimming/
        ```
5. Check the quality of the trimmed reads with multiqc 
    - `multiqc/run_multiqc.sh` to generate the multiqc report and have a overview of all the fragments in one single report. 
        ```
        sbatch run_multiqc.sh results/fastqc_after_trimming results/multiqc_after_trimming
        ```
 
#### 3. Map reads to the reference genome
1.  Produce all required index files for Hisat2
    - `hisat2/run_hisat2_index_files.sh`
2. For each sample separately, map the reads to the reference genome using Hisat2. The correct strandedness setting for this library prep protocol is RF.
3. Convert the resulting sam files to bam format using Samtools 
4. Sort the bam files by genomic coordinates using Samtools
    - `map_reads_summary_convert_sort_bam.sh` the steps 2, 3, and 4 are achieved within this script (for efficiency in storage). The output will be only the sorted bam file. 
5. Index the coordinate sorted bam files using Samtools 
    - `shared/generate_sorted_bam_list.sh` to generate the list of the sorted bam files that will be used to use the slurm array in the script: `index_sorted_bam_files.sh`
        ```
        ./generate_sorted_bam_list
        ```
    - `hisat2/index_sorted_bam_files.sh` it uses the file "sorted_bam_list.tsv" to index the bam files 
6. Check the quality of the mapping through multiqc
    - `hisat2/multiqc_hisat2.sh`to generate a .html report of the quality of the mapping 

#### 4. Count the number of reads per gene
1. Count the number of reads per gene usign featureCounts
    With `featureCounts` a table of counts is generated containing the number of reads per gene in each sample. 
    - `feature_counts/feature_counts.sh` will generate the table of counts and the counting summary. 
2. Analyse the output of featureCounts with multiqc 
    - `feature_counts/multiqc_feature_counts.sh` will generate a .html report of featureCounts. 


To generate a final report of all the previous steps run: 
    - `multiqc/run_multiqc_final_report.sh` you have to pass all the directories 


#### 5-7. Exploratory data analysis, Differential expression analysis and Overrepresentation analysis
1. To remove the first line and the columns containing Chr, Start, End, Strand and Length
    - `feature_counts/format_feature_counts_file.sh` to reformat the table in `feature_count.txt` to correspond to the format expected by DESeq2.
2. To run the differrential expression analysis and the overrepresentation analysis: 
    - download locally the scripts `data_analysis/2_factor_data_analysis.R` and `data_analysis/des_eq_custom_functions.R`, they have to be in teh same directory because the first script is usign the function from the second one. 
    - execute the script `data_analysis/2_factor_data_analysis.R`



## Logs folder 
The `logs` folder generated from running the batch script is structured like this: 
```
logs/
├── errors/
│   ├── fastp/
│   ├── fastqc/
│   ├── feature_counts/
│   ├── feature_counts_not_zipped_gtf/
│   ├── hisat2/
│   └── multiqc/
│
└── outputs/
    ├── fastp/
    ├── fastqc/
    ├── feature_counts/
    ├── feature_counts_not_zipped_gtf/
    ├── hisat2/
    └── multiqc/
```

## The output folder 
The `results` folder, which is generated from running the batch scripts on the cluster, is structure like this: 
```
results/
├── fastp/
│   ├── intermediate_results
├── fastqc/
│   ├── intermediate_results
├── fastqc_after_trimming/
│   └── intermediate_results
│
├── feature_counts/
│   ├── formatted_files/
│   ├── multiqc/
│   ├── feature_count.txt
│   └── feature_count.txt.summary
│
├── hisat2/
│   ├── index_files/
│   ├── index_sorted_bam_files/
│   ├── map_reads_summary_convert_sort_bam_files/
│   └── multiqc/
│
├── multiqc/
│   ├── multiqc_data/
│   └── multiqc_report.html
│
├── multiqc_after_trimming/
│   ├── multiqc_data/
│   └── multiqc_report.html
│
└── multiqc_final_report/
    ├── multiqc_data/
    └── multiqc_report.html
```

The R scripts will generate locally this directories for the output: 
```
results/
└──2_factors/
    ├── 5_exploratory_data_analysis/
    ├── 6_differential_analysis/
    │   ├── LFC_1_padj_0-01/
    │   └── LFC_1_padj_0-05/
    └── 7_overrepresentation_analysis/
        ├── LFC_1_padj_0-01/
        └── LFC_1_padj_0-05/
```

The **directories are created automatically**, unless otherwise specified.

## Dependencies
This pipeline uses containers thorugh apptainer. To work all the `.sh` scripts have to be run on the IBU cluster. 