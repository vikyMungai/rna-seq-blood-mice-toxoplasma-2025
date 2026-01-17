# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# install.packages("ggplot2")
# BiocManager::install("DESeq2")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("clusterProfiler", force = TRUE)


library(DESeq2)

# libraries for the volcano plot 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

# library for the annotation of the genes names 
library(org.Mm.eg.db)

# library of overrepresentation analysis
library(clusterProfiler) # had a problem with the package 'enrichplot'

# library to plot enrichment results 
library(ggplot2)
source("des_eq_custom_functions.R")

### input files 
dir_input <- "/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/feature_counts"
### directories to store plots
dir_step5 <- "/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/2_factors/5_exploratory_data_analysis"
dir_step6 <- "/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/2_factors/6_differential_analysis"
dir_step7 <- "/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/2_factors/7_overrepresentation_analysis"
# LFC and padj thresholds 
LFC_limit = 1 
padj_limit = 0.01

#check if the dir_input exists 
if (!dir.exists(dir_input)) {
  stop(paste0("ERROR: The directory ", dir_input, " do not exist." ))
}

create_not_existing_dir(dir_step5)
create_not_existing_dir(dir_step6)
create_not_existing_dir(dir_step7)


#
# STEP 5: Exploratory data analysis
#

# Directory where the file "formatted_feature_count.txt" is 

setwd(dir_input)

# read counts data 
counts_data <- read.csv("formatted_feature_count.txt", 
                        header = TRUE, 
                        sep = "\t", 
                        row.names="Geneid") 
# row.names="Geneid" is used to have as rownames the gene id

# formatting the col names with the sample names 
new_col_names <- sub(".*(SRR[0-9]+).*", "\\1", colnames(counts_data))
colnames(counts_data) <- new_col_names

# read the sample information table 
setwd("/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/data/raw_data")
col_data <- read.csv("DESeq_colData_2_conditions.csv", 
                     header = TRUE, 
                     sep = ";")
# order the rows by sample 
col_data <- col_data[order(col_data$Sample), ]
# select the column "Sample" as rownames 
col_data <- data.frame(col_data, row.names = 1)

# check if the columns name of the sample in the counts_data matches 
# to the row names in colData
all(colnames(counts_data) %in% rownames(col_data))

# check if the columns name of the sample in the counts_data are in the same order 
# of the row names in colData
all(colnames(counts_data) == rownames(col_data))

# converting columns with the condition of col_data into factors 
col_data$Disease <- factor(col_data$Disease)
col_data$Genotype <- factor(col_data$Genotype)

# two condition separated with interaction 
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_data,
                              design = ~ Disease + Genotype + Genotype:Disease ) 

# by setting the factor level we explicit say which level we want as reference
# otherwise it will be chosen alphabetically 
# for the PCA plot is not important, but it is for the Differential Analysis 
dds$Disease <- relevel(dds$Disease, ref = "Control")
dds$Genotype <- relevel(dds$Genotype, ref = "WT")

# run DESeq
dds <- DESeq(dds)

# the conditions 
resultsNames(dds)

# it is important to not have a correlation between the mean expression and the variance of a gene.  
# in order to weight equally each gene
# blind = TRUE is used to have unbiased PCA, it is not looking at which group the sample is belonging to 


# count data trasformation with rlog
# for PCA and small samples is good to use rlog, because it is more robust than vst 
rld <- rlog(dds, blind=TRUE)

# directory where all the plots of step5 will be stored 
setwd(dir_step5)
pcaData <- plotPCA(rld, intgroup = c("Disease", "Genotype"), returnData=TRUE) # ntop = 500 as default 
percentVar <- round(100 * attr(pcaData, "percentVar"))
sample_names <- rownames(col_data)
ggplot(pcaData, aes(PC1, PC2, color=Disease, shape =Genotype )) +
  geom_point(size=3) +
  geom_text_repel(
    aes(label = sample_names), 
    max.overlaps = Inf,  
    show.legend = FALSE
  ) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_y_continuous(breaks = scales::breaks_width(20)) +
  coord_fixed()
ggsave("plotPCA.pdf")

# Count the knocked-out Ifnar and Ifngr genes to verify that the knockout was successful
# Ifnar1: ENSMUSG00000022967
# Ifnar2: ENSMUSG00000022971
plot_counts_single_gene(dds, "Ifnar1", "ENSMUSG00000022967", "KO_Ifnar1")
plot_counts_single_gene(dds, "Ifnar2", "ENSMUSG00000022971", "KO_Ifnar2")

# Ifngr1: ENSMUSG00000020009
# Ifngr2: ENSMUSG00000022965
plot_counts_single_gene(dds, "Ifngr1", "ENSMUSG00000020009", "KO_Ifngr1")
plot_counts_single_gene(dds, "Ifngr2", "ENSMUSG00000022965", "KO_Ifngr2")


#
# STEP 6: differential analysis 
#

#extract the results from dds with different contrast 
# double_knocked_out vs wild_type (with wild_type as reference) 
res_genotype <- results(dds, contrast = c("Genotype", "DKO", "WT")) 
# infected/case vs control (with control as reference) 
res_disease <- results(dds, contrast = c("Disease", "Case", "Control")) 
# interaction term 
res_interaction <- results(dds, name = "DiseaseCase.GenotypeDKO") 

# summary of the results 
summary(res_genotype)
summary(res_disease)
summary(res_interaction)

# --------------------------------------------------#
#                     VOLCANO PLOTS                 #
# --------------------------------------------------#

# directory specific for the LFC > 1 and padj < 0.05
dir_step6_LFC_1_padj_005 <- paste0(dir_step6, "/LFC_1_padj_0-05")
# directory specific for the LFC > 1 and padj < 0.01 
dir_step6_LFC_1_padj_001 <- paste0(dir_step6, "/LFC_1_padj_0-01")
# selecting the directory based on the adjusted p-value
if(padj_limit == 0.05){
  current_dir_step6 <- dir_step6_LFC_1_padj_005
} else if (padj_limit == 0.01) {
  current_dir_step6 <- dir_step6_LFC_1_padj_001
}


create_not_existing_dir(current_dir_step6)
setwd(current_dir_step6)

# ----------------------------------------------------------------------------------------------
# Volcano plot for double_knocked_out vs wild_type (with wild_type as reference)  
# ----------------------------------------------------------------------------------------------

# Add the column gene_id to results table  
res_genotype <- add_col_gene_id(res_genotype)
# Add a column express_level to specify if they are UP- or DOWN- regulated with a adjusted p-value = padj_limit
res_genotype <- add_col_express_level(res_genotype)


# --- both up-regulated and down-regulated 

# Add the column gene_label 
res_genotype_up_down_regulated <- add_col_gene_label_up_down_regulated(res_genotype)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_genotype_up_down_regulated <- res_genotype_up_down_regulated[which(res_genotype_up_down_regulated$gene_label != ""),]


volcano_ggplot(res_genotype_up_down_regulated, 
               'Volcano plots for DKO vs WT up/down-regulated', 
               'volcano_plot_genotype_up_down.pdf', 
               int_genes_genotype_up_down_regulated)

# --- only up-regulated 

# Add the column gene_label 
res_genotype_up_regulated <- add_col_gene_label_up_regulated(res_genotype)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_genotype_up_regulated <- res_genotype_up_regulated[which(res_genotype_up_regulated$gene_label != ""),]

volcano_ggplot(res_genotype_up_regulated, 
               'Volcano plots for DKO vs WT up-regulated', 
               'volcano_plot_genotype_up.pdf', 
               int_genes_genotype_up_regulated)

# --- only down-regulated 

# Add the column gene_label 
res_genotype_down_regulated <- add_col_gene_label_down_regulated(res_genotype)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_genotype_down_regulated <- res_genotype_down_regulated[which(res_genotype_down_regulated$gene_label != ""),]

volcano_ggplot(res_genotype_down_regulated, 
               'Volcano plots for DKO vs WT down-regulated', 
               'volcano_plot_genotype_down.pdf', 
               int_genes_genotype_down_regulated)



# ----------------------------------------------------------------------------------------------
# Volcano plot for infected/case vs control (with control as reference)   
# ----------------------------------------------------------------------------------------------

# Add the column gene_id to results table  
res_disease <- add_col_gene_id(res_disease)
# Add a column express_level to specify if they are UP- or DOWN- regulated with a adjusted p-value = padj_limit 
res_disease <- add_col_express_level(res_disease)


# --- both up-regulated and down-regulated 

# Add the column gene_label 
res_disease_up_down_regulated <- add_col_gene_label_up_down_regulated(res_disease)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_disease_up_down_regulated <- res_disease_up_down_regulated[which(res_disease_up_down_regulated$gene_label != ""),]


volcano_ggplot(res_disease_up_down_regulated, 
               'Volcano plots for Infected vs Control up/down-regulated', 
               'volcano_plot_disease_up_down.pdf', 
               int_genes_disease_up_down_regulated)

# --- only up-regulated 

# Add the column gene_label 
res_disease_up_regulated <- add_col_gene_label_up_regulated(res_disease)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_disease_up_regulated <- res_disease_up_regulated[which(res_disease_up_regulated$gene_label != ""),]

volcano_ggplot(res_disease_up_regulated, 
               'Volcano plots for Infected vs Control up-regulated', 
               'volcano_plot_disease_up.pdf', 
               int_genes_disease_up_regulated)

# --- only down-regulated 

# Add the column gene_label 
res_disease_down_regulated <- add_col_gene_label_down_regulated(res_disease)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_disease_down_regulated <- res_disease_down_regulated[which(res_disease_down_regulated$gene_label != ""),]

volcano_ggplot(res_disease_down_regulated, 
               'Volcano plots for Infected vs Control down-regulated', 
               'volcano_plot_disease_down.pdf', 
               int_genes_disease_down_regulated)

# ----------------------------------------------------------------------------------------------
# Volcano plot for interaction term 
# ----------------------------------------------------------------------------------------------

# Add the column gene_id to results table  
res_interaction <- add_col_gene_id(res_interaction)
# Add a column express_level to specify if they are UP- or DOWN- regulated with a adjusted p-value = padj_limit 
res_interaction <- add_col_express_level(res_interaction)

# --- both up-regulated and down-regulated 

# Add the column gene_label 
res_interaction_up_down_regulated <- add_col_gene_label_up_down_regulated(res_interaction)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_interaction_up_down_regulated <- res_interaction_up_down_regulated[which(res_interaction_up_down_regulated$gene_label != ""),]


volcano_ggplot(res_interaction_up_down_regulated, 
               'Volcano plots for interaction term up/down-regulated', 
               'volcano_plot_interaction_up_down.pdf', 
               int_genes_interaction_up_down_regulated)

# --- only up-regulated 

# Add the column gene_label 
res_interaction_up_regulated <- add_col_gene_label_up_regulated(res_interaction)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_interaction_up_regulated <- res_interaction_up_regulated[which(res_interaction_up_regulated$gene_label != ""),]

volcano_ggplot(res_interaction_up_regulated, 
               'Volcano plots for interaction term up-regulated', 
               'volcano_plot_interaction_up.pdf', 
               int_genes_interaction_up_regulated)

# --- only down-regulated 

# Add the column gene_label 
res_interaction_down_regulated <- add_col_gene_label_down_regulated(res_interaction)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_interaction_down_regulated <- res_interaction_down_regulated[which(res_interaction_down_regulated$gene_label != ""),]

volcano_ggplot(res_interaction_down_regulated, 
               'Volcano plots for interaction term down-regulated', 
               'volcano_plot_interaction_down.pdf', 
               int_genes_interaction_down_regulated)

#
# investigate the expression level of interesting genes 
# 

# plot counts of gene Gbp5 (ENSMUSG00000105504)
plot_counts_single_gene(dds, "Gbp5", "ENSMUSG00000105504", "Gbp5")


# plot counts of gene Gbp2 (ENSMUSG00000028270)
plot_counts_single_gene(dds, "Gbp2", "ENSMUSG00000028270", "Gbp2")

# plot counts of gene Tent5c (ENSMUSG00000044468)
plot_counts_single_gene(dds, "Tent5c", "ENSMUSG00000044468", "Tent5c")


#
#
# STEP 7: Overrepresentation analysis
#
#

# based on the LFC and padj's thresholds substitute the output_path in the enrichGO function 
dir_step7_LFC_1_padj_005 <- paste0(dir_step7, "/LFC_1_padj_0-05")
dir_step7_LFC_1_padj_001 <- paste0(dir_step7, "/LFC_1_padj_0-01")
# selecting the directory based on the adjusted p-value
if (padj_limit == 0.05){
  current_dir_step7 <- dir_step7_LFC_1_padj_005
} else if (padj_limit == 0.01){
  current_dir_step7 <- dir_step7_LFC_1_padj_001
}

create_not_existing_dir(current_dir_step7)

# -------------------------------------------------------------------------------
# identification of the biological function of the possibile interesting genes 
# -------------------------------------------------------------------------------

# --- in the contrast DKO vs WT 

# analysis upregulated and downregulated genes at the same time 
enrichGO_genotype_up_down_regulated <- enrichGO_plotting(int_genes_genotype_up_down_regulated, 
                                                  res_genotype_up_down_regulated, 
                                                  "DKO vs WT up and down regulated", 
                                                  current_dir_step7)
# analysis only upregulated genes
enrichGO_genotype_up_regulated <- enrichGO_plotting(int_genes_genotype_up_regulated, 
                                                  res_genotype_up_regulated, 
                                                  "DKO vs WT up regulated", 
                                                  current_dir_step7)
# analysis only downregulated genes 
enrichGO_genotype_down_regulated <- enrichGO_plotting(int_genes_genotype_down_regulated, 
                                                  res_genotype_down_regulated, 
                                                  "DKO vs WT down regulated", 
                                                  current_dir_step7)

# --- in the contrast Infected vs Control 

# analysis upregulated and downregulated genes at the same time 
enrichGO_disease_up_down_regulated <- enrichGO_plotting(int_genes_disease_up_down_regulated, 
                                                  res_disease_up_down_regulated, 
                                                  "Infected vs Control up and down regulated", 
                                                  current_dir_step7)
# analysis only upregulated genes
enrichGO_disease_up_regulated <- enrichGO_plotting(int_genes_disease_up_regulated, 
                                                  res_disease_up_regulated, 
                                                  "Infected vs Control up regulated", 
                                                  current_dir_step7)
# analysis only downregulated genes 
enrichGO_disease_down_regulated <- enrichGO_plotting(int_genes_disease_down_regulated, 
                                                  res_disease_down_regulated, 
                                                  "Infected vs Control down regulated", 
                                                  current_dir_step7)


# --- in the contrast interaction term 

# analysis upregulated and downregulated genes at the same time 
enrichGO_interaction_up_down_regulated <- enrichGO_plotting(int_genes_interaction_up_down_regulated, 
                                                  res_interaction_up_down_regulated, 
                                                  "Interaction term up and down regulated", 
                                                  current_dir_step7)
enrichGO_interaction_up_regulated <- enrichGO_plotting(int_genes_interaction_up_regulated, 
                                                  res_interaction_up_regulated, 
                                                  "Interaction term up regulated", 
                                                  current_dir_step7)
enrichGO_interaction_down_regulated <- enrichGO_plotting(int_genes_interaction_down_regulated, 
                                                  res_interaction_down_regulated, 
                                                  "Interaction term down regulated", 
                                                  current_dir_step7)

enrichGO_interaction_up_regulated@result[enrichGO_interaction_up_regulated@result$Description == "erythrocyte differentiation", ]
enrichGO_interaction_up_regulated@result[enrichGO_interaction_up_regulated@result$Description == "erythrocyte differentiation", ]



# TODO: arrivtata fino a qui 
# TODO: 21/12/25 - 23:32 



















##TODO !! THIS COMMENTS BELOW NEEDS TO GO IN THE FUTURE 

# last step of step 6 
# Based on the original publication, select 2-3 genes that are of particular interest and investigate their expression level. 
#You could use, for example, the normalised counts (see DESeq2::counts) 
#where the effect of between-sample differences in sequencing depth has been removed.
# counts function to see the normalised counts 
counts(dds)



# to convert the name of the gene you can use the package that is in step 7 
# there are different genes that can have the same symble, so we should work with them only in the end 
# we have to translate them only when we want to visualise them 

# package org.Mm.eg.

##TODO !! THIS COMMENTS NEEDS TO GO IN THE FUTURE [END]












