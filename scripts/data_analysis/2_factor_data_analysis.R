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

# ----------------------------------------------------------------------------------
# FUNCTIONS 
# ----------------------------------------------------------------------------------

create_not_existing_dir <- function(dir_name){
    #   Create the directory only if it does not exist
    #
    #   Parameter: 
    #     dir_name: absolute path of the directory 
    #
    #   Returns: 
    #     None 
  if (!dir.exists(dir_name)) {
    # create also the parent's folders if necessary 
    dir.create(dir_name, recursive = TRUE)
    warning(paste0("The directory ", dir_name, " was created."))
  }
}

### input files 
dir_input <- "/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/feature_counts"
### directories to store plots
dir_step5 <- "/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/2_factors/5_exploratory_data_analysis"
dir_step6 <- "/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/2_factors/6_differential_analysis"
dir_step7 <- "/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/2_factors/7_overrepresentation_analysis"

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

# LFC and padj thresholds 
LFC_limit = 1 
padj_limit = 0.01

# directory specific for the LFC > 1 and padj < 0.05
dir_step6_LFC_1_padj_005 <- paste0(dir_step6, "/LFC_1_padj_0-05")
# directory specific for the LFC > 1 and padj < 0.01 
dir_step6_LFC_1_padj_001 <- paste0(dir_step6, "/LFC_1_padj_0-01")
# selecting the directory based on the adjusted p-value
if(padj_limit == 0.05){
  current_dir_step6 <- dir_step6_LFC_1_padj_005
} else {
  current_dir_step6 <- dir_step6_LFC_1_padj_001
}

setwd(current_dir_step6)

create_not_existing_dir(current_dir_step6)

# ----------------------------------------------------------------------------------------------
# functions used in step 6
# ----------------------------------------------------------------------------------------------

add_col_gene_id <- function(results_DESeq){
  #   Add the column 'gene_id' with all the gene names 
  #
  #   Parameter: 
  #     results_DESeq: the results of DESeq 
  #
  #   Returns: 
  #     the results_DESeq with the new column added 
  results_DESeq$gene_id <- rownames(results_DESeq)
  return (
    results_DESeq
  )
}


add_col_express_level <- function(results_DESeq, LFC_threshold = LFC_limit, padj_threshold = padj_limit){
  #   Add the column 'express_level' to specify the expression gene level of the 
  #
  #   Parameter: 
  #     results_DESeq (DESeqResults): the results of DESeq 
  #     LFC_threshold (numeric): threshold for LogFoldChange
  #     padj_threshold (numeric): the threshold for the adjusted p-value
  #
  #
  #   Returns: 
  #     the results_DESeq with the new column added 
  
  results_DESeq$express_level <- "Not significant"
  # if log2Foldchange > LFC_threshold and padj < padj_threshold, set as "Upregulated" 
  results_DESeq$express_level[results_DESeq$log2FoldChange > LFC_threshold & results_DESeq$padj < padj_threshold] <- "Upregulated"
  # if log2Foldchange < -LFC_threshold and padj < padj_threshold, set as "Downregulated" 
  results_DESeq$express_level[results_DESeq$log2FoldChange < -LFC_threshold & results_DESeq$padj < padj_threshold] <- "Downregulated"
  
  return (
    results_DESeq
  )
}

add_col_gene_label_up_down_regulated <- function(results_DESeq, padj_threshold = padj_limit){
  #   Add the column 'gene_label', the ENSEMBL gene ID if it is an interesting gene (Upregulated or downregulated), otherwise empty string 
  #
  #   Parameter: 
  #     results_DESeq (DESeqResults): the results of DESeq 
  #     padj_threshold (numeric): the threshold for the adjusted p-value
  #
  #
  #   Returns: 
  #     the results_DESeq with the new column added 
  
  # Set the labels of the possible interesting genes (based on their p-value
  # they also need to be a biological significant, so Down-regulated or Up-regulated)
  results_DESeq$gene_label <- ifelse( results_DESeq$padj < padj_threshold & (results_DESeq$express_level == "Downregulated" | results_DESeq$express_level == "Upregulated" ), results_DESeq$gene_id, "")
  return (
    results_DESeq
  )
}

add_col_gene_label_up_regulated <- function(results_DESeq, padj_threshold = padj_limit){
  #   Add the column 'gene_label', the ENSEMBL gene ID if it is an interesting gene (Up-regulated only), otherwise empty string 
  #
  #   Parameter: 
  #     results_DESeq (DESeqResults): the results of DESeq 
  #     padj_threshold (numeric): the threshold for the adjusted p-value
  #
  #   Returns: 
  #     the results_DESeq with the new column added 
  
  # Set the labels of the possible interesting up-regulated genes (based on their p-value)
  results_DESeq$gene_label <- ifelse( results_DESeq$padj < padj_threshold & (results_DESeq$express_level == "Upregulated" ), results_DESeq$gene_id, "")
  return (
    results_DESeq
  )
}

add_col_gene_label_down_regulated <- function(results_DESeq, padj_threshold = padj_limit){
  #   Add the column 'gene_label', the ENSEMBL gene ID if it is an interesting gene (Down-regulated only), otherwise empty string 
  #
  #   Parameter: 
  #     results_DESeq (DESeqResults): the results of DESeq 
  #     padj_threshold (numeric): a numeric vector with the threshold for the p-value
  #
  #   Returns: 
  #     the results_DESeq with the new column added 
  
  # Set the labels of the possible interesting down-regulated genes (based on their p-value)
  results_DESeq$gene_label <- ifelse(results_DESeq$padj < padj_threshold & (results_DESeq$express_level == "Downregulated"), results_DESeq$gene_id, "")
  return (
    results_DESeq
  )
}


volcano_ggplot <- function(results_DESeq, plot_title, output_file, interesting_genes = c(), log2_fold_change = LFC_limit, padj_value = padj_limit){
  #   Print the volcano ggplot of the results_DESeq with specific settings
  #
  #   Parameter: 
  #     results_DESeq: the results of DESeq 
  #     plot_title: the title of the plot 
  #
  #   Returns: 
  #     character - vector with the list of the interesting genes
  
  
  # TODO: mettere a posto come visualizzare i gene_symbols e non ENSEMBL gene IDs 
  
  # only plotting a subset of the interesting genes 
  top_genes <- head(interesting_genes[order(interesting_genes$padj), "gene_id"], 20)
  # gene_symbols column rappresents the filtering for the labels in the plot 
  results_DESeq$gene_symbols <- ifelse( results_DESeq$gene_id %in% top_genes,
                                        results_DESeq$gene_id, "")
  # TODO: for DB 
  print("for DB")
  print( results_DESeq$gene_symbols[which(results_DESeq$gene_symbols != "")])
  
  # obtain the gene symbols only for the genes I want to plot 
  symbols <- mapIds(org.Mm.eg.db, keys =  top_genes, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  # replace NA with "" for the ggplot 
  symbols[is.na(symbols)] <- ""
  
  all_symbols <- ifelse(results_DESeq$gene_symbols != "", 
                        symbols[results_DESeq$gene_symbols], "")
  
  ggplot(data = results_DESeq, aes(x = log2FoldChange, y = -log10(pvalue), col = express_level)) + 
    geom_vline(xintercept = c(-log2_fold_change, log2_fold_change), col = "gray", linetype = 'dashed') + 
    geom_hline(yintercept = padj_value, col = "gray", linetype = 'dashed') + 
    geom_point(size = 1) +
    scale_color_manual(values = c("blue", "grey", "plum"), 
                       labels = c("Downregulated", "Not significant", "Upregulated")) + 
    coord_cartesian(ylim = c(0, 300), xlim = c(-20, 15)) + 
    labs(color = 'Gene express level', x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
    ggtitle(plot_title) +
    # max overlaps was chosen Inf to not lose labels 
    geom_text_repel(aes(label = all_symbols), size = 2, max.overlaps = Inf, show.legend = FALSE)
  
  ggsave(output_file)
  
  return(
    top_genes
  )
  
}

# summary of the results 
# res_genotype
# res_disease
# res_interaction


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

# investigate the expression level of interesting genes 
# TODO: you have to use the function counts

#
#
# STEP 7 
#
#


enrichGO_plotting <- function(interesting_genes, dds_results, contrast, output_path){
  #   execute the overrepresentation analysis and show the results with different 
  #   plots
  #
  #   Parameter: 
  #     interesting_genes (DESeqResults): selected interest genes in the contrast
  #     dds_results (DESeqResults): all the genes identified in the contrast 
  #     contrast (str): the contrast that will be plotted in the title 
  #     output_path (str): absolute path to save the plots
  #
  #   Returns: 
  #     enrichResult - enrichGO_results 
  #     all the plots saved in one pdf the output_path 
  
  
  
  #filter out the genes that don't have an adjust p-value 
  dds_results_filtered <- dds_results[!is.na(dds_results$padj),]
  
  print("for DB:")
  print(paste("Processing contrast:", contrast))
  print(paste("Interesting genes:", length(interesting_genes$gene_id)))
  print(paste("Universe size:", length(dds_results_filtered$gene_id)))
  
  enrichGO_results <- clusterProfiler::enrichGO(
    gene = interesting_genes$gene_id,
    OrgDb = "org.Mm.eg.db",
    keyType = "ENSEMBL",
    ont = "ALL", 
    universe = dds_results_filtered$gene_id, 
    maxGSSize = 500 # this is the default value 
  )
  
  print("DB: everything is fine so far")
  print(class(enrichGO_results))
  
  file_name <- gsub(" ", "_", contrast)
  pdf(paste0(output_path, "/", file_name, ".pdf"), width = 20, height = 10 )
  
  # plotting with bar plots the enrichment result 
  print(
    barplot(enrichGO_results, 
            showCategory=10, 
            title = paste0("Biological function for interesting genes in contrast ", contrast))
  )
  
  print("DB: printed barplot")
  
  #plotting with dot plot 
  print(
    enrichplot::dotplot(enrichGO_results, 
                        showCategory=10) +  
      ggtitle(paste0("Dotplot for contrast ", contrast))
  )
  
  print("DB: printed dotplot")
  
  # if there is not enriched gene than we don't plot the results for the heatplot 
  if(nrow(enrichGO_results) != 0) {
    # extract the foldChange of the interesting genes 
    foldChange<- interesting_genes$log2FoldChange
    names(foldChange) <- rownames(interesting_genes)
    # heatplot 
    print(
      enrichplot::heatplot(enrichGO_results, 
                           showCategory=5, 
                           foldChange = foldChange ) + 
        ggtitle(paste0("Heatplot for contrast ", contrast))
    )
    print("DB: printed heatplot")
  }
  
  
  dev.off()
  
  return (
    enrichGO_results
  )
}

# based on the LFC and padj's thresholds substitute the output_path in the enrichGO function 
dir_step7_LFC_1_padj_005 <- paste0(dir_step7, "/LFC_1_padj_0-05")
dir_step7_LFC_1_padj_001 <- paste0(dir_step7, "/LFC_1_padj_0-01")
# selecting the directory based on the adjusted p-value
if(padj_limit == 0.05){
  current_dir_step7 <- dir_step7_LFC_1_padj_005
} else {
  current_dir_step7 <- dir_step7_LFC_1_padj_001
}

create_not_existing_dir(current_dir_step7)

# -------------------------------------------------------------------------------
# identification of the biological function of the possibile interesting genes 
# -------------------------------------------------------------------------------

# --- in the contrast DKO vs WT 

# analysis upregulated and downregulated genes at the same time 
enrichGO_genotype <- enrichGO_plotting(int_genes_genotype_up_down_regulated, 
                                                  res_genotype_up_down_regulated, 
                                                  "DKO vs WT up and down regulated", 
                                                  current_dir_step7)
# analysis only upregulated genes
enrichGO_genotype <- enrichGO_plotting(int_genes_genotype_up_regulated, 
                                                  res_genotype_up_regulated, 
                                                  "DKO vs WT up regulated", 
                                                  current_dir_step7)
# analysis only downregulated genes 
enrichGO_genotype <- enrichGO_plotting(int_genes_genotype_down_regulated, 
                                                  res_genotype_down_regulated, 
                                                  "DKO vs WT down regulated", 
                                                  current_dir_step7)

# --- in the contrast Infected vs Control 

# analysis upregulated and downregulated genes at the same time 
enrichGO_disease <- enrichGO_plotting(int_genes_disease_up_down_regulated, 
                                                  res_disease_up_down_regulated, 
                                                  "Infected vs Control up and down regulated", 
                                                  current_dir_step7)
# analysis only upregulated genes
enrichGO_disease <- enrichGO_plotting(int_genes_disease_up_regulated, 
                                                  res_disease_up_regulated, 
                                                  "Infected vs Control up regulated", 
                                                  current_dir_step7)
# analysis only downregulated genes 
enrichGO_disease <- enrichGO_plotting(int_genes_disease_down_regulated, 
                                                  res_disease_down_regulated, 
                                                  "Infected vs Control down regulated", 
                                                  current_dir_step7)


# --- in the contrast interaction term 

# analysis upregulated and downregulated genes at the same time 
enrichGO_interaction <- enrichGO_plotting(int_genes_interaction_up_down_regulated, 
                                                  res_interaction_up_down_regulated, 
                                                  "Interaction term up and down regulated", 
                                                  current_dir_step7)
enrichGO_interaction <- enrichGO_plotting(int_genes_interaction_up_regulated, 
                                                  res_interaction_up_regulated, 
                                                  "Interaction term up regulated", 
                                                  current_dir_step7)
enrichGO_interaction <- enrichGO_plotting(int_genes_interaction_down_regulated, 
                                                  res_interaction_down_regulated, 
                                                  "Interaction term down regulated", 
                                                  current_dir_step7)

# TODO: ho rifatto tutti gli step che avevo fatto anche nel 1 factor contrast, ora devo provare a runnarlo da capo 
# TODO: poi devo aggiungere i pezzi che mi sono segnata 
# TODO: 21/12/25 - 23:32 



















########################################################### OLD VERSION 20/12/25 - 10:58  
# Add the column gene_name to res_genotype 
res_genotype$gene_name <- rownames(res_genotype)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated 
# (log2fc respectively positive or negative)
res_genotype$express_level <- "Not significant"
# if log2Foldchange > 1 and pvalue < 0.05, set as "UP"  
res_genotype$express_level[res_genotype$log2FoldChange > LFC_thresholds[1] & res_genotype$pvalue < p_value_threshold] <- "Upregulated"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN" ( put 1 so they are twice downregulated and significance statistic 0.05)
res_genotype$express_level[res$log2FoldChange < LFC_thresholds[0] & res$pvalue < p_value_threshold] <- "Downregulated"
# Set the labels of the possible interesting genes 
res_genotype$gene_label <- ifelse(res_genotype$gene_name %in% head(res_genotype[order(res_genotype$padj), "gene_name"], 3 ), res_genotype$gene_name, "")


ggplot(data = res_genotype, aes(x = log2FoldChange, y = -log10(pvalue), col = express_level)) + 
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "grey", "plum"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) + 
  coord_cartesian(ylim = c(0, 150), xlim = c(-20, 15)) + 
  labs(color = 'Gene express level',
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  ggtitle('Volcano plots for Genotype DKO vs WT') + 
  geom_text_repel(aes(label = gene_label) , max.overlaps = Inf)

# --- Volcano plot for the disease comparison

# Add the column gene_name to res_disease 
res_disease$gene_name <- rownames(res_disease)
# Add a column to the data frame to specify if they are UP- or DOWN- regulated with a significant p-value 
# (log2fc respectively positive or negative)
res_disease$express_level <- "Not significant"
# if log2Foldchange > 1 and pvalue < 0.05, set as "Upregulated" 
res_disease$express_level[res_disease$log2FoldChange > LFC_thresholds[1] & res_disease$pvalue < p_value_threshold] <- "Upregulated"
# if log2Foldchange < -1 and pvalue < 0.05, set as "Downregulated" ( put 1 so they are twice downregulated and significance statistic 0.05)
res_disease$express_level[res_disease$log2FoldChange < LFC_thresholds[0] & res_disease$pvalue < p_value_threshold] <- "Downregulated"
# Set the labels of the possible interesting genes (based on their p-value)
res_disease$gene_label <- ifelse(res_disease$gene_name %in% head(res_disease[order(res_disease$padj), "gene_name"], 10 ), res_disease$gene_name, "")

 
ggplot(data = res_disease, aes(x = log2FoldChange, y = -log10(pvalue), col = express_level)) + 
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') + 
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "grey", "plum"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) + 
  coord_cartesian(ylim = c(0, 150), xlim = c(-20, 15)) + 
  labs(color = 'Gene express level',
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  ggtitle('Volcano plots for Case vs Control') + 
  geom_text_repel(aes(label = gene_label) , max.overlaps = Inf)


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












