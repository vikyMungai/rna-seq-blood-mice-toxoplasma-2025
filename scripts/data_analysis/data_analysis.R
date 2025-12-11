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


#
# STEP 5: Exploratory data analysis
#

# Directory where the file "formatted_feature_count.txt" is 
setwd("/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/feature_counts")

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
col_data <- read.csv("DESeq_colData_1_condition.csv", 
                     header = TRUE, 
                     sep = "\t")

# order the rows by sample name 
col_data <- col_data[order(col_data$Sample), ]
# select the column "Sample" as rownames 
col_data <- data.frame(col_data, row.names = 1)

# check if the columns name of the sample in the counts_data matches 
# to the row names in colData
all(colnames(counts_data) %in% rownames(col_data))

# check if the columns name of the sample in the counts_data are in the same order 
# of the row names in colData
all(colnames(counts_data) == rownames(col_data))

# converting column Group of col_data into factors before running DESeqDataSetFromMatrix
col_data$Group <- factor(col_data$Group)

# create a DESeqDataSet object with single factor 
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_data,
                              design = ~ Group)

# by setting the factor level we explicit say which level we want as reference
# otherwise it will be chosen alphabetically 
# for the PCA plot is not important, but it is for the Differential Analysis 
# in case the contrast is not specified
dds$Group <- relevel(dds$Group, ref = "Blood_WT_Control")

# run DESeq
dds <- DESeq(dds)

# it is important to not have a correlation between the mean expression and the variance of a gene.  
# in order to weight equally each gene
# blind = TRUE is used to have unbiased PCA, it is not looking at which group the sample is belonging to 

# count data trasformation with rlog
# for PCA and small samples is good to use rlog, because it is more robust than vst 
rld <- rlog(dds, blind=TRUE)

# plotting PCA 
pcaData <- plotPCA(rld, intgroup = "Group", returnData=TRUE) # ntop = 500 as default 
percentVar <- round(100 * attr(pcaData, "percentVar"))
sample_names <- rownames(col_data)
ggplot(pcaData, aes(PC1, PC2, color=Group)) +
  geom_text(
    mapping = aes(label = sample_names), 
    # we are interested in outliers so it is not a problem if we cannot see all the samples names
    check_overlap = TRUE,  
    position = position_nudge(x = NULL, y = 3.5)
  ) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  xlim(-70, 68) + 
  scale_y_continuous(breaks = scales::breaks_width(20)) +
  coord_fixed()


#
# STEP 6: differential analysis 
#

# ----
# extract the results from dds with different contrast 
# ----

# DKO infected vs WT infected (with WT infected as reference) 
res_DKO_Case_vs_WT_Case <- results(dds, contrast=c("Group", "Blood_DKO_Case", "Blood_WT_Case"))
# WT infected vs WT not infected (with WT not infected as reference) 
res_WT_Case_vs_WT_Control <- results(dds, contrast=c("Group","Blood_WT_Case", "Blood_WT_Control"))
# DKO infected vs DKO not infected  (with WT infected as reference) 
res_DKO_Case_vs_DKO_Control <- results(dds, contrast=c("Group", "Blood_DKO_Case", "Blood_DKO_Control"))
# DKO not infected vs WT not infected (with WT not infected as reference) 
res_DKO_Control_vs_WT_Control <- results(dds, contrast=c("Group","Blood_DKO_Control", "Blood_WT_Control"))

# summary of the results 
summary(res_DKO_Case_vs_WT_Case)
summary(res_WT_Case_vs_WT_Control)
summary(res_DKO_Case_vs_DKO_Control)
summary(res_DKO_Control_vs_WT_Control)

# ----
# Volcano plots 
#-----

add_col_gene_name <- function(results_DESeq){
  #   Add the column 'gene_name' with all the gene names 
  #
  #   Parameter: 
  #     results_DESeq: the results of DESeq 
  #
  #   Returns: 
  #     the results_DESeq with the new column added 
  results_DESeq$gene_name <- rownames(results_DESeq)
  return (
    results_DESeq
  )
}

add_col_express_level<- function(results_DESeq, LFC_thresholds = c(-1,1), p_value_threshold = 0.05){
  #   Add the column 'express_level' with all the gene names 
  #
  #   Parameter: 
  #     results_DESeq (DESeqResults): the results of DESeq 
  #     LFC_thresholds (numeric): a numeric vector with the higher and lower thresholds 
  #                               default value: c(-1, 1) 
  #     p_value_threshold (numeric): a numeric vector with the threshold for the p-value
  #                               default value: 0.05 (5%)
  #
  #   Returns: 
  #     the results_DESeq with the new column added 
  
  results_DESeq$express_level <- "Not significant"
  # if log2Foldchange > 1 and pvalue < 0.05, set as "Upregulated" 
  results_DESeq$express_level[results_DESeq$log2FoldChange > LFC_thresholds[2] & results_DESeq$pvalue < p_value_threshold] <- "Upregulated"
  # if log2Foldchange < -1 and pvalue < 0.05, set as "Downregulated" 
  results_DESeq$express_level[results_DESeq$log2FoldChange < LFC_thresholds[1] & results_DESeq$pvalue < p_value_threshold] <- "Downregulated"
  # Set the labels of the possible interesting genes (based on their p-value
  # they also need to be a biological significant, so Downregulated or Upregulated)
  results_DESeq$gene_label <- ifelse(results_DESeq$gene_name %in% head(results_DESeq[order(results_DESeq$padj), "gene_name"], 30 ) & (results_DESeq$express_level == "Downregulated" | results_DESeq$express_level == "Upregulated" ), results_DESeq$gene_name, "")
    return (
    results_DESeq
  )
}

add_col_express_level <- function(results_DESeq, LFC_thresholds = c(-1,1), p_value_threshold = 0.05){
  #   Add the column 'express_level' with all the gene names 
  #
  #   Parameter: 
  #     results_DESeq (DESeqResults): the results of DESeq 
  #     LFC_thresholds (numeric): a numeric vector with the higher and lower thresholds 
  #                               default value: c(-1, 1) 
  #     p_value_threshold (numeric): a numeric vector with the threshold for the p-value
  #                               default value: 0.05 (5%)
  #
  #   Returns: 
  #     the results_DESeq with the new column added 
  
  results_DESeq$express_level <- "Not significant"
  # if log2Foldchange > 1 and pvalue < 0.05, set as "Upregulated" 
  results_DESeq$express_level[results_DESeq$log2FoldChange > LFC_thresholds[2] & results_DESeq$pvalue < p_value_threshold] <- "Upregulated"
  # if log2Foldchange < -1 and pvalue < 0.05, set as "Downregulated" 
  results_DESeq$express_level[results_DESeq$log2FoldChange < LFC_thresholds[1] & results_DESeq$pvalue < p_value_threshold] <- "Downregulated"
  
  return (
    results_DESeq
  )
}

#TODO: remove limit from p-value (put a more strict value for the p-value )
add_col_gene_label_up_down_regulated <- function(results_DESeq, LFC_thresholds = c(-1,1), p_value_threshold = 0.05){
  #   Add the column 'gene_label'. Only the interesting genes will have a label. 
  #
  #   Parameter: 
  #     results_DESeq (DESeqResults): the results of DESeq 
  #     LFC_thresholds (numeric): a numeric vector with the higher and lower thresholds 
  #                               default value: c(-1, 1) 
  #     p_value_threshold (numeric): a numeric vector with the threshold for the p-value
  #                               default value: 0.05 (5%)
  #
  #   Returns: 
  #     the results_DESeq with the new column added 
  
  # Set the labels of the possible interesting genes (based on their p-value
  # they also need to be a biological significant, so Downregulated or Upregulated)
  results_DESeq$gene_label <- ifelse(results_DESeq$gene_name %in% results_DESeq[ which(results_DESeq$padj), "gene_name"] & (results_DESeq$express_level == "Downregulated" | results_DESeq$express_level == "Upregulated" ), results_DESeq$gene_name, "")
  return (
    results_DESeq
  )
}

add_col_gene_label_up_regulated <- function(results_DESeq, LFC_threshold = 1, p_value_threshold = 0.05){
  #   Add the column 'gene_label'. Only the interesting up-regulated genes will have a label. 
  #
  #   Parameter: 
  #     results_DESeq (DESeqResults): the results of DESeq 
  #     LFC_thresholds (numeric): the threshold for the FLC  
  #                               default value:  1 
  #     p_value_threshold (numeric): a numeric vector with the threshold for the p-value
  #                               default value: 0.05 (5%)
  #
  #   Returns: 
  #     the results_DESeq with the new column added 
  
  # Set the labels of the possible interesting genes (based on their p-value
  # they also need to be a biological significant, so Down-regulated or Up-regulated)
  results_DESeq$gene_label <- ifelse(results_DESeq$gene_name %in% head(results_DESeq[order(results_DESeq$padj), "gene_name"], 30 ) & (results_DESeq$express_level == "Upregulated" ), results_DESeq$gene_name, "")
  return (
    results_DESeq
  )
}

add_col_gene_label_up_down_regulated <- function(results_DESeq, LFC_threshold = -1, p_value_threshold = 0.05){
  #   Add the column 'gene_label'. Only the interesting down-regulated genes will have a label. 
  #
  #   Parameter: 
  #     results_DESeq (DESeqResults): the results of DESeq 
  #     LFC_thresholds (numeric): the threshold for the FLC  
  #                               default value:  -1 
  #     p_value_threshold (numeric): a numeric vector with the threshold for the p-value
  #                               default value: 0.05 (5%)
  #
  #   Returns: 
  #     the results_DESeq with the new column added 
  
  # Set the labels of the possible interesting genes (based on their p-value
  # they also need to be a biological significant, so Downregulated or Upregulated)
  results_DESeq$gene_label <- ifelse(results_DESeq$gene_name %in% head(results_DESeq[order(results_DESeq$padj), "gene_name"], 30 ) & (results_DESeq$express_level == "Downregulated"), results_DESeq$gene_name, "")
  return (
    results_DESeq
  )
}



volcano_ggplot <- function(results_DESeq, plot_title, log2_fold_change = 1, p_value = 0.05){
    #   Print the volcano ggplot of the results_DESeq with specific settings
    #
    #   Parameter: 
    #     results_DESeq: the results of DESeq 
    #     plot_title: the title of the plot 
    #
    #   Returns: 
    #     None - show the plot in the 'Plots' window
  
    ggplot(data = results_DESeq, aes(x = log2FoldChange, y = -log10(pvalue), col = express_level)) + 
      geom_vline(xintercept = c(-log2_fold_change, log2_fold_change), col = "gray", linetype = 'dashed') + 
      geom_hline(yintercept = -log10(p_value), col = "gray", linetype = 'dashed') + 
      geom_point(size = 1) +
      scale_color_manual(values = c("blue", "grey", "plum"), 
                         labels = c("Downregulated", "Not significant", "Upregulated")) + 
      coord_cartesian(ylim = c(0, 300), xlim = c(-20, 15)) + 
      labs(color = 'Gene express level',
           x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
      ggtitle(plot_title) +
      # max overlaps was chosen Inf to not lose labels 
      geom_text_repel(aes(label = gene_label), max.overlaps = Inf, show.legend = FALSE)
  

}
   


# --- Volcano plot for the DKO infected vs WT infected

# Add the column gene_name to res_DKO_Case_vs_WT_Case 
res_DKO_Case_vs_WT_Case <- add_col_gene_name(res_DKO_Case_vs_WT_Case)

# Add a column express_level to specify if they are UP- or DOWN- regulated with a significant p-value 
res_DKO_Case_vs_WT_Case <- add_col_express_level(res_DKO_Case_vs_WT_Case)

# Add the column gene_label for both up-regulated and down-regulated 
res_DKO_Case_vs_WT_Case <- add_col_gene_label_up_down_regulated(res_DKO_Case_vs_WT_Case, FLC_thresholds, p_value_threshold)

volcano_ggplot(res_DKO_Case_vs_WT_Case, 'Volcano plots for DKO not infected vs WT not infected')

# the ENSEMBL gene IDs that were filtered with the LFC = 1 and p-value = 0.05 
int_genes_DKO_Case_vs_WT_Case <- res_DKO_Case_vs_WT_Case[which(res_DKO_Case_vs_WT_Case$gene_label != ""),]


# --- Volcano plot for the WT infected vs WT not infected

FLC_thresholds <- c(-1, 1)
p_value_threshold <- 0.05


# Add the column gene_name to res_WT_Case_vs_WT_Control
res_WT_Case_vs_WT_Control <- add_col_gene_name(res_WT_Case_vs_WT_Control)

# Add a column express_level to specify if they are UP- or DOWN- regulated with a significant p-value 
res_WT_Case_vs_WT_Control <- add_col_express_level(res_WT_Case_vs_WT_Control, FLC_thresholds, p_value_threshold)
# Add the column gene_label for both up-regulated and down-regulated 
res_WT_Case_vs_WT_Control <- add_col_gene_label_up_down_regulated(res_WT_Case_vs_WT_Control, FLC_thresholds, p_value_threshold)

volcano_ggplot(res_WT_Case_vs_WT_Control, 'Volcano plots for WT infected vs WT not infected')

# the ENSEMBL gene IDs that were filtered with the LFC = 1 and p-value = 0.05 
int_genes_WT_Case_vs_WT_Control <- res_WT_Case_vs_WT_Control[which(res_WT_Case_vs_WT_Control$gene_label != ""),]

# --- Volcano plot for the DKO infected vs DKO not infected

# Add the column gene_name to res_DKO_Case_vs_DKO_Control
res_DKO_Case_vs_DKO_Control <- add_col_gene_name(res_DKO_Case_vs_DKO_Control)

# Add a column express_level to specify if they are UP- or DOWN- regulated with a significant p-value 
res_DKO_Case_vs_DKO_Control <- add_col_express_level(res_DKO_Case_vs_DKO_Control, FLC_thresholds, p_value_threshold)
# Add the column gene_label for both up-regulated and down-regulated 
res_DKO_Case_vs_DKO_Control <- add_col_gene_label_up_down_regulated(res_DKO_Case_vs_DKO_Control, FLC_thresholds, p_value_threshold)

volcano_ggplot(res_DKO_Case_vs_DKO_Control, 'Volcano plots for DKO infected vs DKO not infected')

# the ENSEMBL gene IDs that were filtered with the LFC = 1 and p-value = 0.05 
int_genes_DKO_Case_vs_DKO_Control <- res_DKO_Case_vs_DKO_Control[which(res_DKO_Case_vs_DKO_Control$gene_label != ""),]

# --- Volcano plot for the DKO not infected vs WT not infected

# Add the column gene_name to res_DKO_Control_vs_WT_Control
res_DKO_Control_vs_WT_Control <- add_col_gene_name(res_DKO_Control_vs_WT_Control)

# Add a column express_level to specify if they are UP- or DOWN- regulated with a significant p-value 
res_DKO_Control_vs_WT_Control <- add_col_express_level(res_DKO_Control_vs_WT_Control, FLC_thresholds, p_value_threshold)
# Add the column # Add the column gene_label for both up-regulated and down-regulated 
res_DKO_Case_vs_DKO_Control <- add_col_gene_label_up_down_regulated(res_DKO_Case_vs_DKO_Control, FLC_thresholds, p_value_threshold)

volcano_ggplot(res_DKO_Control_vs_WT_Control, 'Volcano plots for DKO not infected vs WT not infected')

# the ENSEMBL gene IDs that were filtered with the LFC = 1 and p-value = 0.05 
int_genes_DKO_Control_vs_WT_Control <- res_DKO_Control_vs_WT_Control[which(res_DKO_Control_vs_WT_Control$gene_label != ""),]


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
  
  enrichGO_results <- clusterProfiler::enrichGO(
    gene = interesting_genes$gene_name,
    OrgDb = "org.Mm.eg.db",
    keyType = "ENSEMBL",
    ont = "BP", # put all of them (ALL or a vector )
    universe = dds_results$gene_name, #filter out the genes that don't have a adjust p-value before passing them to the universe 
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
  
  if(nrow(enrichGO_results) == 0) {
    # if there is not enriched gene than we don't plot the results for the heatplot 
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

output_path <- "/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/overrepresentation_analysis"

# identification of the biological function of the possibile interesting genes 

# --- in the contrast DKO infected vs WT infected
enrichGO_DKO_Case_vs_WT_Case <- enrichGO_plotting(int_genes_DKO_Case_vs_WT_Case, 
                                                  res_DKO_Case_vs_WT_Case, 
                                                  "DKO infected vs WT infected", 
                                                  output_path)

# --- in the contrast WT infected vs WT not infected
enrichGO_WT_Case_vs_WT_Control <- enrichGO_plotting(int_genes_WT_Case_vs_WT_Control, 
                                                  res_WT_Case_vs_WT_Control, 
                                                  "WT infected vs WT not infected", 
                                                  output_path)

# --- in the contrast DKO infected vs DKO not infected
enrichGO_DKO_Case_vs_DKO_Control <- enrichGO_plotting(int_genes_DKO_Case_vs_DKO_Control, 
                                                      res_DKO_Case_vs_DKO_Control, 
                                                    "DKO infected vs DKO not infected", 
                                                    output_path)

# --- in the contrast DKO not infected vs WT not infected
enrichGO_DKO_Control_vs_WT_Control <- enrichGO_plotting(int_genes_DKO_Control_vs_WT_Control, 
                                                        res_DKO_Control_vs_WT_Control, 
                                                        "DKO not infected vs WT not infected", 
                                                        output_path)








