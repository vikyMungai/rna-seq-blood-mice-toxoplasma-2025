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

### directories to store plots
dir_step5 <- "/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/5_exploratory_data_analysis"
dir_step6 <- "/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/6_differential_analysis"
dir_step7 <- "/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/7_overrepresentation_analysis"

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

# directory where all the plots of step5 will be stored 
setwd(dir_step5)

# plotting PCA 
pcaData <- plotPCA(rld, intgroup = "Group", returnData=TRUE) # ntop = 500 as default 
percentVar <- round(100 * attr(pcaData, "percentVar"))
sample_names <- rownames(col_data)
ggplot(pcaData, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  geom_text_repel(
    aes(label = sample_names), 
    max.overlaps = Inf,  
    show.legend = FALSE
  ) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  xlim(-70, 68) + 
  scale_y_continuous(breaks = scales::breaks_width(20)) +
  coord_fixed()
ggsave("plotPCA.pdf")

# Count the knocked-out Ifnar and Ifngr genes to verify that the knockout was successful
# Ifnar1: ENSMUSG00000022967
# Ifnar2: ENSMUSG00000022971
# dds_gene_symbols <- 
d <- plotCounts(dds, gene="ENSMUSG00000022967", intgroup="Group", 
                returnData=TRUE)
ggplot(d, aes(x=Group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks = scales::log_breaks(n = 6, base = 10)) + 
  xlab("Genotype and Condition") + 
  ylab("Normalised counts") + 
  ggtitle("Ifnar1 - ENSMUSG00000022967 counts") + 
  theme(plot.title = element_text(margin = margin(t = 10, b = 10)))
ggsave("Normalised_counts_KO_Ifnar1.pdf")
d <- plotCounts(dds, gene="ENSMUSG00000022971", intgroup="Group", 
                returnData=TRUE)
ggplot(d, aes(x=Group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks = scales::log_breaks(n = 6, base = 10)) + 
  xlab("Genotype and Condition") + 
  ylab("Normalised counts") + 
  ggtitle("Ifnar2 - ENSMUSG00000022971 counts") + 
  theme(plot.title = element_text(margin = margin(t = 10, b = 10)))
ggsave("Normalised_counts_KO_Ifnar2.pdf")

# Ifngr1: ENSMUSG00000020009
# Ifngr2: ENSMUSG00000022965
# dds_gene_symbols <- 
d <- plotCounts(dds, gene="ENSMUSG00000020009", intgroup="Group", 
                returnData=TRUE)
ggplot(d, aes(x=Group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks = scales::log_breaks(n = 6, base = 10)) + 
  xlab("Genotype and Condition") + 
  ylab("Normalised counts") + 
  ggtitle("Ifngr1 - ENSMUSG00000020009 counts") + 
  theme(plot.title = element_text(margin = margin(t = 10, b = 10)))
ggsave("Normalised_counts_KO_Ifngr1.pdf")
d <- plotCounts(dds, gene="ENSMUSG00000022965", intgroup="Group", 
                returnData=TRUE)
ggplot(d, aes(x=Group, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks = scales::log_breaks(n = 6, base = 10)) + 
  xlab("Genotype and Condition") + 
  ylab("Normalised counts") + 
  ggtitle("Ifngr2 - ENSMUSG00000020009 counts") + 
  theme(plot.title = element_text(margin = margin(t = 10, b = 10)))
ggsave("Normalised_counts_KO_Ifngr2.pdf")


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
setwd(dir_step6_LFC_1_padj_001)

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
   


# ----------------------------------------------------------------------------------------------
# Volcano plot for the DKO infected vs WT infected
# ----------------------------------------------------------------------------------------------

# Add the column gene_id to res_DKO_Case_vs_WT_Case 
res_DKO_Case_vs_WT_Case <- add_col_gene_id(res_DKO_Case_vs_WT_Case)
# Add a column express_level to specify if they are UP- or DOWN- regulated with a adjusted p-value = padj_limit
res_DKO_Case_vs_WT_Case <- add_col_express_level(res_DKO_Case_vs_WT_Case)


# --- both up-regulated and down-regulated 

# Add the column gene_label 
res_DKO_Case_vs_WT_Case_up_down_regulated <- add_col_gene_label_up_down_regulated(res_DKO_Case_vs_WT_Case)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_DKO_Case_vs_WT_Case_up_down_regulated <- res_DKO_Case_vs_WT_Case_up_down_regulated[which(res_DKO_Case_vs_WT_Case_up_down_regulated$gene_label != ""),]


volcano_ggplot(res_DKO_Case_vs_WT_Case_up_down_regulated, 
               'Volcano plots for DKO infected vs WT infected up/down-regulated', 
               'volcano_plot_DKO_Case_vs_WT_Case_up_down.pdf', 
               int_genes_DKO_Case_vs_WT_Case_up_down_regulated)



# --- only up-regulated 

# Add the column gene_label 
res_DKO_Case_vs_WT_Case_up_regulated <- add_col_gene_label_up_regulated(res_DKO_Case_vs_WT_Case)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_DKO_Case_vs_WT_Case_up_regulated <- res_DKO_Case_vs_WT_Case_up_regulated[which(res_DKO_Case_vs_WT_Case_up_regulated$gene_label != ""),]

volcano_ggplot(res_DKO_Case_vs_WT_Case_up_regulated, 
               'Volcano plots for DKO infected vs WT infected up-regulated', 
               'volcano_plot_DKO_Case_vs_WT_Case_up.pdf', 
               int_genes_DKO_Case_vs_WT_Case_up_regulated)

# --- only down-regulated 

# Add the column gene_label 
res_DKO_Case_vs_WT_Case_down_regulated <- add_col_gene_label_down_regulated(res_DKO_Case_vs_WT_Case)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_DKO_Case_vs_WT_Case_down_regulated <- res_DKO_Case_vs_WT_Case_down_regulated[which(res_DKO_Case_vs_WT_Case_down_regulated$gene_label != ""),]

volcano_ggplot(res_DKO_Case_vs_WT_Case_down_regulated, 
               'Volcano plots for DKO infected vs WT infected down-regulated', 
               'volcano_plot_DKO_Case_vs_WT_Case_down.pdf', 
               int_genes_DKO_Case_vs_WT_Case_down_regulated)



# ----------------------------------------------------------------------------------------------
# Volcano plot for the WT infected vs WT not infected
# ----------------------------------------------------------------------------------------------

# Add the column gene_id to res_WT_Case_vs_WT_Control
res_WT_Case_vs_WT_Control <- add_col_gene_id(res_WT_Case_vs_WT_Control)

# Add a column express_level to specify if they are UP- or DOWN- regulated with adjusted padj_limit  
res_WT_Case_vs_WT_Control <- add_col_express_level(res_WT_Case_vs_WT_Control)

# --- both up-regulated and down-regulated 

# Add the column gene_label
res_WT_Case_vs_WT_Control_up_down_regulated <- add_col_gene_label_up_down_regulated(res_WT_Case_vs_WT_Control)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_WT_Case_vs_WT_Control_up_down_regulated <- res_WT_Case_vs_WT_Control_up_down_regulated[which(res_WT_Case_vs_WT_Control_up_down_regulated$gene_label != ""),]

volcano_ggplot(res_WT_Case_vs_WT_Control_up_down_regulated, 
               'Volcano plots for WT infected vs WT not infected', 
               'volcano_plot_WT_Case_vs_WT_Control_up_down.pdf', 
               int_genes_WT_Case_vs_WT_Control_up_down_regulated)


# --- only up-regulated 

# Add the column gene_label
res_WT_Case_vs_WT_Control_up_regulated <- add_col_gene_label_up_regulated(res_WT_Case_vs_WT_Control)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_WT_Case_vs_WT_Control_up_regulated <- res_WT_Case_vs_WT_Control_up_regulated[which(res_WT_Case_vs_WT_Control_up_regulated$gene_label != ""),]

volcano_ggplot(res_WT_Case_vs_WT_Control_up_regulated, 
               'Volcano plots for WT infected vs WT not infected up-regulated', 
               'volcano_plot_WT_Case_vs_WT_Control_up.pdf', 
               int_genes_WT_Case_vs_WT_Control_up_regulated)


# --- only down-regulated 

# Add the column gene_label
res_WT_Case_vs_WT_Control_down_regulated <- add_col_gene_label_down_regulated(res_WT_Case_vs_WT_Control)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_WT_Case_vs_WT_Control_down_regulated <- res_WT_Case_vs_WT_Control_down_regulated[which(res_WT_Case_vs_WT_Control_down_regulated$gene_label != ""),]

volcano_ggplot(res_WT_Case_vs_WT_Control_down_regulated, 
               'Volcano plots for WT infected vs WT not infected down-regulated', 
               'volcano_plot_WT_Case_vs_WT_Control_down.pdf', 
               int_genes_WT_Case_vs_WT_Control_down_regulated)


# ----------------------------------------------------------------------------------------------
# Volcano plot for the DKO infected vs DKO not infected
# ---------------------------------------------------------------------------------------------

# Add the column gene_id to res_DKO_Case_vs_DKO_Control
res_DKO_Case_vs_DKO_Control <- add_col_gene_id(res_DKO_Case_vs_DKO_Control)

# Add a column express_level to specify if they are UP- or DOWN- regulated with adjusted padj_limit
res_DKO_Case_vs_DKO_Control <- add_col_express_level(res_DKO_Case_vs_DKO_Control)


# --- both up-regulated and down-regulated 

# Add the column gene_label 
res_DKO_Case_vs_DKO_Control_up_down_regulated <- add_col_gene_label_up_down_regulated(res_DKO_Case_vs_DKO_Control)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_DKO_Case_vs_DKO_Control_up_down_regulated <- res_DKO_Case_vs_DKO_Control_up_down_regulated[which(res_DKO_Case_vs_DKO_Control_up_down_regulated$gene_label != ""),]

volcano_ggplot(res_DKO_Case_vs_DKO_Control_up_down_regulated, 
               'Volcano plots for DKO infected vs DKO not infected up/down-regulated', 
               'volcano_plot_DKO_Case_vs_DKO_Control_up_down.pdf', 
               int_genes_DKO_Case_vs_DKO_Control_up_down_regulated)


# --- only up-regulated 

# Add the column gene_label 
res_DKO_Case_vs_DKO_Control_up_regulated <- add_col_gene_label_up_regulated(res_DKO_Case_vs_DKO_Control)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_DKO_Case_vs_DKO_Control_up_regulated <- res_DKO_Case_vs_DKO_Control_up_regulated[which(res_DKO_Case_vs_DKO_Control_up_regulated$gene_label != ""),]

volcano_ggplot(res_DKO_Case_vs_DKO_Control_up_regulated, 
               'Volcano plots for DKO infected vs DKO not infected up-regulated', 
               'volcano_plot_DKO_Case_vs_DKO_Control_up.pdf', 
               int_genes_DKO_Case_vs_DKO_Control_up_regulated)

# --- only down-regulated 

# Add the column gene_label 
res_DKO_Case_vs_DKO_Control_down_regulated <- add_col_gene_label_down_regulated(res_DKO_Case_vs_DKO_Control)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_DKO_Case_vs_DKO_Control_down_regulated <- res_DKO_Case_vs_DKO_Control_down_regulated[which(res_DKO_Case_vs_DKO_Control_down_regulated$gene_label != ""),]

volcano_ggplot(res_DKO_Case_vs_DKO_Control_down_regulated, 
               'Volcano plots for DKO infected vs DKO not infected down-regulated', 
               'volcano_plot_DKO_Case_vs_DKO_Control_down.pdf', 
               int_genes_DKO_Case_vs_DKO_Control_down_regulated)


# ----------------------------------------------------------------------------------------------
# Volcano plot for the DKO not infected vs WT not infected
# ----------------------------------------------------------------------------------------------

# Add the column gene_id to res_DKO_Control_vs_WT_Control
res_DKO_Control_vs_WT_Control <- add_col_gene_id(res_DKO_Control_vs_WT_Control)

# Add a column express_level to specify if they are UP- or DOWN- regulated with an adjusted padj_limit
res_DKO_Control_vs_WT_Control <- add_col_express_level(res_DKO_Control_vs_WT_Control)


# --- both up-regulated and down-regulated 

# Add the column gene_label
res_DKO_Control_vs_WT_Control_up_down_regulated <- add_col_gene_label_up_down_regulated(res_DKO_Control_vs_WT_Control)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_DKO_Control_vs_WT_Control_up_down_regulated <- res_DKO_Control_vs_WT_Control_up_down_regulated[which(res_DKO_Control_vs_WT_Control_up_down_regulated$gene_label != ""),]

volcano_ggplot(res_DKO_Control_vs_WT_Control_up_down_regulated, 
               'Volcano plots for DKO not infected vs WT not infected up/down-regulated', 
               'volcano_plot_DKO_Control_vs_WT_Control_up_down.pdf',
               int_genes_DKO_Control_vs_WT_Control_up_down_regulated)


# --- only up-regulated 

# Add the column gene_label
res_DKO_Control_vs_WT_Control_up_regulated <- add_col_gene_label_up_regulated(res_DKO_Control_vs_WT_Control)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_DKO_Control_vs_WT_Control_up_regulated <- res_DKO_Control_vs_WT_Control_up_regulated[which(res_DKO_Control_vs_WT_Control_up_regulated$gene_label != ""),]

volcano_ggplot(res_DKO_Control_vs_WT_Control_up_regulated, 
               'Volcano plots for DKO not infected vs WT not infected up-regulated', 
               'volcano_plot_DKO_Control_vs_WT_Control_up.pdf', 
               int_genes_DKO_Control_vs_WT_Control_up_regulated)


# --- only down-regulated 

# Add the column gene_label
res_DKO_Control_vs_WT_Control_down_regulated <- add_col_gene_label_down_regulated(res_DKO_Control_vs_WT_Control)

# the ENSEMBL gene IDs that were filtered with the LFC_limit and padj_limit 
int_genes_DKO_Control_vs_WT_Control_down_regulated <- res_DKO_Control_vs_WT_Control_down_regulated[which(res_DKO_Control_vs_WT_Control_down_regulated$gene_label != ""),]

volcano_ggplot(res_DKO_Control_vs_WT_Control_down_regulated, 
               'Volcano plots for DKO not infected vs WT not infected down-regulated', 
               'volcano_plot_DKO_Control_vs_WT_Control_down.pdf', 
               int_genes_DKO_Control_vs_WT_Control_down_regulated)


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

dir_step7_LFC_1_padj_005 <- paste0(dir_step7, "/LFC_1_padj_0-05")
dir_step7_LFC_1_padj_001 <- paste0(dir_step7, "/LFC_1_padj_0-01")


# identification of the biological function of the possibile interesting genes 
# analysis upregulated and downregulated genes at the same time 

# --- in the contrast DKO infected vs WT infected
enrichGO_DKO_Case_vs_WT_Case <- enrichGO_plotting(int_genes_DKO_Case_vs_WT_Case_up_down_regulated, 
                                                  res_DKO_Case_vs_WT_Case_up_down_regulated, 
                                                  "DKO infected vs WT infected up and down regulated", 
                                                  dir_step7_LFC_1_padj_001)
enrichGO_DKO_Case_vs_WT_Case <- enrichGO_plotting(int_genes_DKO_Case_vs_WT_Case_up_regulated, 
                                                  res_DKO_Case_vs_WT_Case_up_regulated, 
                                                  "DKO infected vs WT infected up regulated", 
                                                  dir_step7_LFC_1_padj_001)
enrichGO_DKO_Case_vs_WT_Case <- enrichGO_plotting(int_genes_DKO_Case_vs_WT_Case_down_regulated, 
                                                  res_DKO_Case_vs_WT_Case_down_regulated, 
                                                  "DKO infected vs WT infected down regulated", 
                                                  dir_step7_LFC_1_padj_001)

# --- in the contrast WT infected vs WT not infected
enrichGO_WT_Case_vs_WT_Control <- enrichGO_plotting(int_genes_WT_Case_vs_WT_Control_up_down_regulated, 
                                                  res_WT_Case_vs_WT_Control_up_down_regulated, 
                                                  "WT infected vs WT not infected up and down regulated", 
                                                  dir_step7_LFC_1_padj_001)
enrichGO_WT_Case_vs_WT_Control <- enrichGO_plotting(int_genes_WT_Case_vs_WT_Control_up_regulated, 
                                                    res_WT_Case_vs_WT_Control_up_regulated, 
                                                    "WT infected vs WT not infected up regulated", 
                                                    dir_step7_LFC_1_padj_001)
enrichGO_WT_Case_vs_WT_Control <- enrichGO_plotting(int_genes_WT_Case_vs_WT_Control_down_regulated, 
                                                    res_WT_Case_vs_WT_Control_down_regulated, 
                                                    "WT infected vs WT not infected down regulated", 
                                                    dir_step7_LFC_1_padj_001)

# --- in the contrast DKO infected vs DKO not infected
enrichGO_DKO_Case_vs_DKO_Control <- enrichGO_plotting(int_genes_DKO_Case_vs_DKO_Control_up_down_regulated, 
                                                      res_DKO_Case_vs_DKO_Control_up_down_regulated, 
                                                    "DKO infected vs DKO not infected up and down regulated", 
                                                    dir_step7_LFC_1_padj_001)
enrichGO_DKO_Case_vs_DKO_Control <- enrichGO_plotting(int_genes_DKO_Case_vs_DKO_Control_up_regulated, 
                                                      res_DKO_Case_vs_DKO_Control_up_regulated, 
                                                      "DKO infected vs DKO not infected up regulated", 
                                                      dir_step7_LFC_1_padj_001)
enrichGO_DKO_Case_vs_DKO_Control <- enrichGO_plotting(int_genes_DKO_Case_vs_DKO_Control_down_regulated, 
                                                      res_DKO_Case_vs_DKO_Control_down_regulated, 
                                                      "DKO infected vs DKO not infected down regulated", 
                                                      dir_step7_LFC_1_padj_001)

# --- in the contrast DKO not infected vs WT not infected
enrichGO_DKO_Control_vs_WT_Control <- enrichGO_plotting(int_genes_DKO_Control_vs_WT_Control_up_down_regulated, 
                                                        res_DKO_Control_vs_WT_Control_up_down_regulated, 
                                                        "DKO not infected vs WT not infected up and down regulated", 
                                                        dir_step7_LFC_1_padj_001)
enrichGO_DKO_Control_vs_WT_Control <- enrichGO_plotting(int_genes_DKO_Control_vs_WT_Control_up_regulated, 
                                                        res_DKO_Control_vs_WT_Control_up_regulated, 
                                                        "DKO not infected vs WT not infected up regulated", 
                                                        dir_step7_LFC_1_padj_001)
enrichGO_DKO_Control_vs_WT_Control <- enrichGO_plotting(int_genes_DKO_Control_vs_WT_Control_down_regulated, 
                                                        res_DKO_Control_vs_WT_Control_down_regulated, 
                                                        "DKO not infected vs WT not infected down regulated", 
                                                        dir_step7_LFC_1_padj_001)











