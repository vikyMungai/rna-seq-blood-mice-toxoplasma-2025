if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("ggplot2")
BiocManager::install("DESeq2")

# browseVignettes("DESeq2")

library(ggplot2)

library(DESeq2)

setwd("/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/results/feature_counts")

# read counts data 
counts_data <- read.csv("formatted_feature_count.txt", 
                        header = TRUE, 
                        sep = "\t", 
                        row.names="Geneid") 
# row.names="Geneid" is used to have as columns only the samples names 

new_col_names <- c("SRR7821949", "SRR7821950", "SRR7821951", "SRR7821952", 
                   "SRR7821953", "SRR7821954", "SRR7821955", "SRR7821956", "SRR7821957", 
                   "SRR7821968", "SRR7821969", "SRR7821970", "SRR7821971", "SRR7821972", 
                   "SRR7821973")
colnames(counts_data) <- new_col_names

# read the sample information table 
setwd("/Users/vittoriamungai/Desktop/RNA_sequencing/RNA_project/data/raw_data")
col_data <- read.csv("README.csv", 
                     header = TRUE, 
                     sep = "\t")

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

# converting column Group of col_data into factors 
col_data$Group <- factor(col_data$Group)

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = col_data,
                              design = ~ Group)
# set the factor level 
# we explicit say which level we want as reference
# otherwise it will choose alphabetically 
dds$Group <- relevel(dds$Group, ref = "Blood_WT_Control")

# If there are technical replicates (not biological replicates) then you should 
# collapse them before running the differential analysis  

# run DESeq
dds <- DESeq(dds)

# lists the coeffiecients
resultsNames(dds)


vsd <- vst(dds, blind = TRUE) # DESeqTransform
head(assay(vsd), 3)


pcaData <- plotPCA(vsd, intgroup = "Group", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Group)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


#------ see results from dds (not requested)
res <- results(dds) 
# res <- results(dds, contrast=c("condition","Blood_WT_Control","Blood_DKO_Control"))
# res <- results(dds, contrast=c("condition","Blood_WT_Case","Blood_DKO_Case"))
res

summary(res)

# see what was compared in the result 
res_WT_Control_vs_DKO_Control <- results(dds, contrast=c("condition","Blood_WT_Control","Blood_DKO_Control"))
res_WT_Control_vs_Blood_WT_Case <- results(dds, contrast=c("condition","Blood_WT_Control","Blood_WT_Case"))
res_WT_Control_vs_DKO_Case <- results(dds, contrast=c("condition","Blood_WT_Control","Blood_DKO_Case"))

plotMA(res)
# ----


