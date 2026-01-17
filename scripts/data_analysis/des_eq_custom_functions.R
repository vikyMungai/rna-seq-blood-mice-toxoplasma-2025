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

plot_counts_single_gene<- function(dds, gene_name, gene_id, file_tag){
  #   Plot the normalized count for a specific gene comparing the counts 
  #   in  the two possibile genotypes (WT and DKO) and the two treatments 
  #   (Control and Infected) 
  #
  #   Parameter: 
  #     dds (DESeqDataSet): counts from FeatureCounts  
  #     gene_name: gene symbol from ENSEMBl associated to gene id 
  #     gene_id: ENSEMBL gene id 
  #     file_tag: end of the file name where to save the plot. 
  #     The name of the file will be fromatted: ""Normalised_counts_[file_tag].pdf""
  #     
  #
  #   Returns: 
  #     None 
  d <- plotCounts(dds, gene=gene_id, 
                  intgroup=c("Disease", "Genotype"), 
                  returnData=TRUE)
  ggplot(d, aes(x=Genotype, y=count, col = Disease)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) + 
    scale_y_log10(breaks = scales::log_breaks(n = 6, base = 10)) + 
    stat_summary(fun = mean, 
                 geom = "line", 
                 aes(group = Disease), 
                 colour = "black") + 
    scale_color_manual(values = c("steelblue3", "plum"), 
                       labels = c("Control", "Infected")) + 
    xlab("Genotype") + 
    ylab("Normalised counts") + 
    ggtitle(paste0(gene_name, " - ", gene_id, " counts")) + 
    theme(plot.title = element_text(margin = margin(t = 10, b = 10)))
  ggsave( paste0("Normalised_counts_", file_tag, ".pdf"))
}


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
  
  
  # only plotting a subset of the interesting genes 
  top_genes <- head(interesting_genes[order(interesting_genes$padj), "gene_id"], 40)
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
  
  ggplot(data = results_DESeq, aes(x = log2FoldChange, y = -log10(padj), col = express_level)) + 
    geom_vline(xintercept = c(-log2_fold_change, log2_fold_change), col = "gray", linetype = 'dashed') + 
    geom_hline(yintercept = -log10(padj_value), col = "gray", linetype = 'dashed') + 
    geom_point(size = 1) +
    scale_color_manual(values = c("steelblue3", "grey", "plum"), 
                       labels = c("Downregulated", "Not significant", "Upregulated")) + 
    coord_cartesian( xlim = c(-20, 15)) + 
    labs(color = 'Gene express level', x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
    ggtitle(plot_title) +
    # max overlaps was chosen Inf to not lose labels 
    geom_text_repel(aes(label = all_symbols), size = 3, max.overlaps = Inf, show.legend = FALSE)
  
  ggsave(output_file)
  
  return(
    top_genes
  )
  
}

# ----------------------------------------------------------------------------------------------
# functions used in step 7
# ----------------------------------------------------------------------------------------------

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
            orderBy = "p.adjust",
            title = paste0("Barplot for interesting genes in contrast ", contrast))
  )
  
  
  print("DB: printed barplot")
  
  #plotting with dot plot 
  print(
    enrichplot::dotplot(enrichGO_results, 
                        showCategory=10, 
                        x = 'GeneRatio', 
                        color = 'p.adjust', 
                        orderBy = 'p.adjust' ) +  
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


