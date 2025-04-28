#!/usr/bin/env Rscript

#' DESeq2 Differential Expression Analysis
#' 
#' This script runs the DESeq2 analysis pipeline for RNA-seq data.
#'
#' @author David Lin
#' @date April, 2025

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# Source configuration and utility functions
source("config/analysis_config.R")
source("R/deseq2_functions.R")

# Create output directories if they don't exist
if (!dir.exists(CONFIG$output_dir)) {
  dir.create(CONFIG$output_dir, recursive = TRUE)
}
if (!dir.exists(CONFIG$figures_dir)) {
  dir.create(CONFIG$figures_dir, recursive = TRUE)
}
if (!dir.exists(CONFIG$tables_dir)) {
  dir.create(CONFIG$tables_dir, recursive = TRUE)
}

# Set up parameters
count_file <- file.path(CONFIG$counts_dir, "MBsimple_counts.txt")
reps <- as.integer(CONFIG$deseq2$design)
condition_groups <- CONFIG$deseq2$conditions
comparisons <- CONFIG$deseq2$comparisons
output_files <- lapply(CONFIG$deseq2$output_files, function(file) {
  file.path(CONFIG$tables_dir, file)
})

# Log analysis start
cat("Starting DESeq2 analysis...\n")
cat("Count file:", count_file, "\n")
cat("Number of comparisons:", length(comparisons), "\n")

# Run DESeq2 pipeline
results <- tryCatch({
  run_deseq_pipeline(count_file, condition_groups, reps, comparisons, output_files)
}, error = function(e) {
  cat("Error in DESeq2 analysis:", e$message, "\n")
  return(NULL)
})

# If analysis was successful, create additional visualizations
if (!is.null(results)) {
  # Load count data
  counts <- load_count_data(count_file)
  
  # Create DESeq2 dataset
  dds <- create_deseq_dataset(counts, condition_groups, reps)
  dds <- DESeq(dds)
  
  # Create PCA plot
  vsd <- vst(dds, blind = FALSE)
  pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    theme_minimal() +
    ggtitle("PCA of RNA-seq Samples")
  
  # Save PCA plot
  ggsave(file.path(CONFIG$figures_dir, "pca_plot.pdf"), pca_plot, width = 8, height = 6)
  
  # Create heatmap of sample distances
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(vsd)
  colnames(sampleDistMatrix) <- colnames(vsd)
  
  # Save heatmap
  pdf(file.path(CONFIG$figures_dir, "sample_distance_heatmap.pdf"), width = 8, height = 6)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           main = "Sample Distance Heatmap")
  dev.off()
  
  # Create MA plots for each comparison
  for (i in seq_along(comparisons)) {
    # Convert data.frame to DESeqResults object for plotMA
    res_data <- results[[i]]
    res_obj <- DESeqResults(
      DataFrame(
        baseMean = res_data$baseMean,
        log2FoldChange = res_data$log2FoldChange,
        pvalue = res_data$PValue,
        padj = res_data$FDR
      )
    )
    
    # Create MA plot
    pdf(file.path(CONFIG$figures_dir, paste0("ma_plot_", 
                                           gsub(" ", "_", paste(comparisons[[i]], collapse = "_vs_")), 
                                           ".pdf")), 
        width = 8, height = 6)
    plotMA(res_obj, main = paste("MA Plot:", comparisons[[i]][1], "vs", comparisons[[i]][2]))
    dev.off()
  }
  
  cat("DESeq2 analysis completed successfully.\n")
  cat("Results saved to:", CONFIG$tables_dir, "\n")
  cat("Figures saved to:", CONFIG$figures_dir, "\n")
} else {
  cat("DESeq2 analysis failed.\n")
}