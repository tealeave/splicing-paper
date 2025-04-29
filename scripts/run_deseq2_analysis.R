#!/usr/bin/env Rscript

#' DESeq2 Differential Expression Analysis
#' 
#' This script runs the DESeq2 analysis pipeline for RNA-seq data.
#'
#' @author David Lin
#' @date April 28, 2025

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
dir.create(CONFIG$deseq2_output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(CONFIG$deseq2_figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(CONFIG$deseq2_tables_dir, showWarnings = FALSE, recursive = TRUE)

# Set up parameters
count_file <- file.path(CONFIG$counts_dir, "MBsimple_counts.txt") # Assuming MBsimple for now
reps <- as.integer(CONFIG$deseq2$design)
condition_groups <- CONFIG$deseq2$conditions
comparisons <- CONFIG$deseq2$comparisons
output_filenames <- CONFIG$deseq2$output_files

# Log analysis start
cat("Starting DESeq2 analysis...\n")
cat("Count file:", count_file, "\n")
cat("Number of comparisons:", length(comparisons), "\n")

# Run DESeq2 pipeline
results <- tryCatch({
  # Pass the specific tables directory
  run_deseq_pipeline(count_file, condition_groups, reps, comparisons, output_filenames, CONFIG$deseq2_tables_dir)
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
  
  # Save PCA plot to the correct figures directory
  pca_plot_path <- file.path(CONFIG$deseq2_figures_dir, "pca_plot.pdf")
  ggsave(pca_plot_path, pca_plot, width = 8, height = 6)
  cat("  Saved PCA plot to:", pca_plot_path, "\n")
  
  # Create heatmap of sample distances
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(vsd)
  colnames(sampleDistMatrix) <- colnames(vsd)
  
  # Save heatmap to the correct figures directory
  heatmap_path <- file.path(CONFIG$deseq2_figures_dir, "sample_distance_heatmap.pdf")
  pdf(heatmap_path, width = 8, height = 6)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           main = "Sample Distance Heatmap")
  dev.off()
  cat("  Saved heatmap to:", heatmap_path, "\n")
  
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
    
    # Save MA plot to the correct figures directory
    ma_plot_filename <- paste0("ma_plot_", 
                             gsub(" ", "_", paste(comparisons[[i]], collapse = "_vs_")), 
                             ".pdf")
    ma_plot_path <- file.path(CONFIG$deseq2_figures_dir, ma_plot_filename)
    pdf(ma_plot_path, width = 8, height = 6)
    plotMA(res_obj, main = paste("MA Plot:", comparisons[[i]][1], "vs", comparisons[[i]][2]))
    dev.off()
    cat("  Saved MA plot to:", ma_plot_path, "\n")
  }
  
  cat("DESeq2 analysis completed successfully.\n")
  cat("Results saved to:", CONFIG$deseq2_tables_dir, "\n")
  cat("Figures saved to:", CONFIG$deseq2_figures_dir, "\n")
} else {
  cat("DESeq2 analysis failed.\n")
}