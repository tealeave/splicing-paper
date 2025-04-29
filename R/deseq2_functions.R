#' DESeq2 Analysis Functions
#' 
#' This file contains utility functions for RNA-seq differential expression analysis
#' using DESeq2.
#'
#' @author David Lin
#' @date April, 2025

#' Load count data from file
#'
#' @param file_path Path to the count data file
#' @param sep Separator used in the file (default: tab)
#' @return A data frame containing the count data
#' @export
load_count_data <- function(file_path, sep = "\t") {
  counts <- read.table(file_path, header = TRUE, sep = sep, row.names = 1)
  
  # Convert real numbers to integers (DESeq2 requirement)
  numeric_idx <- sapply(counts, mode) == 'numeric'
  counts[numeric_idx] <- round(counts[numeric_idx], 0)
  
  return(counts)
}

#' Create DESeq2 dataset from count data
#'
#' @param counts Count data matrix
#' @param condition_groups List of condition groups
#' @param reps Number of replicates per condition
#' @return DESeq2 dataset
#' @export
create_deseq_dataset <- function(counts, condition_groups, reps) {
  # Create condition factor
  condition_values <- unlist(
    mapply(
      function(cond, rep_count) rep(cond, rep_count),
      condition_groups,
      reps,
      SIMPLIFY = FALSE
    )
  )
  
  # Ensure we have the right number of samples
  samples <- colnames(counts)
  if (length(samples) != length(condition_values)) {
    stop("Number of samples in count data (", length(samples), 
         ") does not match number of samples in condition groups (", 
         length(condition_values), ")")
  }
  
  # Create colData with proper sample names
  condition <- factor(condition_values)
  col_data <- data.frame(row.names = samples)
  col_data$condition <- condition
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = col_data,
    design = ~condition
  )
  
  return(dds)
}

#' Run DESeq2 analysis and extract results
#'
#' @param dds DESeq2 dataset
#' @param comparison Vector with two condition names to compare
#' @return Data frame with DESeq2 results
#' @export
run_deseq_analysis <- function(dds, comparison) {
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results for specific comparison
  res <- results(dds, contrast = c("condition", comparison[1], comparison[2]))
  
  # Convert to data frame
  data <- data.frame(res)
  
  # Rename columns
  names(data)[names(data) == "pvalue"] <- "PValue"
  names(data)[names(data) == "padj"] <- "FDR"
  
  # Add additional columns
  data$foldChange <- 2 ^ data$log2FoldChange
  data$PAdj <- p.adjust(data$PValue, method = "hochberg")
  
  # Sort by p-value
  data <- data[with(data, order(PValue, -foldChange)), ]
  
  # Compute false discovery counts
  data$falsePos <- 1:nrow(data) * data$FDR
  
  # Rearrange columns
  data <- data[c(1, 7, 2, 3, 4, 5, 8, 6, 9)]
  
  return(data)
}

#' Save DESeq2 results to CSV file
#'
#' @param data DESeq2 results data frame
#' @param output_dir Directory to save the file
#' @param filename Name of the output CSV file
#' @export
save_deseq_results <- function(data, output_dir, filename) {
  # Ensure output directory exists
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Construct full file path
  file_path <- file.path(output_dir, filename)
  
  # Save the data
  write.csv(data, file = file_path, row.names = TRUE, quote = FALSE)
  cat("  Saved DESeq2 results to:", file_path, "\n")
}

#' Run complete DESeq2 analysis pipeline
#'
#' @param count_file Path to count data file
#' @param condition_groups List of condition names
#' @param reps Number of replicates per condition
#' @param comparisons List of comparisons to make
#' @param output_filenames List of output filenames (not full paths)
#' @param output_dir Directory to save result tables
#' @return List of result data frames
#' @export
run_deseq_pipeline <- function(count_file, condition_groups, reps, comparisons, output_filenames, output_dir) {
  # Load count data
  counts <- load_count_data(count_file)
  
  # Create DESeq2 dataset
  dds <- create_deseq_dataset(counts, condition_groups, reps)
  
  # Run analysis for each comparison
  results_list <- list()
  
  for (i in seq_along(comparisons)) {
    cat("  Running DESeq2 comparison:", paste(comparisons[[i]], collapse=" vs "), "\n")
    # Run analysis
    results <- run_deseq_analysis(dds, comparisons[[i]])
    
    # Save results using the specified directory and filename
    save_deseq_results(results, output_dir, output_filenames[[i]])
    
    # Store results
    results_list[[i]] <- results
  }
  
  return(results_list)
}