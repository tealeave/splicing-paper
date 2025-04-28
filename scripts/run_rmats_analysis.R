#!/usr/bin/env Rscript

#' rMATS/MASER Alternative Splicing Analysis
#' 
#' This script runs the rMATS/MASER analysis pipeline for alternative splicing data.
#'
#' @author David Lin
#' @date April 28, 2025

# Load required libraries
suppressPackageStartupMessages({
  library(maser)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Source configuration and utility functions
source("config/analysis_config.R")
source("R/rmats_functions.R")

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
rmats_dir <- CONFIG$rmats_dir
comparisons <- CONFIG$rmats$comparisons
event_types <- CONFIG$rmats$event_types
figures_dir <- file.path(CONFIG$figures_dir, "splicing")
tables_dir <- file.path(CONFIG$tables_dir, "splicing")

# Create subdirectories
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}
if (!dir.exists(tables_dir)) {
  dir.create(tables_dir, recursive = TRUE)
}

# Log analysis start
cat("Starting rMATS/MASER analysis...\n")
cat("rMATS directory:", rmats_dir, "\n")
cat("Number of comparisons:", length(comparisons), "\n")
cat("Event types:", paste(event_types, collapse = ", "), "\n")

# Run rMATS/MASER pipeline
results <- tryCatch({
  run_rmats_pipeline(rmats_dir, comparisons, event_types, figures_dir)
}, error = function(e) {
  cat("Error in rMATS/MASER analysis:", e$message, "\n")
  return(NULL)
})

# Check if analysis was successful
if (!is.null(results) && length(results) > 0) {
  cat("rMATS/MASER analysis completed successfully.\n")
  cat("Results saved to:", tables_dir, "\n")
  cat("Figures saved to:", figures_dir, "\n")
} else {
  cat("rMATS/MASER analysis failed or no results were generated.\n")
}