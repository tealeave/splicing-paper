#!/usr/bin/env Rscript

#' RNA Splicing Analysis Pipeline
#' 
#' This script runs the complete RNA splicing analysis pipeline,
#' including differential expression analysis with DESeq2 and
#' alternative splicing analysis with rMATS/MASER.
#'
#' @author David Lin
#' @date April 28, 2025

# Record start time
start_time <- Sys.time()

# Load required libraries
suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
})

# Source configuration
source("config/analysis_config.R")

# Create all necessary directories
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(CONFIG$figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(CONFIG$tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$figures_dir, "splicing"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$tables_dir, "splicing"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$output_dir, "reports"), showWarnings = FALSE, recursive = TRUE)

# Create log directory
log_dir <- file.path(CONFIG$output_dir, "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# Set up logging
log_file <- file.path(log_dir, paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
con <- file(log_file, open = "w")
sink(con, type = "output")
sink(con, type = "message")

# Print header
cat("=======================================================\n")
cat("RNA Splicing Analysis Pipeline\n")
cat("Started:", format(start_time), "\n")
cat("=======================================================\n\n")

# Define pipeline steps
cat("Step 1: Running DESeq2 differential expression analysis...\n")
tryCatch({
  source("scripts/run_deseq2_analysis.R")
  cat("DESeq2 analysis completed successfully.\n\n")
}, error = function(e) {
  cat("Error in DESeq2 analysis:", e$message, "\n\n")
})

cat("Step 2: Running rMATS/MASER alternative splicing analysis...\n")
tryCatch({
  source("scripts/run_rmats_analysis.R")
  cat("rMATS/MASER analysis completed successfully.\n\n")
}, error = function(e) {
  cat("Error in rMATS/MASER analysis:", e$message, "\n\n")
})

# Generate integrated report
cat("Step 3: Generating integrated report...\n")
tryCatch({
  # Create reports directory using absolute path
  reports_dir <- normalizePath(file.path(CONFIG$output_dir, "reports"), mustWork = FALSE)
  dir.create(reports_dir, showWarnings = TRUE, recursive = TRUE)
  
  # Create the directory manually if it doesn't exist
  if (!dir.exists(reports_dir)) {
    system(paste("mkdir -p", reports_dir))
    Sys.sleep(1)  # Wait a moment to ensure directory is created
  }
  
  # Double-check that the directory exists
  if (!dir.exists(reports_dir)) {
    stop(paste("Failed to create directory:", reports_dir))
  } else {
    cat("Reports directory created successfully:", reports_dir, "\n")
  }
  
  # Convert relative paths in CONFIG to absolute paths
  project_dir <- normalizePath(".")
  abs_config <- CONFIG
  
  # Convert relative paths to absolute
  if (grepl("^\\./", abs_config$counts_dir)) {
    abs_config$counts_dir <- file.path(project_dir, substring(abs_config$counts_dir, 3))
  }
  if (grepl("^\\./", abs_config$rmats_dir)) {
    abs_config$rmats_dir <- file.path(project_dir, substring(abs_config$rmats_dir, 3))
  }
  if (grepl("^\\./", abs_config$output_dir)) {
    abs_config$output_dir <- file.path(project_dir, substring(abs_config$output_dir, 3))
  }
  if (grepl("^\\./", abs_config$figures_dir)) {
    abs_config$figures_dir <- file.path(project_dir, substring(abs_config$figures_dir, 3))
  }
  if (grepl("^\\./", abs_config$tables_dir)) {
    abs_config$tables_dir <- file.path(project_dir, substring(abs_config$tables_dir, 3))
  }
  
  # Create Rmarkdown report
  rmarkdown::render(
    input = normalizePath("reports/integrated_report.Rmd"),
    output_file = file.path(reports_dir, "integrated_report.html"),
    params = list(
      config = abs_config,
      date = Sys.Date()
    )
  )
  cat("Report generation completed successfully.\n\n")
}, error = function(e) {
  cat("Error in report generation:", e$message, "\n\n")
})

# Record end time and calculate duration
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

# Print footer
cat("=======================================================\n")
cat("Pipeline completed\n")
cat("Started:", format(start_time), "\n")
cat("Ended:", format(end_time), "\n")
cat("Duration:", round(as.numeric(duration), 2), "minutes\n")
cat("=======================================================\n")

# Close log file
sink(type = "output")
sink(type = "message")
close(con)

# Print final message to console
cat("Pipeline completed in", round(as.numeric(duration), 2), "minutes\n")
cat("Log file:", log_file, "\n")
cat("Results directory:", CONFIG$output_dir, "\n")