#!/usr/bin/env Rscript

#' RNA Splicing Analysis Pipeline
#' 
#' This script runs the complete RNA splicing analysis pipeline,
#' including differential expression analysis with DESeq2 and
#' alternative splicing analysis with rMATS/MASER.
#'
#' @author David Lin
#' @date April, 2025

# Record start time
start_time <- Sys.time()

# Load required libraries
suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
})

# Source configuration
source("config/analysis_config.R")

# Create log directory
log_dir <- file.path(CONFIG$output_dir, "logs")
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}

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
  # Create Rmarkdown report
  rmarkdown::render(
    "reports/integrated_report.Rmd",
    output_file = file.path(CONFIG$output_dir, "integrated_report.html"),
    params = list(
      config = CONFIG,
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