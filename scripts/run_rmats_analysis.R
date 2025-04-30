#!/usr/bin/env Rscript

#' rMATS/MASER Alternative Splicing Analysis
#' 
#' This script runs the rMATS/MASER analysis pipeline for alternative splicing data.
#'
#' @author David Lin

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

# --- Setup Output Directories and Logging ---
BASE_OUTPUT_DIR <- CONFIG$rmats_maser_output_dir
LOGS_DIR <- file.path(BASE_OUTPUT_DIR, "logs")
PLOTS_DIR <- file.path(BASE_OUTPUT_DIR, "plots") # Consolidated figures dir
TABLES_DIR <- file.path(BASE_OUTPUT_DIR, "tables") # Consolidated tables dir

# Create output directories if they don't exist
dir.create(BASE_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOGS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

# Create a timestamped log file name
log_file_name <- paste0("run_rmats_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
log_file_path <- file.path(LOGS_DIR, log_file_name)

# Redirect console output (stdout and stderr/messages) to the log file
log_con <- file(log_file_path, open = "wt")
sink(log_con, type = "output") # Redirects print(), cat(), etc.
sink(log_con, type = "message") # Redirects message(), warning()

# Ensure sink is closed even if the script stops with an error
on.exit({
  message("Closing log file...")
  sink(type = "message")
  sink(type = "output")
  close(log_con)
  message("Log file closed: ", log_file_path)
})

message("Logging started. Output directed to: ", log_file_path)
message("Base output directory: ", BASE_OUTPUT_DIR)

# --- Parameters --- 
# Set up parameters from config
rmats_dir <- CONFIG$rmats_dir
comparisons <- CONFIG$rmats$comparisons
event_types <- CONFIG$rmats$event_types
fdr_cutoff <- CONFIG$rmats$fdr_cutoff
deltaPSI_cutoff <- CONFIG$rmats$deltaPSI_cutoff
avg_reads_filter <- CONFIG$rmats$avg_reads_filter

# Log analysis start
message("Starting rMATS/MASER analysis...")
message("rMATS directory: ", rmats_dir)
message("Number of comparisons: ", length(comparisons))
message("Event types: ", paste(event_types, collapse = ", "))
message("FDR cutoff: ", fdr_cutoff)
message("Delta PSI cutoff: ", deltaPSI_cutoff)
message("Avg Reads Filter: ", avg_reads_filter)

# --- Run Pipeline --- 
# Run rMATS/MASER pipeline
results <- tryCatch({
  run_rmats_pipeline(rmats_dir = rmats_dir, 
                     comparisons = comparisons, 
                     event_types = event_types, 
                     plots_output_dir = PLOTS_DIR,      # Use consolidated dir
                     tables_output_dir = TABLES_DIR,    # Use consolidated dir
                     fdr_cutoff = fdr_cutoff, 
                     deltaPSI_cutoff = deltaPSI_cutoff, 
                     avg_reads_filter = avg_reads_filter)
}, error = function(e) {
  message("Error in rMATS/MASER analysis pipeline: ", e$message)
  # Print stack trace if possible (might depend on R environment options)
  # print(sys.calls())
  return(NULL)
})

# --- Wrap Up --- 
# Check if analysis was successful
if (!is.null(results) && length(results) > 0) {
  message("rMATS/MASER analysis completed successfully.")
  message("Summary and event tables saved to: ", TABLES_DIR)
  message("Plots saved to: ", PLOTS_DIR)
} else {
  message("rMATS/MASER analysis failed or no results were generated.")
}

message("Script finished.")
# The on.exit() function will handle closing the sink automatically here