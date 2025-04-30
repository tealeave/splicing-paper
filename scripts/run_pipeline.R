#!/usr/bin/env Rscript

#' RNA Splicing Analysis Pipeline
#' 
#' This script runs the complete RNA splicing analysis pipeline,
#' including differential expression analysis with DESeq2 and
#' alternative splicing analysis with rMATS/MASER.
#'
#' @author David Lin

# Record start time
start_time <- Sys.time()

# Load required libraries (optional, as sub-scripts load their own)
# suppressPackageStartupMessages({
#   library(targets)
#   library(tarchetypes)
# })

# Source configuration
source("config/analysis_config.R")

# --- Directory Creation ---
# Create base output, logs, and reports directories from config
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(CONFIG$logs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(CONFIG$reports_dir, showWarnings = FALSE, recursive = TRUE)
# Note: Sub-analysis directories (deseq2/, rmats_maser/) are created within their respective scripts/functions

# --- Logging Setup ---
log_file <- file.path(CONFIG$logs_dir, paste0("pipeline_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
con <- file(log_file, open = "w")
sink(con, type = "output")
sink(con, type = "message")

# Print header
cat("=======================================================\n")
cat("RNA Splicing Analysis Pipeline\n")
cat("Started:", format(start_time), "\n")
cat("Using Configuration:", normalizePath("config/analysis_config.R"), "\n")
cat("Output Directory:", normalizePath(CONFIG$output_dir), "\n")
cat("=======================================================\n\n")

# --- Helper function to create absolute paths ---
#' Create Absolute Path
#'
#' Converts a relative path to an absolute path based on a base directory.
#' If the path is already absolute, it returns it unchanged.
#'
#' @param path The path to convert (character string).
#' @param base_dir The base directory to use if the path is relative (character string).
#' @return An absolute path (character string).
make_absolute <- function(path, base_dir) {
  # Check if path is already absolute
  # Simple check: starts with / or ~ on Unix-like, or drive letter on Windows
  is_absolute <- grepl("^(/|~|[A-Za-z]:)", path) 
  if (is_absolute) {
    return(normalizePath(path, mustWork = FALSE))
  } else {
    return(normalizePath(file.path(base_dir, path), mustWork = FALSE))
  }
}


# # --- Run DESeq2 Analysis ---
# cat("Step 1: Running DESeq2 differential expression analysis...\n")
# tryCatch({
#   source("scripts/run_deseq2_analysis.R")
#   # Logging is now done within run_deseq2_analysis.R
#   # cat("DESeq2 analysis completed successfully.\n\n")
# }, error = function(e) {
#   cat("Error in DESeq2 analysis step:", e$message, "\n\n")
# })

# --- Run rMATS/MASER Analysis ---
cat("Step 2: Running rMATS/MASER alternative splicing analysis...\n")
tryCatch({
  source("scripts/run_rmats_analysis.R")
  # Logging is now done within run_rmats_analysis.R
  # cat("rMATS/MASER analysis completed successfully.\n\n")
}, error = function(e) {
  cat("Error in rMATS/MASER analysis step:", e$message, "\n\n")
})

# --- Generate Integrated Report ---
cat("\nStep 3: Generating integrated report...\n")
tryCatch({
  # Ensure reports directory exists (already created above, but double-check)
  reports_dir_abs <- normalizePath(CONFIG$reports_dir, mustWork = FALSE)
  if (!dir.exists(reports_dir_abs)) {
     dir.create(reports_dir_abs, showWarnings = TRUE, recursive = TRUE)
     if (!dir.exists(reports_dir_abs)) stop(paste("Failed to create reports directory:", reports_dir_abs))
  }
  cat("Using reports directory:", reports_dir_abs, "\n")
  
  # Convert relative paths in CONFIG to absolute paths for the report
  project_dir <- normalizePath(".")
  abs_config <- list(
    counts_dir = make_absolute(CONFIG$counts_dir, project_dir),
    rmats_dir = make_absolute(CONFIG$rmats_dir, project_dir),
    output_dir = make_absolute(CONFIG$output_dir, project_dir),
    logs_dir = make_absolute(CONFIG$logs_dir, project_dir),
    reports_dir = make_absolute(CONFIG$reports_dir, project_dir),
    deseq2 = list(
      design        = CONFIG$deseq2$design,
      conditions    = CONFIG$deseq2$conditions,
      comparisons   = CONFIG$deseq2$comparisons,
      output_files  = CONFIG$deseq2$output_files,
      counts_file   = make_absolute(CONFIG$deseq2$counts_file, project_dir),
      tables_dir    = make_absolute(CONFIG$deseq2_tables_dir, project_dir),
      figures_dir   = make_absolute(CONFIG$deseq2_figures_dir, project_dir),
      fdr_cutoff    = CONFIG$deseq2$fdr_cutoff
    ),
    rmats = list(
      event_types          = CONFIG$rmats$event_types,
      comparisons          = CONFIG$rmats$comparisons,
      summary_tables_dir   = make_absolute(CONFIG$rmats_maser_summary_tables_dir, project_dir),
      figures_dir          = make_absolute(CONFIG$rmats_maser_figures_dir, project_dir),
      event_tables_dir     = make_absolute(CONFIG$rmats_maser_event_tables_dir, project_dir),
      fdr_cutoff           = CONFIG$rmats$fdr_cutoff,
      deltaPSI_cutoff      = CONFIG$rmats$deltaPSI_cutoff,
      avg_reads_filter     = CONFIG$rmats$avg_reads_filter
    )
  )
  
  # Define input and output paths for rendering
  report_input_path <- normalizePath("reports/integrated_report.Rmd")
  report_output_path <- file.path(reports_dir_abs, "integrated_report.html")
  
  cat("Rendering report:", report_input_path, "\n")
  cat("Outputting report to:", report_output_path, "\n")
  
  # Create Rmarkdown report
  rmarkdown::render(
    input = report_input_path,
    output_file = report_output_path,
    params = list(
      config = abs_config, # Pass the config with absolute paths
      date = Sys.Date()
    ),
    envir = new.env() # Render in a clean environment
  )
  cat("Report generation completed successfully.\n\n")
}, error = function(e) {
  cat("Error in report generation:", e$message, "\n\n")
})

# --- Pipeline Completion ---
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

# Print footer
cat("=======================================================\n")
cat("Pipeline completed\n")
cat("Started:", format(start_time), "\n")
cat("Ended:", format(end_time), "\n")
cat("Duration:", round(as.numeric(duration), 2), "minutes\n")
cat("Output located in:", normalizePath(CONFIG$output_dir), "\n")
cat("=======================================================\n")

# Close log file
sink(type = "output")
sink(type = "message")
close(con)

# Print final message to console
cat("Pipeline completed in", round(as.numeric(duration), 2), "minutes\n")
cat("Log file:", log_file, "\n")
cat("Results directory:", normalizePath(CONFIG$output_dir), "\n")