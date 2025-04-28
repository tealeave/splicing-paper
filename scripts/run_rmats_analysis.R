#!/usr/bin/env Rscript

#' rMATS/MASER Alternative Splicing Analysis
#' 
#' This script runs the rMATS/MASER analysis pipeline for alternative splicing data.
#'
#' @author David Lin
#' @date April, 2025

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

# If analysis was successful, create summary tables
if (!is.null(results)) {
  # Create summary tables for each comparison
  for (comparison in names(results)) {
    maser_obj <- results[[comparison]]
    
    # Create summary table for all event types
    summary_table <- data.frame(
      EventType = character(),
      TotalEvents = integer(),
      SignificantEvents = integer(),
      UpregulatedEvents = integer(),
      DownregulatedEvents = integer(),
      stringsAsFactors = FALSE
    )
    
    for (event_type in event_types) {
      # Get events
      events <- tryCatch({
        getSig(maser_obj, event_type, fdr = 0.05, dpsi = 0.1)
      }, error = function(e) {
        warning(paste("Error getting significant events for", event_type, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(events)) {
        # Try to access metadata
        metadata <- tryCatch({
          slot(events, event_type)@metadata
        }, error = function(e) {
          warning(paste("Error accessing metadata for", event_type, ":", e$message))
          return(NULL)
        })
        
        if (!is.null(metadata) && nrow(metadata) > 0) {
          # Count events
          total_events <- nrow(metadata)
          
          # Check if FDR and dPSI columns exist
          if (all(c("FDR", "dPSI") %in% colnames(metadata))) {
            sig_events <- sum(metadata$FDR < 0.05 & abs(metadata$dPSI) > 0.1)
            up_events <- sum(metadata$FDR < 0.05 & metadata$dPSI > 0.1)
            down_events <- sum(metadata$FDR < 0.05 & metadata$dPSI < -0.1)
          } else {
            sig_events <- NA
            up_events <- NA
            down_events <- NA
          }
          
          # Add to summary table
          summary_table <- rbind(summary_table, data.frame(
            EventType = event_type,
            TotalEvents = total_events,
            SignificantEvents = sig_events,
            UpregulatedEvents = up_events,
            DownregulatedEvents = down_events,
            stringsAsFactors = FALSE
          ))
          
          # Save detailed table
          write.csv(metadata,
                    file = file.path(tables_dir, paste0(comparison, "_", event_type, "_events.csv")),
                    row.names = FALSE)
        }
      }
    }
    
    # Save summary table
    write.csv(summary_table, 
              file = file.path(tables_dir, paste0(comparison, "_summary.csv")),
              row.names = FALSE)
  }
  
  cat("rMATS/MASER analysis completed successfully.\n")
  cat("Results saved to:", tables_dir, "\n")
  cat("Figures saved to:", figures_dir, "\n")
} else {
  cat("rMATS/MASER analysis failed.\n")
}