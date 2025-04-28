#' rMATS Analysis Functions
#' 
#' This file contains utility functions for alternative splicing analysis
#' using rMATS and MASER.
#'
#' @author David Lin
#' @date April, 2025

#' Load rMATS results for a specific event type
#'
#' @param rmats_dir Directory containing rMATS results
#' @param comparison Comparison name (e.g., "468_0_v_468_12")
#' @param event_type Event type (e.g., "SE", "RI", "A3SS", "A5SS", "MXE")
#' @param jc_or_jcec Whether to use JC or JCEC results (default: "JCEC")
#' @return Data frame with rMATS results
#' @export
load_rmats_results <- function(rmats_dir, comparison, event_type, jc_or_jcec = "JCEC") {
  file_path <- file.path(rmats_dir, comparison, paste0(event_type, ".MATS.", jc_or_jcec, ".txt"))
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  results <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  return(results)
}

#' Create MASER object from rMATS results
#'
#' @param rmats_dir Directory containing rMATS results
#' @param comparison Comparison name
#' @param event_types Vector of event types to include
#' @return MASER object
#' @export
create_maser_object <- function(rmats_dir, comparison, event_types = c("SE", "RI", "A3SS", "A5SS", "MXE")) {
  # This function requires the MASER package
  if (!requireNamespace("maser", quietly = TRUE)) {
    stop("MASER package is required for this function")
  }
  
  # Get the full path to the rMATS results
  path <- file.path(rmats_dir, comparison)
  
  # Check if the directory exists
  if (!dir.exists(path)) {
    warning(paste("Directory does not exist:", path))
    return(NULL)
  }
  
  # Check if the required files exist
  files_exist <- TRUE
  for (event_type in event_types) {
    file_path <- file.path(path, paste0(event_type, ".MATS.JCEC.txt"))
    if (!file.exists(file_path)) {
      warning(paste("File does not exist:", file_path))
      files_exist <- FALSE
    }
  }
  
  if (!files_exist) {
    return(NULL)
  }
  
  # Create MASER object with basic constructor
  maser_obj <- tryCatch({
    obj <- maser::maser(path)
    return(obj)
  }, error = function(e) {
    warning(paste("Error creating MASER object:", e$message))
    return(NULL)
  })
  
  return(maser_obj)
}

#' Get significant events from MASER object
#'
#' @param maser_obj MASER object
#' @param event_type Event type
#' @param fdr FDR cutoff (default: 0.05)
#' @param dpsi Delta PSI cutoff (default: 0.1)
#' @return Filtered MASER object with significant events
#' @export
getSig <- function(maser_obj, event_type, fdr = 0.05, dpsi = 0.1) {
  # Check if the MASER object is NULL
  if (is.null(maser_obj)) {
    warning("MASER object is NULL")
    return(NULL)
  }
  
  # Check if the event type exists in the MASER object
  event_exists <- tryCatch({
    # Try to access the event type slot
    slot(maser_obj, event_type)
    TRUE
  }, error = function(e) {
    warning(paste("Event type", event_type, "not found in MASER object"))
    FALSE
  })
  
  if (!event_exists) {
    return(NULL)
  }
  
  # Get the event data
  event_data <- slot(maser_obj, event_type)
  
  # Check if the event data has metadata
  if (!("metadata" %in% slotNames(event_data))) {
    warning(paste("Event type", event_type, "has no metadata"))
    return(NULL)
  }
  
  # Check if the metadata has FDR and dPSI columns
  metadata <- event_data@metadata
  if (!all(c("FDR", "IncLevelDifference") %in% colnames(metadata))) {
    warning(paste("Event type", event_type, "metadata missing required columns"))
    return(NULL)
  }
  
  # Filter significant events
  sig_idx <- metadata$FDR < fdr & abs(metadata$IncLevelDifference) > dpsi
  
  # If no significant events, return NULL
  if (!any(sig_idx)) {
    warning(paste("No significant events found for", event_type))
    return(NULL)
  }
  
  # Create a new object with only significant events
  sig_obj <- maser_obj
  sig_metadata <- metadata[sig_idx, ]
  slot(slot(sig_obj, event_type), "metadata") <- sig_metadata
  
  return(sig_obj)
}

#' Run complete rMATS/MASER analysis pipeline
#'
#' @param rmats_dir Directory containing rMATS results
#' @param comparisons List of comparisons to analyze
#' @param event_types Vector of event types to include
#' @param output_dir Directory to save output figures (tables go to ../../tables/splicing)
#' @return List of MASER objects
#' @export
run_rmats_pipeline <- function(rmats_dir, comparisons, event_types = c("SE", "RI", "A3SS", "A5SS", "MXE"), 
                              output_dir = "./output/figures") {
  # Define base output directories
  base_output_dir <- dirname(output_dir)
  figures_output_dir <- file.path(base_output_dir, "figures", "splicing")
  tables_output_dir <- file.path(base_output_dir, "tables", "splicing")
  
  # Create output directories if they don't exist
  dir.create(figures_output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tables_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Run analysis for each comparison
  results_list <- list()
  
  for (comparison in comparisons) {
    cat("Processing comparison:", comparison, "\n")
    
    # Create MASER object
    maser_obj <- create_maser_object(rmats_dir, comparison, event_types)
    
    if (is.null(maser_obj)) {
      warning(paste("Failed to create MASER object for comparison:", comparison))
      next
    }
    
    # Store results
    results_list[[comparison]] <- maser_obj
    
    # Create summary data
    summary_table <- data.frame(
      EventType = character(),
      TotalEvents = integer(),
      SignificantEvents = integer(),
      UpregulatedEvents = integer(),
      DownregulatedEvents = integer(),
      stringsAsFactors = FALSE
    )
    
    # Flag to track if we found any events for this comparison
    found_events_for_comparison <- FALSE
    
    # Process each event type
    for (event_type in event_types) {
      cat("  Processing event type:", event_type, "\n")
      
      # Check if the event type exists as a slot
      if (!(event_type %in% slotNames(maser_obj))) {
        warning(paste("Event type slot", event_type, "not found in MASER object for comparison:", comparison))
        next
      }
      
      # Get event data
      event_data_slot <- tryCatch({
        slot(maser_obj, event_type)
      }, error = function(e) {
        warning(paste("Error accessing slot", event_type, "for comparison", comparison, ":", e$message))
        return(NULL)
      })
      
      if (is.null(event_data_slot)) {
        next
      }
      
      # Check if metadata exists within the event slot
      if (!("metadata" %in% slotNames(event_data_slot))) {
        warning(paste("Event type", event_type, "has no metadata slot for comparison:", comparison))
        next
      }
      
      metadata <- event_data_slot@metadata
      
      # Check if metadata is empty
      if (is.null(metadata) || nrow(metadata) == 0) {
        cat("    No metadata found for event type", event_type, "\n")
        next
      }
      
      # Check if required columns exist (FDR and IncLevelDifference)
      required_cols <- c("FDR", "IncLevelDifference")
      if (!all(required_cols %in% colnames(metadata))) {
        warning(paste("Event type", event_type, "metadata missing required columns (FDR, IncLevelDifference) for comparison:", comparison))
        next
      }
      
      # Mark that we found events for this comparison
      found_events_for_comparison <- TRUE
      
      # Count events
      total_events <- nrow(metadata)
      # Ensure FDR and IncLevelDifference are numeric and handle NAs
      metadata$FDR <- as.numeric(metadata$FDR)
      metadata$IncLevelDifference <- as.numeric(metadata$IncLevelDifference)
      
      sig_idx <- !is.na(metadata$FDR) & metadata$FDR < 0.05 & !is.na(metadata$IncLevelDifference) & abs(metadata$IncLevelDifference) > 0.1
      up_idx <- !is.na(metadata$FDR) & metadata$FDR < 0.05 & !is.na(metadata$IncLevelDifference) & metadata$IncLevelDifference > 0.1
      down_idx <- !is.na(metadata$FDR) & metadata$FDR < 0.05 & !is.na(metadata$IncLevelDifference) & metadata$IncLevelDifference < -0.1
      
      sig_events <- sum(sig_idx)
      up_events <- sum(up_idx)
      down_events <- sum(down_idx)
      
      cat("    Total events:", total_events, "Significant:", sig_events, "\n")
      
      # Add to summary table
      summary_table <- rbind(summary_table, data.frame(
        EventType = event_type,
        TotalEvents = total_events,
        SignificantEvents = sig_events,
        UpregulatedEvents = up_events,
        DownregulatedEvents = down_events,
        stringsAsFactors = FALSE
      ))
      
      # Create volcano plot
      tryCatch({
        # Create output file path
        plot_file <- file.path(figures_output_dir, paste0(comparison, "_", event_type, "_volcano.pdf"))
        
        # Create plot
        pdf(plot_file, width = 8, height = 6)
        
        # Plot the data manually
        plot_title <- paste(comparison, event_type, "Volcano Plot")
        # Handle potential -Inf values from log10(0)
        log_fdr <- -log10(metadata$FDR)
        log_fdr[is.infinite(log_fdr)] <- max(log_fdr[is.finite(log_fdr)], 0) + 1 # Replace Inf with a value slightly larger than max finite
        
        plot(metadata$IncLevelDifference, log_fdr,
             xlab = "Delta PSI (Inclusion Level Difference)", ylab = "-log10(FDR)",
             main = plot_title,
             pch = 16, col = "gray",
             xlim = c(-1, 1), ylim = c(0, max(log_fdr[is.finite(log_fdr)], 1)) # Adjust ylim based on data
            )
        
        # Highlight significant events
        if (any(sig_idx)) {
          points(metadata$IncLevelDifference[sig_idx], 
                 log_fdr[sig_idx],
                 pch = 16, col = "red")
        }
        
        # Add horizontal line at FDR = 0.05
        abline(h = -log10(0.05), lty = 2, col = "blue")
        
        # Add vertical lines at dPSI = ±0.1
        abline(v = c(-0.1, 0.1), lty = 2, col = "blue")
        
        dev.off()
        
        cat("    Created volcano plot for", event_type, "\n")
      }, error = function(e) {
        warning(paste("Error creating volcano plot for", event_type, "in comparison", comparison, ":", e$message))
        # Ensure device is closed if error occurs during plotting
        if (names(dev.cur()) == "pdf") dev.off()
      })
      
      # Save event data
      tryCatch({
        write.csv(metadata, 
                  file = file.path(tables_output_dir, 
                                  paste0(comparison, "_", event_type, "_events.csv")),
                  row.names = FALSE)
        cat("    Saved event data for", event_type, "\n")
      }, error = function(e) {
        warning(paste("Error saving event data for", event_type, "in comparison", comparison, ":", e$message))
      })
    }
    
    # If no events were found for the entire comparison, add a placeholder row
    if (!found_events_for_comparison) {
      summary_table <- rbind(summary_table, data.frame(
        EventType = "None",
        TotalEvents = 0,
        SignificantEvents = 0,
        UpregulatedEvents = 0,
        DownregulatedEvents = 0,
        stringsAsFactors = FALSE
      ))
      cat("No processable events found for comparison:", comparison, "\n")
    }
    
    # Save summary table
    tryCatch({
      write.csv(summary_table, 
                file = file.path(tables_output_dir, 
                                paste0(comparison, "_summary.csv")),
                row.names = FALSE)
      cat("Saved summary table for comparison", comparison, "\n")
    }, error = function(e) {
      warning(paste("Error saving summary table for comparison", comparison, ":", e$message))
    })
  }
  
  return(results_list)
}