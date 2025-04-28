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
  
  # Extract conditions from comparison name
  comparison_parts <- strsplit(comparison, "_v_")[[1]]
  
  # Create MASER object - using the path only, not the additional parameters
  # that were causing the error
  maser_obj <- tryCatch({
    maser::maser(path = file.path(rmats_dir, comparison))
  }, error = function(e) {
    # If the standard constructor fails, try with the older version syntax
    maser::maser(
      path = file.path(rmats_dir, comparison),
      ftype = "JCEC"
    )
  })
  
  return(maser_obj)
}

#' Filter significant events from MASER object
#'
#' @param maser_obj MASER object
#' @param event_type Event type
#' @param fdr_cutoff FDR cutoff (default: 0.05)
#' @param dpsi_cutoff Delta PSI cutoff (default: 0.1)
#' @return Filtered MASER object
#' @export
filter_significant_events <- function(maser_obj, event_type, fdr_cutoff = 0.05, dpsi_cutoff = 0.1) {
  # This function requires the MASER package
  if (!requireNamespace("maser", quietly = TRUE)) {
    stop("MASER package is required for this function")
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
    return(maser_obj)
  }
  
  # Try to filter events with error handling
  filtered <- tryCatch({
    # First try to filter by coverage
    filtered_coverage <- tryCatch({
      if (exists("filterByCoverage", where = asNamespace("maser"), mode = "function")) {
        maser::filterByCoverage(maser_obj, event_type)
      } else {
        maser_obj
      }
    }, error = function(e) {
      warning(paste("Error filtering by coverage:", e$message))
      return(maser_obj)
    })
    
    # Then try to filter by p-value
    filtered_pval <- tryCatch({
      if (exists("filterByPval", where = asNamespace("maser"), mode = "function")) {
        maser::filterByPval(filtered_coverage, event_type, fdr = fdr_cutoff)
      } else {
        # Manual filtering by FDR
        event_data <- slot(filtered_coverage, event_type)
        if ("metadata" %in% slotNames(event_data) && "FDR" %in% colnames(event_data@metadata)) {
          sig_idx <- event_data@metadata$FDR < fdr_cutoff
          if (any(sig_idx)) {
            event_data@metadata <- event_data@metadata[sig_idx, ]
          }
          slot(filtered_coverage, event_type) <- event_data
        }
        filtered_coverage
      }
    }, error = function(e) {
      warning(paste("Error filtering by p-value:", e$message))
      return(filtered_coverage)
    })
    
    # Finally try to filter by delta PSI
    filtered_dpsi <- tryCatch({
      if (exists("filterByDeltaPSI", where = asNamespace("maser"), mode = "function")) {
        maser::filterByDeltaPSI(filtered_pval, event_type, dpsi = dpsi_cutoff)
      } else {
        # Manual filtering by delta PSI
        event_data <- slot(filtered_pval, event_type)
        if ("metadata" %in% slotNames(event_data) && "dPSI" %in% colnames(event_data@metadata)) {
          sig_idx <- abs(event_data@metadata$dPSI) > dpsi_cutoff
          if (any(sig_idx)) {
            event_data@metadata <- event_data@metadata[sig_idx, ]
          }
          slot(filtered_pval, event_type) <- event_data
        }
        filtered_pval
      }
    }, error = function(e) {
      warning(paste("Error filtering by delta PSI:", e$message))
      return(filtered_pval)
    })
    
    return(filtered_dpsi)
  }, error = function(e) {
    warning(paste("Error filtering events:", e$message))
    return(maser_obj)
  })
  
  return(filtered)
}

#' Create volcano plot for alternative splicing events
#'
#' @param maser_obj MASER object
#' @param event_type Event type
#' @param title Plot title
#' @param output_file Output file path (if NULL, plot is displayed but not saved)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @export
create_volcano_plot <- function(maser_obj, event_type, title = NULL, 
                               output_file = NULL, width = 8, height = 6) {
  # This function requires the MASER package
  if (!requireNamespace("maser", quietly = TRUE)) {
    stop("MASER package is required for this function")
  }
  
  # Create volcano plot
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
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
    # Create an empty plot with a message if the event type doesn't exist
    plot(1, 1, type = "n", xlab = "Delta PSI", ylab = "-log10(FDR)", 
         main = paste("No", event_type, "events found"))
    text(1, 1, "No events found", cex = 1.5)
    
    if (!is.null(output_file)) {
      dev.off()
    }
    return(NULL)
  }
  
  # Try to create a volcano plot using different methods
  tryCatch({
    # First try using plotVolcano if it exists
    if (exists("plotVolcano", where = asNamespace("maser"), mode = "function")) {
      p <- maser::plotVolcano(maser_obj, event_type, 
                             fdr = 0.05, dpsi = 0.1, 
                             xlim = c(-1, 1), ylim = c(0, 10))
    } else {
      # If plotVolcano doesn't exist, try using volcano
      if (exists("volcano", where = asNamespace("maser"), mode = "function")) {
        p <- maser::volcano(maser_obj, event_type, 
                           fdr = 0.05, dpsi = 0.1)
      } else {
        # If neither function exists, create a basic volcano plot manually
        # Get event data
        event_data <- slot(maser_obj, event_type)@metadata
        
        # Create a basic volcano plot
        plot(event_data$dPSI, -log10(event_data$FDR), 
             xlab = "Delta PSI", ylab = "-log10(FDR)",
             main = paste(event_type, "Volcano Plot"),
             pch = 16, col = "gray")
        
        # Highlight significant events
        sig_idx <- event_data$FDR < 0.05 & abs(event_data$dPSI) > 0.1
        if (any(sig_idx)) {
          points(event_data$dPSI[sig_idx], -log10(event_data$FDR[sig_idx]), 
                 pch = 16, col = "red")
        }
        
        # Add a horizontal line at FDR = 0.05
        abline(h = -log10(0.05), lty = 2, col = "blue")
        
        # Add vertical lines at dPSI = ±0.1
        abline(v = c(-0.1, 0.1), lty = 2, col = "blue")
        
        p <- NULL  # No ggplot object to return
      }
    }
    
    # Add title if provided and p is a ggplot object
    if (!is.null(title) && !is.null(p) && "ggplot" %in% class(p)) {
      p <- p + ggplot2::ggtitle(title)
      print(p)
    }
    
  }, error = function(e) {
    # If all methods fail, create a simple error message plot
    plot(1, 1, type = "n", xlab = "Delta PSI", ylab = "-log10(FDR)", 
         main = paste("Error creating", event_type, "volcano plot"))
    text(1, 1, paste("Error:", e$message), cex = 1.2)
  })
  
  # Close device if saving
  if (!is.null(output_file)) {
    dev.off()
  }
  
  return(p)
}

#' Run complete rMATS/MASER analysis pipeline
#'
#' @param rmats_dir Directory containing rMATS results
#' @param comparisons List of comparisons to analyze
#' @param event_types Vector of event types to include
#' @param output_dir Directory to save output files
#' @return List of MASER objects
#' @export
run_rmats_pipeline <- function(rmats_dir, comparisons, event_types = c("SE", "RI", "A3SS", "A5SS", "MXE"), 
                              output_dir = "./output/figures") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Run analysis for each comparison
  results_list <- list()
  
  for (comparison in comparisons) {
    # Create MASER object
    maser_obj <- create_maser_object(rmats_dir, comparison, event_types)
    
    # Analyze each event type
    for (event_type in event_types) {
      # Filter significant events
      filtered <- filter_significant_events(maser_obj, event_type)
      
      # Create volcano plot
      plot_title <- paste(comparison, event_type, "Volcano Plot")
      output_file <- file.path(output_dir, paste0(comparison, "_", event_type, "_volcano.pdf"))
      create_volcano_plot(filtered, event_type, title = plot_title, output_file = output_file)
    }
    
    # Store results
    results_list[[comparison]] <- maser_obj
  }
  
  return(results_list)
}