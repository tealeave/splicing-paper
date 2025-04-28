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
  
  # Create MASER object
  maser_obj <- maser::maser(
    path = file.path(rmats_dir, comparison),
    cond1 = strsplit(comparison, "_v_")[[1]][1],
    cond2 = strsplit(comparison, "_v_")[[1]][2],
    ftype = "JCEC"
  )
  
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
  
  # Filter events
  filtered <- maser::filterByCoverage(maser_obj, event_type)
  filtered <- maser::filterByPval(filtered, event_type, fdr = fdr_cutoff)
  filtered <- maser::filterByDeltaPSI(filtered, event_type, dpsi = dpsi_cutoff)
  
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
  
  # Create plot
  p <- maser::plotVolcano(maser_obj, event_type, 
                         fdr = 0.05, dpsi = 0.1, 
                         xlim = c(-1, 1), ylim = c(0, 10))
  
  # Add title if provided
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }
  
  # Display plot
  print(p)
  
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