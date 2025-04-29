#' rMATS Analysis Functions
#' 
#' This file contains utility functions for alternative splicing analysis
#' using rMATS and MASER.
#'
#' @author David Lin
#' @date April 28, 2025

# Ensure necessary packages are loaded (can be done in the main script too)
# suppressPackageStartupMessages({
#   library(maser)
#   library(ggplot2)
# })

#' Create MASER object from rMATS results
#'
#' @param rmats_dir Directory containing rMATS results
#' @param comparison Comparison name (e.g., "468_0_v_468_12")
#' @param event_types Vector of event types to include (currently unused by maser constructor)
#' @return MASER object or NULL if creation fails
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
  
  # Extract condition names from comparison
  conds <- strsplit(comparison, "_v_")[[1]]
  if (length(conds) != 2) {
    message("Invalid comparison format: ", comparison)
    return(NULL)
  }
  cond1_label <- conds[1]
  cond2_label <- conds[2]
  
  # Create MASER object using the constructor from Maser_tests.R
  maser_obj <- tryCatch({
    # Use the constructor format from Maser_tests.R
    obj <- maser::maser(path, c(cond1_label, cond2_label), ftype = "JCEC")
    message("Successfully created MASER object for ", comparison)
    return(obj)
  }, error = function(e) {
    message("Error creating MASER object for ", comparison, ": ", e$message)
    # Fallback attempt removed for simplicity, ensure path/format is correct upstream
    return(NULL)
  })
  
  return(maser_obj)
}

#' Create volcano plot for alternative splicing events using maser::volcano
#'
#' @param maser_obj MASER object (should contain data for the event type)
#' @param event_type Event type (e.g., "SE")
#' @param fdr FDR cutoff (default: 0.05)
#' @param deltaPSI Delta PSI cutoff (default: 0.1, corresponds to IncLevelDifference)
#' @param title Plot title
#' @param output_dir Directory to save the plot
#' @param filename Name for the output PDF file
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @export
create_volcano_plot <- function(maser_obj, event_type, fdr = 0.05, deltaPSI = 0.1, title = NULL, 
                               output_dir, filename, width = 7, height = 7) { # Adjusted default size
  # This function requires the MASER package
  if (!requireNamespace("maser", quietly = TRUE)) {
    stop("MASER package is required for this function")
  }
  # Check if the MASER object is valid
  if (is.null(maser_obj) || !inherits(maser_obj, "Maser")) {
      warning("Invalid MASER object provided for plotting.")
      return(NULL)
  }

  # Check if the event type exists as a slot or if data is present
  # maser::volcano handles missing types internally, but we can add a check
  event_summary <- tryCatch({ summary(maser_obj, type = event_type) }, error = function(e) NULL)
  if (is.null(event_summary) || nrow(event_summary) == 0) {
      message("    No data found for event type ", event_type, " in the provided maser object. Skipping plot.")
      return(NULL)
  }

  # Ensure output directory exists (though typically done by caller)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, filename)
  
  # Create volcano plot using maser::volcano
  plot_obj <- tryCatch({
    # Use the official volcano function
    p <- maser::volcano(maser_obj, type = event_type, fdr = fdr, deltaPSI = deltaPSI)
    
    if (!is.null(p)) {
        # Apply theme customizations
        p <- p +
          ggplot2::theme(
            axis.title = ggplot2::element_text(size = 14),
            axis.text = ggplot2::element_text(size = 12),
            legend.key = ggplot2::element_rect(fill = 'bisque'),
            legend.text = ggplot2::element_text(size = 12),
            legend.title = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5)
          ) +
          ggplot2::xlim(-1.0, 1.0) # Consistent x-axis limits

        # Add title if provided
        if (!is.null(title)) {
          p <- p + ggplot2::ggtitle(title)
        }
        
        # Save using ggsave
        ggplot2::ggsave(filename = output_file, plot = p, device = "pdf", width = width, height = height)
        message("    Saved volcano plot to: ", output_file)
        return(p) # Return the plot object
    } else {
        message("    maser::volcano returned NULL for ", event_type, ". Plot not saved.")
        return(NULL)
    }
    
  }, error = function(e) {
    message("    Error creating volcano plot with maser::volcano for ", event_type, ": ", e$message)
    # Optionally create a placeholder plot if needed, but returning NULL is cleaner
    # plot(1, 1, type = "n", main = paste("Error plotting", event_type)); text(1, 1, "Plotting failed")
    return(NULL)
  })
  
  return(plot_obj)
}

#' Run complete rMATS/MASER analysis pipeline using official functions
#'
#' @param rmats_dir Directory containing rMATS results
#' @param comparisons List of comparisons to analyze
#' @param event_types Vector of event types to include
#' @param plots_output_dir Directory to save volcano plots (passed from caller)
#' @param tables_output_dir Directory to save summary tables and RDS objects (passed from caller)
#' @param fdr_cutoff FDR cutoff for significance (default: 0.05)
#' @param deltaPSI_cutoff Delta PSI cutoff for significance (default: 0.1)
#' @param avg_reads_filter Minimum average reads for filterByCoverage (default: 5)
#' @return List of MASER objects containing top significant events for each comparison
#' @export
run_rmats_pipeline <- function(rmats_dir, comparisons, event_types = c("SE", "RI", "A3SS", "A5SS", "MXE"), 
                              plots_output_dir, tables_output_dir, # Use dirs passed from caller
                              fdr_cutoff = 0.05, deltaPSI_cutoff = 0.1, avg_reads_filter = 5) {
  
  # Output directories are created by the calling script (run_rmats_analysis.R)
  
  # Run analysis for each comparison
  results_list <- list()
  
  for (comparison in comparisons) {
    message("Processing comparison: ", comparison)
    
    # 1. Create MASER object
    maser_obj_orig <- create_maser_object(rmats_dir, comparison, event_types)
    
    if (is.null(maser_obj_orig)) {
      message("  Skipping comparison due to MASER object creation failure: ", comparison)
      next
    }
    
    # 2. Filter by Coverage
    maser_obj_filt <- tryCatch({
       message("  Applying filterByCoverage with avg_reads = ", avg_reads_filter)
       obj <- maser::filterByCoverage(maser_obj_orig, avg_reads = avg_reads_filter)
       message("  Events after filtering:")
       print(obj) # Print summary after filtering
       obj
    }, error = function(e){
        message("  filterByCoverage failed for ", comparison, ": ", e$message, ". Using original object.")
        maser_obj_orig # Use original if filtering fails
    })
    
    # 3. Get Top Significant Events
    top_events_obj <- tryCatch({
        message("  Getting topEvents with fdr < ", fdr_cutoff, " and deltaPSI > ", deltaPSI_cutoff)
        obj <- maser::topEvents(maser_obj_filt, fdr = fdr_cutoff, deltaPSI = deltaPSI_cutoff)
        message("  Top significant events identified:")
        print(obj) # Print summary of top events
        obj
    }, error = function(e){
        message("  topEvents failed for ", comparison, ": ", e$message)
        NULL
    })
    
    # Store the top_events_obj (a Maser object) in the results list
    results_list[[comparison]] <- top_events_obj 
    
    # 4. Save Top Events Object
    if (!is.null(top_events_obj)) {
        rds_filename <- paste0(comparison, "_top_events.rds")
        rds_filepath <- file.path(tables_output_dir, rds_filename)
        tryCatch({
            saveRDS(top_events_obj, file = rds_filepath)
            message("  Saved top events MASER object to: ", rds_filepath)
        }, error = function(e) {
            message("  Error saving top events RDS for comparison ", comparison, ": ", e$message)
        })
    } else {
        message("  No top events object generated to save.")
    }

    # 5. Generate and Save Summary Table (based on filtered and top events)
    summary_table <- data.frame(
      EventType = character(),
      TotalEventsFiltered = integer(),
      SignificantEvents = integer(),
      UpregulatedEvents = integer(), 
      DownregulatedEvents = integer(),
      stringsAsFactors = FALSE
    )
    
    # Get summaries from the filtered object and the top events object
    summary_filt_list <- tryCatch({ summary(maser_obj_filt) }, error = function(e) { message(" Error getting summary for filtered obj: ", e$message); NULL })
    summary_top_list <- if (!is.null(top_events_obj)) tryCatch({ summary(top_events_obj) }, error = function(e) { message(" Error getting summary for top obj: ", e$message); NULL }) else NULL

    processed_any_type = FALSE
    for (event_type in event_types) {
        message("  Summarizing event type: ", event_type)
        
        total_filt_count <- 0
        if (!is.null(summary_filt_list) && event_type %in% names(summary_filt_list)) {
            total_filt_count <- nrow(summary_filt_list[[event_type]])
        }
        
        sig_count <- 0
        up_count <- 0
        down_count <- 0
        if (!is.null(summary_top_list) && event_type %in% names(summary_top_list)) {
            top_df <- summary_top_list[[event_type]]
            if (!is.null(top_df) && nrow(top_df) > 0 && "IncLevelDifference" %in% colnames(top_df)) {
                sig_count <- nrow(top_df)
                # Ensure numeric for comparison
                inc_level_diff_num <- suppressWarnings(as.numeric(as.character(top_df$IncLevelDifference)))
                up_count <- sum(inc_level_diff_num > deltaPSI_cutoff, na.rm = TRUE)
                down_count <- sum(inc_level_diff_num < -deltaPSI_cutoff, na.rm = TRUE)
            }
        }
        
        message("    Total (filtered): ", total_filt_count, ", Significant: ", sig_count, ", Up: ", up_count, ", Down: ", down_count)

        # Add row only if there were filtered events for this type
        if (total_filt_count > 0 || sig_count > 0) {
             summary_table <- rbind(summary_table, data.frame(
                EventType = event_type,
                TotalEventsFiltered = total_filt_count,
                SignificantEvents = sig_count,
                UpregulatedEvents = up_count,
                DownregulatedEvents = down_count,
                stringsAsFactors = FALSE
            ))
            processed_any_type = TRUE
        }
    }
    
    # Add placeholder if no types had data
    if (!processed_any_type) {
         summary_table <- rbind(summary_table, data.frame(
            EventType = "None", TotalEventsFiltered = 0, SignificantEvents = 0, UpregulatedEvents = 0, DownregulatedEvents = 0, stringsAsFactors = FALSE
         ))
         message("  No processable event types found for summary in comparison: ", comparison)
    }

    # Save summary table
    summary_filename <- paste0(comparison, "_summary.csv")
    summary_filepath <- file.path(tables_output_dir, summary_filename)
    tryCatch({
      write.csv(summary_table, file = summary_filepath, row.names = FALSE)
      message("  Saved summary table to: ", summary_filepath)
    }, error = function(e) {
      message("  Error saving summary table for comparison ", comparison, ": ", e$message)
    })

    # 6. Generate and Save Volcano Plots (using filtered object)
    message("  Generating volcano plots...")
    for (event_type in event_types) {
      plot_title <- paste(comparison, event_type, "(FDR<",fdr_cutoff, ", |dPSI|>", deltaPSI_cutoff, ")")
      plot_filename <- paste0(comparison, "_", event_type, "_volcano.pdf")
      # Pass the filtered object for plotting all points, highlighting significant ones
      create_volcano_plot(maser_obj_filt, event_type, fdr = fdr_cutoff, deltaPSI = deltaPSI_cutoff, 
                          title = plot_title, output_dir = plots_output_dir, filename = plot_filename)
    }
  }
  
  message("Pipeline finished for all comparisons.")
  return(results_list) # Return the list of top_events Maser objects
}