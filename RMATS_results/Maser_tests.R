# Apply R best practices and manage output directories

# --- Configuration ---
BASE_OUTPUT_DIR <- "./output/test"
PLOTS_DIR <- file.path(BASE_OUTPUT_DIR, "plots")
TABLES_DIR <- file.path(BASE_OUTPUT_DIR, "tables")
LOGS_DIR <- file.path(BASE_OUTPUT_DIR, "logs")

# Create output directories if they don't exist
dir.create(BASE_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOGS_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Setup Logging ---
# Create a timestamped log file name
log_file_name <- paste0("maser_test_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
log_file_path <- file.path(LOGS_DIR, log_file_name)

# Redirect console output (stdout and stderr/messages) to the log file
# Use append = FALSE to overwrite if run again quickly, TRUE to append
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

# --- Load Libraries ---
# Ensure necessary packages are installed before running
# BiocManager::install("maser")
# BiocManager::install("rtracklayer")
# install.packages("ggplot2") # Needed for ggsave and themes

library(maser)
library(rtracklayer)
library(ggplot2) # For saving plots and themes

# --- Parameters ---
# Input path for rMATS results
rmats_results_path <- "./RMATS_results/468_0_v_468_12"
# Condition names
condition_names <- c("MB468_0mins", "MB468_720mins")
# rMATS file type (JCEC or JC)
fmats_file_type <- "JCEC"
# Filtering parameters
avg_reads_threshold <- 5
fdr_threshold <- 0.05
delta_psi_threshold <- 0.1

# --- Analysis ---

# 1. Load rMATS data using maser
message("Loading rMATS data from: ", rmats_results_path)
maser_obj <- maser(rmats_results_path, condition_names, ftype = fmats_file_type)
print(maser_obj)

# 2. Summarize and save initial data
message("Summarizing initial data...")
initial_summary <- summary(maser_obj) # Get summary for all types
# Save summary tables for each splicing type
for (splicing_type in names(initial_summary)) {
  output_file <- file.path(TABLES_DIR, paste0("initial_summary_", splicing_type, ".tsv"))
  write.table(initial_summary[[splicing_type]], file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved initial summary for ", splicing_type, " to: ", output_file)
}
# Print head of SE summary as an example
print(head(summary(maser_obj, type = "SE")[, 1:8]))


# 3. Filter by coverage
message("Filtering events by average read coverage (threshold: ", avg_reads_threshold, ")")
maser_obj_filt <- filterByCoverage(maser_obj, avg_reads = avg_reads_threshold)
message("Events after filtering:")
print(maser_obj_filt)


# 4. Identify top significant events
message("Identifying top significant events (FDR < ", fdr_threshold, ", deltaPSI > ", delta_psi_threshold, ")")
maser_obj_top <- topEvents(maser_obj_filt, fdr = fdr_threshold, deltaPSI = delta_psi_threshold)
message("Top significant events:")
print(maser_obj_top)
# Save top events summary
top_events_summary_file <- file.path(TABLES_DIR, "top_events_summary.rds")
saveRDS(maser_obj_top, file = top_events_summary_file)
message("Saved top events maser object to: ", top_events_summary_file)
# Optionally save as table if needed for external tools
# top_events_df <- slot(maser_obj_top, "sig_events") # Access the data frame
# write.csv(top_events_df, file.path(TABLES_DIR, "top_events_summary.csv"), row.names = FALSE)


# 5. Generate and save Volcano plots for all event types (using filtered data)
splicing_event_types <- c("SE", "RI", "MXE", "A5SS", "A3SS")
message("Generating and saving volcano plots for filtered data...")

for (event_type in splicing_event_types) {
  message("  Processing: ", event_type)
  plot_title <- paste(condition_names[1], "vs", condition_names[2], "-", event_type)
  plot_filename <- file.path(PLOTS_DIR, paste0("volcano_", event_type, "_filtered.pdf"))

  # Generate plot object
  p <- volcano(maser_obj_filt, fdr = fdr_threshold, deltaPSI = delta_psi_threshold, type = event_type)

  # Check if plot is not NULL (volcano returns NULL if no data for type)
  if (!is.null(p)) {
    # Customize theme (optional, apply desired theme adjustments here)
    p <- p +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.key = element_rect(fill = 'bisque'),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), # Often preferred
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5) # Center title
      ) +
      ggtitle(plot_title) # Add title

    # Save the plot
    ggsave(filename = plot_filename, plot = p, device = "pdf", width = 7, height = 7)
    message("    Saved plot to: ", plot_filename)
  } else {
    message("    No data to plot for type: ", event_type)
  }
}


# 6. Generate and save Volcano plots for TOP significant events only
message("Generating and saving volcano plots for top significant events...")

for (event_type in splicing_event_types) {
    message("  Processing: ", event_type)
    plot_title <- paste("Top Significant:", condition_names[1], "vs", condition_names[2], "-", event_type)
    plot_filename <- file.path(PLOTS_DIR, paste0("volcano_", event_type, "_top_significant.pdf"))

    # Generate plot object using the topEvents object
    # The volcano function should handle cases where no events of this type exist in maser_obj_top
    p <- volcano(maser_obj_top, fdr = fdr_threshold, deltaPSI = delta_psi_threshold, type = event_type)

    # Check if plot object was successfully created
    if (!is.null(p)) {
        # Customize theme
        p <- p +
            theme(
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            legend.key = element_rect(fill = 'bisque'),
            legend.text = element_text(size = 12),
            legend.title = element_blank(),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
            ) +
            ggtitle(plot_title) +
            xlim(-1.0, 1.0) # Adjust xlim if needed, e.g., xlim(-0.75, 0.75)

        # Save the plot
        ggsave(filename = plot_filename, plot = p, device = "pdf", width = 7, height = 7)
        message("    Saved plot to: ", plot_filename)
    } else {
            message("    No significant events of type '", event_type, "' found in topEvents object or plot could not be generated.")
    }
}


# 7. Generate and save Splicing Distribution plot
message("Generating and saving splicing distribution plot...")
splicing_dist_plot_filename <- file.path(PLOTS_DIR, "splicing_distribution.pdf")
# Need to capture the plot object if possible, or use pdf()/dev.off()
pdf(splicing_dist_plot_filename, width = 7, height = 5)
splicingDistribution(maser_obj_filt, fdr = fdr_threshold, deltaPSI = delta_psi_threshold)
dev.off()
message("Saved splicing distribution plot to: ", splicing_dist_plot_filename)


# 8. Generate and save Dotplot for top SE events
message("Generating and saving dotplot for top 30 SE events...")
dotplot_se_filename <- file.path(PLOTS_DIR, "dotplot_SE_top30_significant.pdf")

# Check if there are significant SE events by checking the summary of the top object
top_se_plot <- dotplot(maser_obj_top, fdr = fdr_threshold, deltaPSI = delta_psi_threshold, type = "SE")

# Save the plot if it exists
if (!is.null(top_se_plot)) {
    ggsave(filename = dotplot_se_filename, plot = top_se_plot, device = "pdf", width = 7, height = 5)
    message("Saved dotplot for top SE events to: ", dotplot_se_filename)
} else {
    message("No significant SE events found in the topEvents object or plot could not be generated.")
}

