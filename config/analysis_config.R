# RNA Splicing Analysis Configuration
# Created: April 28, 2025

# Base output directory
BASE_OUTPUT_DIR <- "./output"

# Define paths based on the base directory
CONFIG <- list(
  # Input data
  counts_dir = "./Counts",
  rmats_dir = "./RMATS_results",
  
  # Base Output directories
  output_dir = BASE_OUTPUT_DIR,
  logs_dir = file.path(BASE_OUTPUT_DIR, "logs"),
  reports_dir = file.path(BASE_OUTPUT_DIR, "reports"),
  
  # DESeq2 specific paths
  deseq2_output_dir = file.path(BASE_OUTPUT_DIR, "deseq2"),
  deseq2_figures_dir = file.path(BASE_OUTPUT_DIR, "deseq2", "figures"),
  deseq2_tables_dir = file.path(BASE_OUTPUT_DIR, "deseq2", "tables"),
  
  # rMATS/MASER specific paths
  rmats_maser_output_dir = file.path(BASE_OUTPUT_DIR, "rmats_maser"),
  rmats_maser_figures_dir = file.path(BASE_OUTPUT_DIR, "rmats_maser", "figures", "volcano_plots"),
  rmats_maser_summary_tables_dir = file.path(BASE_OUTPUT_DIR, "rmats_maser", "tables", "event_summaries"),
  rmats_maser_event_tables_dir = file.path(BASE_OUTPUT_DIR, "rmats_maser", "tables", "significant_events"),
  
  # Analysis parameters
  deseq2 = list(
    design = c("3", "3", "3", "3"),  # Number of replicates per condition
    conditions = list(
      MB468_000 = "MB468_000",
      MB468_30 = "MB468_30",
      MB468_120 = "MB468_120",
      MB468_720 = "MB468_720"
    ),
    comparisons = list(
      c("MB468_30", "MB468_000"),
      c("MB468_120", "MB468_000"),
      c("MB468_720", "MB468_000")
    ),
    output_files = c(
      "ComMB_0_30.csv",
      "ComMB_0_120.csv",
      "ComMB_0_720.csv"
    ),
    tables_dir   = file.path(BASE_OUTPUT_DIR, "deseq2", "tables"),
    counts_file = "./Counts/MBsimple_counts.txt",
    fdr_cutoff = 0.05
  ),
  
  # RMATS analysis parameters
  rmats = list(
    event_types = c("SE", "RI", "A3SS", "A5SS", "MXE"),
    comparisons = list(
      "468_0_v_468_12",
      "468_0_v_R8_0",
      "468_12_v_R8_12",
      "R8_0_v_R8_12"
    ),
    # Add filtering parameters used in rmats_functions.R
    fdr_cutoff = 0.05,
    deltaPSI_cutoff = 0.1,
    avg_reads_filter = 5
  )
)