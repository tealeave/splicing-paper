# RNA Splicing Analysis Configuration
# Created: April 28, 2025

# Data paths
CONFIG <- list(
  # Input data
  counts_dir = "./Counts",
  rmats_dir = "./RMATS_results",
  
  # Output directories
  output_dir = "./output",
  figures_dir = "./output/figures",
  tables_dir = "./output/tables",
  
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
    output_files = list(
      "ComMB_0_30.csv",
      "ComMB_0_120.csv",
      "ComMB_0_720.csv"
    )
  ),
  
  # RMATS analysis parameters
  rmats = list(
    event_types = c("SE", "RI", "A3SS", "A5SS", "MXE"),
    comparisons = list(
      "468_0_v_468_12",
      "468_0_v_R8_0",
      "468_12_v_R8_12",
      "R8_0_v_R8_12"
    )
  )
)