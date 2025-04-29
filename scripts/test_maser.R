#!/usr/bin/env Rscript

# Test script to verify MASER is working properly with rMATS results
library(maser)
library(ggplot2)




# Use absolute path to rMATS results
project_dir <- normalizePath("/home/tealeave/projects/splicing/splicing-paper")
rmats_dir <- file.path(project_dir, "RMATS_results")
comparisons <- c("468_0_v_468_12", "468_0_v_R8_0", "468_12_v_R8_12", "R8_0_v_R8_12")
event_types <- c("SE", "RI", "A3SS", "A5SS", "MXE")

# project_dir/output/test/logs
# Create output directories if they don't exist
output_dir <- file.path(project_dir, "output", "test")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
logs_dir <- file.path(output_dir, "logs")
if (!dir.exists(logs_dir)) {
  dir.create(logs_dir, recursive = TRUE)
}
# Create a timestamped log file name
log_file_name <- paste0("maser_test_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
log_file_path <- file.path(logs_dir, log_file_name)

# Redirect console output (stdout and stderr/messages) to the log file
# open a connection for appending
zz <- file(log_file_path, open = "a")

# redirect standard output (cat, print, etc):
sink(zz, append = TRUE, split = TRUE)

# redirect messages (message(), warnings, errors) as well:
sink(zz, append = TRUE, type = "message")

# now both cat() and message() go into maser_test.log
message("🚀 Starting MASER test at ", Sys.time())

cat("MASER version:", as.character(packageVersion("maser")), "\n")

# Function to test loading rMATS results
test_maser_loading <- function(comparison) {
  cat("\n===== Testing comparison:", comparison, "=====\n")
  
  # Get the full path to the rMATS results
  path <- file.path(rmats_dir, comparison)
  cat("Path:", path, "\n")
  
  # Check if the directory exists
  if (!dir.exists(path)) {
    cat("ERROR: Directory does not exist\n")
    return(NULL)
  }
  
  # Check if the required files exist
  for (event_type in event_types) {
    file_path <- file.path(path, paste0(event_type, ".MATS.JCEC.txt"))
    if (file.exists(file_path)) {
      cat("File exists:", file_path, "\n")
      # Check file size
      file_size <- file.info(file_path)$size
      cat("  Size:", file_size, "bytes\n")
      
      # Check if file has content
      if (file_size > 0) {
        # Read the first few lines to check structure
        lines <- readLines(file_path, n = 5)
        cat("  Header:", lines[1], "\n")
      } else {
        cat("  WARNING: File is empty\n")
      }
    } else {
      cat("WARNING: File does not exist:", file_path, "\n")
    }
  }
  
  # Try to create MASER object with different constructor approaches
  cat("\nAttempting to create MASER object...\n")
  
  # Method 1: Basic constructor
  maser_obj <- tryCatch({
    obj <- maser(path)
    cat("SUCCESS: MASER object created with basic constructor\n")
    obj
  }, error = function(e) {
    cat("ERROR with basic constructor:", e$message, "\n")
    NULL
  })
  
  # If we have a MASER object, try to access some data
  if (!is.null(maser_obj)) {
    cat("\nTesting data access...\n")
    
    # Check the class and structure of the MASER object
    cat("MASER object class:", class(maser_obj), "\n")
    cat("MASER object slot names:", paste(slotNames(maser_obj), collapse=", "), "\n")
    
    # Try to get summary
    tryCatch({
      summary_data <- summary(maser_obj)
      cat("Summary available\n")
      print(head(summary_data))
    }, error = function(e) {
      cat("ERROR getting summary:", e$message, "\n")
    })
    
    # Try to access each event type directly
    for (event_type in event_types) {
      tryCatch({
        # Try to access the event type directly as a slot
        if (event_type %in% slotNames(maser_obj)) {
          event_data <- slot(maser_obj, event_type)
          if (!is.null(event_data)) {
            # Try to access metadata if it exists
            if ("metadata" %in% slotNames(event_data)) {
              cat("Event type", event_type, "has", nrow(event_data@metadata), "events\n")
            } else {
              cat("Event type", event_type, "exists but has no metadata slot\n")
            }
          } else {
            cat("Event type", event_type, "slot exists but is NULL\n")
          }
        } else {
          cat("Event type", event_type, "slot does not exist\n")
        }
      }, error = function(e) {
        cat("ERROR accessing event type", event_type, ":", e$message, "\n")
      })
    }
    
    # Try to filter events for one event type
    if (length(event_types) > 0) {
      test_event_type <- event_types[1]
      cat("\nTesting filtering for event type:", test_event_type, "\n")
      
      # Test if filterByCoverage exists and try to use it
      if (exists("filterByCoverage", where = asNamespace("maser"), mode = "function")) {
        tryCatch({
          filtered <- filterByCoverage(maser_obj, test_event_type)
          cat("filterByCoverage successful\n")
        }, error = function(e) {
          cat("ERROR with filterByCoverage:", e$message, "\n")
        })
      } else {
        cat("filterByCoverage function does not exist in this version of MASER\n")
      }
      
      # Test if filterByPval exists and try to use it
      if (exists("filterByPval", where = asNamespace("maser"), mode = "function")) {
        tryCatch({
          filtered <- filterByPval(maser_obj, test_event_type, fdr = 0.05)
          cat("filterByPval successful\n")
        }, error = function(e) {
          cat("ERROR with filterByPval:", e$message, "\n")
        })
      } else {
        cat("filterByPval function does not exist in this version of MASER\n")
      }
      
      # Test if filterByDeltaPSI exists and try to use it
      if (exists("filterByDeltaPSI", where = asNamespace("maser"), mode = "function")) {
        tryCatch({
          filtered <- filterByDeltaPSI(maser_obj, test_event_type, dpsi = 0.1)
          cat("filterByDeltaPSI successful\n")
        }, error = function(e) {
          cat("ERROR with filterByDeltaPSI:", e$message, "\n")
        })
      } else {
        cat("filterByDeltaPSI function does not exist in this version of MASER\n")
      }
      
      # Try to use topEvents if it exists
      if (exists("topEvents", where = asNamespace("maser"), mode = "function")) {
        tryCatch({
          top_events <- topEvents(maser_obj, fdr = 0.05, deltaPSI = 0.1)
          cat("topEvents successful\n")
        }, error = function(e) {
          cat("ERROR with topEvents:", e$message, "\n")
        })
      } else {
        cat("topEvents function does not exist in this version of MASER\n")
      }
      
      # Try to use volcano if it exists
      if (exists("volcano", where = asNamespace("maser"), mode = "function")) {
        tryCatch({
          # Just test if the function exists, don't actually create the plot
          cat("volcano function exists\n")
        }, error = function(e) {
          cat("ERROR with volcano:", e$message, "\n")
        })
      } else {
        cat("volcano function does not exist in this version of MASER\n")
      }
    }
  }
  
  return(maser_obj)
}

# Test each comparison
for (comparison in comparisons) {
  test_maser_loading(comparison)
}

cat("\nMASER testing completed\n")
