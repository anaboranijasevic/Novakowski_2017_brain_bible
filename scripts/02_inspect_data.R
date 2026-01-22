# Script: 02_inspect_data.R
# Purpose: Load, clean, and inspect the Nowakowski dataset

library(readr)
library(dplyr)

# --- CONFIGURATION ---
meta_path <- "data/meta.tsv"
expr_path <- "data/exprMatrix.tsv.gz"

# ==============================================================================
# PART 1: LOAD METADATA (The "Clinical Chart")
# ==============================================================================
print("--- 1. Loading Metadata ---")

# We use read_tsv with quote="" to handle messy text errors
metadata <- read_tsv(meta_path, quote = "", show_col_types = FALSE) %>%
  as.data.frame()

# Set the cell IDs (first column) as row names
rownames(metadata) <- metadata[[1]]
metadata <- metadata[, -1] # Remove the now-redundant ID column

print(paste("Metadata loaded:", nrow(metadata), "cells processed."))
print("Column names available:")
print(colnames(metadata))


# ==============================================================================
# PART 2: LOAD EXPRESSION DATA (The "Genes")
# ==============================================================================
print("--- 2. Loading Expression Matrix (This takes time...) ---")

# We use read.table for the matrix because it handles row.names easier for base R
expr_data <- read.table(expr_path, header = TRUE, row.names = 1, sep = "\t")

print(paste("Expression data loaded:", nrow(expr_data), "genes x", ncol(expr_data), "cells."))


# ==============================================================================
# PART 3: SANITY CHECKS (The "Supervisor's Check")
# ==============================================================================
print("--- 3. Running Integrity Checks ---")

# CHECK A: Do the cell names match?
# This checks if the patients in the chart match the patients in the sequencing machine.
cells_in_meta <- rownames(metadata)
cells_in_data <- colnames(expr_data)

# Calculate intersection
common_cells <- intersect(cells_in_meta, cells_in_data)
print(paste("Number of overlapping cells:", length(common_cells)))

if (length(common_cells) < nrow(metadata)) {
  print("WARNING: Some metadata cells are missing expression data (or vice versa).")
} else {
  print("SUCCESS: All metadata cells found in expression data.")
}


# CHECK B: Data Type (Raw vs Normalized)
print("--- Data Peek (Top Left Corner) ---")
print(expr_data[1:5, 1:5])

# Automatic detection logic
first_val <- expr_data[1,1]
if (first_val %% 1 == 0) {
  print("OBSERVATION: Data appears to be INTEGERS (Raw Counts).")
} else {
  print("OBSERVATION: Data appears to be DECIMALS (Normalized TPM/CPM).")
}

