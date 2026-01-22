# Script: 03_create_object.R
# Purpose: Clean gene names and build the final Seurat Object

library(Seurat)
library(dplyr)

# Check if inputs exist in the environment (Safety Check)
if (!exists("expr_data") || !exists("metadata")) {
  stop("Error: 'expr_data' or 'metadata' not found. Please run the data loading script first.")
}

# --- 1. Clean the Gene Names ---
# Currently they look like "SOX2|SOX2". We want just "SOX2".
print("Cleaning gene names...")

# Split the name at the "|" and keep only the first part
# We use 'unname' to remove the old names from the vector for a cleaner object
clean_genes <- unname(sapply(strsplit(rownames(expr_data), "\\|"), `[`, 1))

# Handle duplicates (e.g. if two variants of SOX2 exist)
# make.unique adds .1, .2 to the end if duplicates are found
# Note: Ensure expr_data is the matrix, not the whole object yet
rownames(expr_data) <- make.unique(clean_genes)

print("Gene names cleaned. Example:")
print(head(rownames(expr_data)))

# --- 2. Create the Seurat Object ---
print("Creating Seurat Object...")

# Create the object using the cleaned matrix and metadata
cortex <- CreateSeuratObject(
  counts = expr_data,
  meta.data = metadata,
  project = "Nowakowski_Brain"
)

# --- 3. Save the Checkpoint ---
print("Saving checkpoint...")

# Safety: Create the directory if it doesn't exist yet
if (!dir.exists("processed")) {
  dir.create("processed")
}

# Save the processed object
saveRDS(cortex, file = "processed/nowakowski_seurat.rds")

print("SUCCESS! Object created and saved to 'processed/nowakowski_seurat.rds'")
print(cortex)
