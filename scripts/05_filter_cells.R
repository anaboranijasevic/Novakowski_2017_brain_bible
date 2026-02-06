# Script: 05_filter_cells.R
# Purpose: Remove low quality cells and doublets

library(Seurat)
library(dplyr)

# 1. Load the object if needed
# cortex <- readRDS("processed/nowakowski_seurat.rds")

print("Dimensions before filtering:")
print(dim(cortex))

# 2. Define Cutoffs
# nFeature_RNA > 1000: Remove empty droplets (low quality)
# nFeature_RNA < 10000: Remove doublets (cells with insanely high gene counts)
# percent.mt < 5: Standard cutoff (Keep this even if currently 0, it won't hurt)

cortex_clean <- subset(cortex, subset = nFeature_RNA > 1000 & 
                         nFeature_RNA < 10000 & 
                         percent.mt < 5)

print("Dimensions after filtering:")
print(dim(cortex_clean))

# 3. Save the clean object
saveRDS(cortex_clean, "processed/nowakowski_seurat_filtered.rds")
print("Filtered object saved successfully.")