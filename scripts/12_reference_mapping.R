# Script: 11_reference_mapping.R
# Purpose: Compare our clusters to the original Nowakowski paper labels

library(Seurat)
library(ggplot2)

# 1. Load data
# cortex_clean <- readRDS("processed/nowakowski_seurat_clustered.rds")

# 2. Check metadata for original labels
# We look for columns that look like "WGCNA", "Cluster", or "CellType"
print("Available Metadata Columns:")
print(colnames(cortex_clean@meta.data))

# 3. Visualize Side-by-Side
# We assume the column is named 'WGCNAcluster' or similar based on the paper.
# If you see a different relevant name in the printout above, change it here!
# (Note: Standard Nowakowski metadata usually has 'WGCNAcluster')

p1 <- DimPlot(cortex_clean, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("Our Clusters (Res 0.5)")

# TRYING to plot the original labels. 
# If 'WGCNAcluster' doesn't exist, this line might fail. 
# Check the console output from step 2 first!
p2 <- DimPlot(cortex_clean, reduction = "umap", group.by = "WGCNAcluster", label = TRUE, repel = TRUE) + 
  ggtitle("Original Paper Labels")

# Combine and Plot
print(p1 + p2)
ggsave("processed/comparison_umap.png", width = 14, height = 6)