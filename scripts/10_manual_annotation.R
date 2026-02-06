# Script: 10_manual_annotation.R
# Purpose: Find marker genes for manual cell type identification

library(Seurat)
library(dplyr)

# 1. Load data if needed
# cortex_clean <- readRDS("processed/nowakowski_seurat_clustered.rds")

# 2. Find Markers for ALL clusters
# This compares each cluster against the rest of the dataset.
print("Finding markers (this takes 2-3 minutes)...")
markers <- FindAllMarkers(cortex_clean, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)

# 3. Save the full table for your records
write.csv(markers, "processed/cluster_markers_manual.csv")

# 4. Display the Top 5 Genes per Cluster
# This is what we need to look at right now to name them.
top_genes <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

print("Top 5 Markers per Cluster:")
print(top_genes, n = 60) # Force print all rows