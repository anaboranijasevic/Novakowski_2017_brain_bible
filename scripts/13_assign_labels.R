# Script: 13_assign_labels.R
# Purpose: Permanently rename clusters to biological cell types

library(Seurat)
library(ggplot2)

# 1. Load the clustered object
# cortex_clean <- readRDS("processed/nowakowski_seurat_clustered.rds")

# 2. Define the Identity Mapping
# The order MUST match the cluster numbers (0, 1, 2...) exactly.
# Based on our Manual Annotation + Reference Mapping confirmation:
new_ids <- c(
  "Excitatory Neurons (Deep)",  # Cluster 0
  "MGE Interneurons",           # Cluster 1
  "Excitatory Neurons (Upper)", # Cluster 2
  "Dividing Progenitors",       # Cluster 3
  "Radial Glia",                # Cluster 4
  "CGE Interneurons",           # Cluster 5
  "Excitatory Neurons (Deep)",  # Cluster 6 (Merging with 0)
  "Meninges/Vascular",          # Cluster 7
  "Microglia",                  # Cluster 8
  "Unknown/Outliers",           # Cluster 9
  "Pericytes",                  # Cluster 10
  "OPC"                         # Cluster 11
)

# 3. Apply the Names
names(new_ids) <- levels(cortex_clean)
cortex_clean <- RenameIdents(cortex_clean, new_ids)

# 4. Stash the names in metadata
# This allows you to switch back and forth between "seurat_clusters" (numbers) and "cell_type" (names)
cortex_clean$cell_type <- Idents(cortex_clean)

# 5. Generate the Final Annotated UMAP
print("Generating Final Annotated Map...")

# We use 'label = TRUE' to write the names on the map
p1 <- DimPlot(cortex_clean, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Final Annotated Cell Types") +
  NoLegend() # The labels are enough

# 6. Save Plot and Object
print(p1)
ggsave("processed/final_annotated_umap.png", plot = p1, width = 10, height = 8)

saveRDS(cortex_clean, "processed/nowakowski_final_annotated.rds")
print("SUCCESS! Analysis Complete. Final object saved as 'processed/nowakowski_final_annotated.rds'")