# Script: 11_assign_labels.R
# Purpose: Rename the numbered clusters to real biological names

library(Seurat)
library(ggplot2)

# 1. Load data
# cortex_clean <- readRDS("processed/nowakowski_seurat_clustered.rds")

# 2. Define the list of names (Based on our analysis)
# The order MUST match the cluster numbers (0, 1, 2...)
new_cluster_ids <- c(
  "Excitatory Neurons (Deep)", # 0
  "MGE Interneurons",          # 1
  "Excitatory Neurons (Upper)",# 2
  "Dividing Progenitors",      # 3
  "Radial Glia",               # 4
  "CGE Interneurons",          # 5
  "Excitatory Neurons (Deep)", # 6 (Merging with 0 is an option, or keep distinct)
  "Meninges/Vascular",         # 7
  "Microglia",                 # 8
  "Unknown/Outliers",          # 9 (The markers were vague lncRNAs)
  "Pericytes",                 # 10
  "OPC"                        # 11
)

# 3. Rename the identity
names(new_cluster_ids) <- levels(cortex_clean)
cortex_clean <- RenameIdents(cortex_clean, new_cluster_ids)

# 4. Save the "Active Identity" to a new metadata column
cortex_clean$cell_type <- Idents(cortex_clean)

# 5. Visualize the Final Annotated Map
print("Generating Annotated UMAP...")
p1 <- DimPlot(cortex_clean, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Annotated Cell Types") +
  NoLegend() # The labels on the plot are enough

print(p1)
ggsave("processed/final_annotated_umap.png", width = 10, height = 8)

# 6. Save the Final Object
saveRDS(cortex_clean, "processed/nowakowski_final_annotated.rds")
print("SUCCESS! Analysis Complete.")