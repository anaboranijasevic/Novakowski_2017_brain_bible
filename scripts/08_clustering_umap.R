# Script: 08_clustering_analysis.R
# Purpose: Determine optimal clustering resolution, visualize stability, and generate UMAP.

library(Seurat)
library(dplyr)
library(ggplot2)
library(glue)
library(clustree)

# 1. Load the PCA object
# (Uncomment the line below if starting fresh in a new session)
# cortex_clean <- readRDS("processed/nowakowski_seurat_pca.rds")

# --- Step 1: Find Neighbors ---
# We use the top 20 PCs (determined from the Elbow Plot previously)
print("Finding Neighbors...")
cortex_clean <- FindNeighbors(cortex_clean, dims = 1:20)

# --- Step 2: Explore Resolutions (0.1 to 1.0) ---
# We calculate clusters at many resolutions to see how stable they are.
print("Calculating clusters for resolutions 0.1 to 1.0...")
cortex_clean <- FindClusters(cortex_clean, resolution = seq(from = 0.1, to = 1.0, by = 0.1))

# --- Step 3: Generate & Display Clustree ---
print("Generating Clustree Plot...")
clustree_plot <- clustree::clustree(cortex_clean, 
                                    prefix = "RNA_snn_res.", 
                                    highlight_core = TRUE, 
                                    node_size = 5, 
                                    edge_arrow = TRUE, 
                                    edge_arrow_ends = "first", 
                                    edge_width = 0.5) + 
  ggtitle(glue('Clustree Analysis: {ncol(cortex_clean)} cells'))

# SAVE the plot to file
ggsave("processed/clustree_plot.png", plot = clustree_plot, width = 12, height = 10)

# DISPLAY the plot to screen (Force it to appear)
print(clustree_plot)
print("Clustree displayed and saved to 'processed/clustree_plot.png'")

# --- Step 4: Set Final Identity (Resolution 0.5) ---
# Based on our analysis, 0.5 is the most stable resolution.
print("Setting active identity to resolution 0.5...")
Idents(cortex_clean) <- "RNA_snn_res.0.5"

# --- Step 5: Run UMAP & Display ---
print("Running UMAP...")
cortex_clean <- RunUMAP(cortex_clean, dims = 1:20)

print("Generating UMAP Plot...")
umap_plot <- DimPlot(cortex_clean, reduction = "umap", label = TRUE, label.size = 5) +
  ggtitle("Final UMAP Clusters (Res 0.5)")

# SAVE the plot to file
ggsave("processed/umap_clusters.png", plot = umap_plot, width = 8, height = 6)

# DISPLAY the plot to screen (Force it to appear)
print(umap_plot)
print("UMAP displayed and saved to 'processed/umap_clusters.png'")

# --- Step 6: Save the Gold Master Object ---
# We stash the object identity so it defaults to res.0.5 when you load it next time.
cortex_clean[["seurat_clusters"]] <- Idents(cortex_clean)

saveRDS(cortex_clean, "processed/nowakowski_seurat_clustered.rds")
print("SUCCESS! Final object saved. Ready for marker identification.")
                     