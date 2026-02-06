# Script: 07_pca_visualization.R
# Purpose: Visualize PCA and determine how many PCs to use

library(Seurat)
library(ggplot2)

# 1. Load data (if needed)
# cortex_clean <- readRDS("processed/nowakowski_seurat_pca.rds")

# 2. Visualize the first 2 Principal Components
# This maps the 56,000 genes down to just 2 dimensions (X and Y)
print("Generating PCA Plot...")
p1 <- DimPlot(cortex_clean, reduction = "pca") + 
  ggtitle("PCA Plot (PC1 vs PC2)")
print(p1)

# Look at the genes driving the first dimension
print("Genes driving PC 1:")
print(Loadings(cortex_clean[["pca"]], dims = 1)[1:5, ])

# 3. The Elbow Plot
# This helps us decide how many PCs to use for clustering.
# We look for the "Elbow" where the curve flattens out.
print("Generating Elbow Plot...")
p2 <- ElbowPlot(cortex_clean, ndims = 50) +
  ggtitle("Elbow Plot: How many PCs to use?") +
  geom_vline(xintercept = 20, linetype="dashed", color="red") # Reference line

print(p2)

# Save plots
ggsave("processed/pca_plot.png", plot = p1, width = 6, height = 5)
ggsave("processed/elbow_plot.png", plot = p2, width = 6, height = 5)

print("Plots saved to 'processed/' folder.")