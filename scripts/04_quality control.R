# Script: 04_quality_control.R
# Purpose: Calculate QC metrics and visualize them

library(Seurat)
library(ggplot2)

# 1. Load the object you created in the last step
# (If it's not currently in your environment)
# cortex <- readRDS("processed/nowakowski_seurat.rds")

print("Calculating mitochondrial percentage...")

# Calculate percentage of counts from mitochondrial genes (start with "MT-")
cortex[["percent.mt"]] <- PercentageFeatureSet(cortex, pattern = "^MT-")

# 2. Visualize QC Metrics
print("Generating QC Violin Plot...")

# We plot 3 things:
# nFeature_RNA: Genes per cell
# nCount_RNA:   Molecules per cell
# percent.mt:   Mitochondrial percentage (measure of cell health)
VlnPlot(cortex, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3,
        pt.size = 0.1) 

# 3. Save the plot
ggsave("processed/qc_violin_plot.png", width = 10, height = 6)
print("Plot saved to 'processed/qc_violin_plot.png'")