# Script: 06_preprocessing.R
# Purpose: Normalize, Find Variable Genes, Scale, and Run PCA

library(Seurat)
library(dplyr)

# 1. Load the filtered object
# (If you just finished the last step, it's already called 'cortex_clean')
# If you restarted R, uncomment the line below:
# cortex_clean <- readRDS("processed/nowakowski_seurat_filtered.rds")

# 2. Log-Normalize
# Since your data is already CPM (counts per million), this effectively
# just performs the log-transformation: log(count + 1)
print("Normalizing data...")
cortex_clean <- NormalizeData(cortex_clean, 
                              normalization.method = "LogNormalize", 
                              scale.factor = 10000)

# 3. Find Variable Features
# We identify the top 2000 genes that vary the most across cells
print("Finding variable features...")
cortex_clean <- FindVariableFeatures(cortex_clean, 
                                     selection.method = "vst", 
                                     nfeatures = 2000)

# Quick check: What are the most variable genes? (Usually markers like SOX2, VIM)
top10 <- head(VariableFeatures(cortex_clean), 10)
print(paste("Top variable genes:", paste(top10, collapse = ", ")))


# 4. Scale Data
# This centers the data (mean = 0) so PCA works correctly.
# By default, it only scales the top 2000 variable genes to save memory.
print("Scaling data...")
cortex_clean <- ScaleData(cortex_clean)

# 5. Run PCA
print("Running PCA...")
cortex_clean <- RunPCA(cortex_clean, 
                       features = VariableFeatures(object = cortex_clean))

# 6. Save the checkpoints
# We save this object because PCA takes time, and we don't want to redo it.
saveRDS(cortex_clean, "processed/nowakowski_seurat_pca.rds")

print("SUCCESS! Preprocessing complete.")