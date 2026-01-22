# Script: 01_setup.R
# Purpose: Load libraries for single-cell analysis

library(Seurat)
library(tidyverse)

print("Setup loaded successfully!")
# 1. Define the path to your data
# (Replace 'data/exprMatrix.tsv.gz' with the actual name of the file you downloaded)
data_path <- "data/meta.tsv" 

# 2. Check if the computer can find it
if (file.exists(data_path)) {
  print("File found! Ready to load.")
} else {
  print("WARNING: File not found. Check the spelling.")
}

