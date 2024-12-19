# Load libraries and data
library(Seurat)
library(Matrix)
library(ggplot2)
library(cowplot)

input_dir <- "output/"
output_dir <- "output/"

# Load QC data
alldata <- readRDS(file.path(input_dir, "data_after_qc.rds"))

# PCA and UMAP
alldata <- ScaleData(alldata, vars.to.regress = c("percent_mito", "nFeature_RNA"))
alldata <- RunPCA(alldata, npcs = 30)
alldata <- RunUMAP(alldata, dims = 1:30)

# Save plots
save_dim_reduction_plots(alldata, output_dir)
