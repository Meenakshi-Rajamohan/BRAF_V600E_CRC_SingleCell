# Load necessary libraries
library(Seurat)
library(Matrix)
library(ggplot2)
library(cowplot)

# Paths
input_dir <- "data/"
output_dir <- "output/"

# Load data
alldata <- readRDS(file.path(input_dir, "1_Data.rds"))

# QC Metrics
source("scripts/utils.R")  # Load helper functions
alldata <- calculate_qc_metrics(alldata)

# Save QC plots
save_violin_plot(alldata, output_dir, "QC_pre_filter.png")
alldata <- filter_data(alldata)

# Save processed data
saveRDS(alldata, file.path(output_dir, "data_after_qc.rds"))
