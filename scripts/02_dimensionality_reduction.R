library(Seurat)
library(Matrix)
library(knitr)
library(cowplot)
library(ggplot2)
library(scran)

setwd("/Users/merajam/Documents/CRC/2_Analysis/2_DR/")


suppressWarnings(suppressMessages(alldata <- FindVariableFeatures(alldata, selection.method =
                                                                    "vst", nfeatures = 2000, verbose = FALSE, assay = "RNA")))

top20 <- head(VariableFeatures(alldata), 20)
png("DR_1_high_variable_genes.png", units="in", width=10, height=10, res=300)
LabelPoints(plot = VariableFeaturePlot(alldata), points = top20, repel = TRUE)
dev.off()

## -------------------------------------------------------------------------------------------
print("Z-score transformation !")
alldata <- ScaleData(alldata, vars.to.regress = c("percent_mito", "nFeature_RNA"), assay = "RNA")

## -------------------------------------------------------------------------------------------
print("PCA !")
alldata <- RunPCA(alldata, npcs = 50, verbose = F)

## ---- fig.asp=.28---------------------------------------------------------------------------
png("DR_2_PCs.png", units="in", width=20, height=10, res=300)
plot_grid(ncol = 3,
          DimPlot(alldata, reduction = "pca", group.by = "orig.ident", dims = 1:2),
          DimPlot(alldata, reduction = "pca", group.by = "orig.ident", dims = 3:4),
          DimPlot(alldata, reduction = "pca", group.by = "orig.ident", dims = 5:6))
dev.off()

## ----fig.asp=.5-----------------------------------------------------------------------------
png("DR_3_PCs.png", units="in", width=20, height=10, res=300)
VizDimLoadings(alldata, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
dev.off()

## ----fig.asp=.3-----------------------------------------------------------------------------
png("DR_4_variance_PCs.png", units="in", width=20, height=10, res=300)
ElbowPlot(alldata, reduction = "pca", ndims = 50)
dev.off()

## ----fig.asp=1------------------------------------------------------------------------------
print("tSNE !")
alldata <- RunTSNE(alldata, reduction = "pca", dims = 1:30,
                   perplexity = 30,
                   max_iter = 1000,
                   theta = 0.5,
                   eta = 200,
                   num_threads = 0)

## ----fig.asp=.28----------------------------------------------------------------------------
png("DR_5_tsne.png", units="in", width=10, height=10, res=300)
plot_grid(ncol = 1, DimPlot(alldata, reduction = "tsne", group.by = "orig.ident"))
dev.off()

## -------------------------------------------------------------------------------------------
print("RunUMAP !")
alldata <- RunUMAP(alldata, reduction = "pca", dims = 1:30,
                   n.components = 2,
                   n.neighbors = 30,
                   n.epochs = 200,
                   min.dist = 0.3,
                   learning.rate = 1,
                   spread = 1)

# Save the processed data
saveRDS(alldata, "CRC_merged_afterDR.rds")
save.image("CRC_merged_DR.RData")
