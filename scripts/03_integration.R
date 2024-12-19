library(Seurat)
library(Matrix)
library(knitr)
library(cowplot)
library(ggplot2)
library(scran)
library(harmony)
library(pheatmap)

setwd("/Users/merajam/Documents/CRC/2_Analysis/3_Integration/")

khaliq <- readRDS("/Users/merajam/Documents/CRC/2_Analysis/2_DR/khaliq/CRC_merged_afterDR.rds")
tan <- readRDS("/Users/merajam/Documents/CRC/2_Analysis/2_DR/tan/CRC_merged_afterDR.rds")

khaliq$dataset <- "khaliq"
tan$dataset <- "tan"

DefaultAssay(khaliq) <- "RNA"
DefaultAssay(tan) <- "RNA"

options(future.globals.maxSize = 200000000 * 1024^2)

alldata.list <- lapply(X = c(khaliq, tan), FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- SCTransform(object = x, verbose = FALSE, return.only.var.genes = FALSE)
})

alldata.features <- SelectIntegrationFeatures(object.list = alldata.list, nfeatures = 2000)
alldata.list <- PrepSCTIntegration(object.list = alldata.list, anchor.features = alldata.features)

print("Finding anchors")
alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, normalization.method = "SCT",
                                          dims = 1:25, anchor.features = alldata.features)

alldata.integrated <- IntegrateData(anchorset = alldata.anchors, normalization.method = "SCT",
                                    dims = 1:25)

alldata.integrated <- ScaleData(alldata.integrated, verbose = FALSE)
alldata.integrated <- RunPCA(alldata.integrated, pcs = 30, verbose = FALSE, reduction.name = "PCA_on_CCA")
alldata.integrated <- RunUMAP(alldata.integrated, reduction = "PCA_on_CCA", dims = 1:25,
                              reduction.name = "UMAP_on_CCA")

# Save integrated data
saveRDS(alldata.integrated, "alldata.integrated.rds")
save.image("CRC_integrated.RData")
