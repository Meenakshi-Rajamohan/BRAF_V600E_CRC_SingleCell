library(Seurat)
library(Matrix)
library(DoubletFinder)
library(knitr)
library(cowplot)
library(ggplot2)
library(scran)
setwd("/Users/merajam/Documents/CRC/2_Analysis/1_QC/BRAF/")
alldata <- readRDS("/Users/merajam/Documents/CRC/1_Data.rds")
# Way1: Doing it using Seurat function
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")
# Way2: Doing it manually
total_counts_per_cell <- colSums(alldata@assays$RNA@counts)
mito_genes <- rownames(alldata)[grep("^MT-", rownames(alldata))]
alldata$percent_mito <- colSums(alldata@assays$RNA@counts[mito_genes, ])/total_counts_per_cell
head(mito_genes, 10)
# Way1: Doing it using Seurat function
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")
# Way2: Doing it manually
ribo_genes <- rownames(alldata)[grep("^RP[SL]", rownames(alldata))]
head(ribo_genes, 10)
alldata$percent_ribo <- colSums(alldata@assays$RNA@counts[ribo_genes, ])/total_counts_per_cell
# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")
alldata <- PercentageFeatureSet(alldata, "PECAM1|PF4", col.name = "percent_plat")
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo",
           "percent_hb","percent_plat")
png("QC_1_sc_preQC_CRC_BRAF.png", units="in", width=20, height=10, res=300)
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + NoLegend()
dev.off()
png("QC_2_sc_FeatureScatter_CRC_BRAF.png" , units="in", width=10, height=8, res=300)
FeatureScatter(alldata, "nCount_RNA" , "nFeature_RNA", group.by = "orig.ident", pt.size = .5)
dev.off()
print("Before QC")
dim(alldata)
table(alldata$orig.ident)
png("QC_3_most_expressed.png", units="in", width=20, height=10, res=300)
par(mar=c(4,8,2,1))
C <- alldata@assays$RNA@counts
C <- Matrix::t( Matrix::t(C) / Matrix::colSums(C) ) * 100
most_expressed <- order(apply(C,1,median),decreasing = T)[20:1]
boxplot( as.matrix(t(C[most_expressed,])),cex=.1, las=1, xlab="% total count per
cell",col=scales::hue_pal()(20)[20:1],horizontal=TRUE)
dev.off()
selected_c <- WhichCells(alldata, expression = nFeature_RNA > 300)
selected_f <- rownames(alldata)[Matrix::rowSums(alldata) > 3]
data.filt <- subset(alldata, features = selected_f, cells = selected_c)
dim(data.filt)
selected_mito <- WhichCells(data.filt, expression = percent_mito < 0.2)
#selected_ribo <- WhichCells(data.filt, expression = percent_ribo > 0.05)
# and subset the object to only keep those cells
data.filt <- subset(data.filt, cells = selected_mito)
#data.filt <- subset(data.filt, cells = selected_ribo)
dim(data.filt)
table(data.filt$orig.ident)
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito","percent_ribo")
png("QC_4_sc_PostQC_CRC_BRAF.png", units="in", width=20, height=10, res=300)
VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 4) +
  NoLegend()
dev.off()
data.filt <- data.filt[!grepl("MALAT1", rownames(data.filt)), ]
# Filter Mitocondrial
data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]
# Filter Ribossomal gene (optional if that is a problem on your data) data.filt
data.filt <- data.filt[ ! grepl('^RP[SL]', rownames(data.filt)), ]
# Filter Hemoglobin gene (optional if that is a problem on your data)
data.filt <- data.filt[!grepl("^HB[^(P)]", rownames(data.filt)), ]
dim(data.filt)
chrY.gene <- read.csv("/Users/merajam/Documents/CRC/1_Data/genes.table_Y.tsv", sep ="\t")
genes.table <- chrY.gene[chrY.gene$external_gene_name %in% rownames(data.filt),]
chrY.gene = genes.table$external_gene_name[genes.table$chromosome_name == "Y"]
data.filt$pct_chrY = colSums(data.filt@assays$RNA@counts[chrY.gene,
])/colSums(data.filt@assays$RNA@counts)
png("QC_4_sample_sex_BRAF.png", units="in", width=15, height=7, res=300)
VlnPlot(data.filt, group.by = "orig.ident" ,features = c("XIST", "pct_chrY"))
dev.off()
data.filt = NormalizeData(data.filt)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data.filt <- CellCycleScoring(object = data.filt, s.features = s.genes, g2m.features = g2m.genes)
gc()
png("QC_5_cellcycle_score_BRAF.png", units="in", width=15, height=7, res=300)
VlnPlot(data.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
        ncol = 4, pt.size = 0.1)
dev.off()
dim(data.filt)
data.filt = FindVariableFeatures(data.filt, verbose = F)
data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "percent_mito"), verbose =
                        F)
gc()
data.filt = RunPCA(data.filt, verbose = F, npcs = 20)
gc()
data.filt = RunUMAP(data.filt, dims = 1:10, verbose = F)
gc()
sweep.res <- paramSweep_v3(data.filt)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn_data.filt <- find.pK(sweep.stats)
gc()
png("QC_6_paramSweep.bcmvn_data.filt_CRC_BRAF.png", units="in", width=10, height=7, res=300)
barplot(bcmvn_data.filt$BCmetric, names.arg = bcmvn_data.filt$pK, las=2)
dev.off()
pK = bcmvn_data.filt$pK[which.max(bcmvn_data.filt$BCmetric)]
pK = as.numeric(levels(pK))[pK]
# define the expected number of doublet cellscells.
nExp <- round(ncol(data.filt) * 0.04) # expect 4% doublets
gc()
data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = pK, nExp = nExp, PCs = 1:10)
gc()
# name of the DF prediction can htange, so extract the correct column name.
DF.name = colnames(data.filt@meta.data)[grepl("DF.classification",
                                              colnames(data.filt@meta.data))]
png("QC_7_doublet_data.filt_CRC_BRAF.png", units="in", width=15, height=7, res=300)
cowplot::plot_grid(ncol = 2,
                   DimPlot(data.filt, group.by = "orig.ident") + NoAxes(),
                   DimPlot(data.filt, group.by = DF.name) + NoAxes())
dev.off()
data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
gc()
png("QC_7a_removed_doublet_data.filt_CRC_BRAF.png", units="in", width=15, height=7, res=300)
cowplot::plot_grid(ncol = 2,
                   DimPlot(data.filt, group.by = "orig.ident") + NoAxes(),
                   DimPlot(data.filt, group.by = DF.name) + NoAxes())
dev.off()
dim(data.filt)
png("QC_9_most_expressed_after_QC_CRC_BRAF.png", units="in", width=20, height=10, res=300)
par(mar=c(4,8,2,1))
C <- data.filt@assays$RNA@counts
C <- Matrix::t( Matrix::t(C) / Matrix::colSums(C) ) * 100
most_expressed <- order(apply(C,1,median),decreasing = T)[20:1]
boxplot( as.matrix(t(C[most_expressed,])),cex=.1, las=1, xlab="% total count per
cell",col=scales::hue_pal()(20)[20:1],horizontal=TRUE)
dev.off()
gc()
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito","percent_ribo")
png("QC_10_sc_preQC_CRC_BRAF_sc_PostQC_CRC_BRAF.png", units="in", width=20, height=10,
    res=300)
VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 4) +
  NoLegend()
dev.off()
saveRDS(data.filt, "CRC_BRAF_afterQC_BRAF.rds")
save.image("CRC_BRAF_QC.RData")
#savehistory("CRC_BRAF_QC.R")
setwd("/Users/merajam/Documents/CRC/2_Analysis/1_QC/NonBRAF/")
alldata <- readRDS("/Users/merajam/Documents/CRC/1_Data/NonBRAF.rds")
# Way1: Doing it using Seurat function
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")
# Way2: Doing it manually
total_counts_per_cell <- colSums(alldata@assays$RNA@counts)
mito_genes <- rownames(alldata)[grep("^MT-", rownames(alldata))]
alldata$percent_mito <- colSums(alldata@assays$RNA@counts[mito_genes, ])/total_counts_per_cell
head(mito_genes, 10)
# Way1: Doing it using Seurat function
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")
# Way2: Doing it manually
ribo_genes <- rownames(alldata)[grep("^RP[SL]", rownames(alldata))]
head(ribo_genes, 10)
alldata$percent_ribo <- colSums(alldata@assays$RNA@counts[ribo_genes, ])/total_counts_per_cell
# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")
alldata <- PercentageFeatureSet(alldata, "PECAM1|PF4", col.name = "percent_plat")
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo",
           "percent_hb","percent_plat")
png("QC_1_sc_preQC_CRC_NonBRAF.png", units="in", width=20, height=10, res=300)
VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + NoLegend()
dev.off()
png("QC_2_sc_FeatureScatter_CRC_NonBRAF.png" , units="in", width=10, height=8, res=300)
FeatureScatter(alldata, "nCount_RNA" , "nFeature_RNA", group.by = "orig.ident", pt.size = .5)
dev.off()
print("Before QC")
dim(alldata)
table(alldata$orig.ident)
png("QC_3_most_expressed.png", units="in", width=20, height=10, res=300)
par(mar=c(4,8,2,1))
C <- alldata@assays$RNA@counts
C <- Matrix::t( Matrix::t(C) / Matrix::colSums(C) ) * 100
most_expressed <- order(apply(C,1,median),decreasing = T)[20:1]
boxplot( as.matrix(t(C[most_expressed,])),cex=.1, las=1, xlab="% total count per
cell",col=scales::hue_pal()(20)[20:1],horizontal=TRUE)
dev.off()
selected_c <- WhichCells(alldata, expression = nFeature_RNA > 300)
selected_f <- rownames(alldata)[Matrix::rowSums(alldata) > 3]
data.filt <- subset(alldata, features = selected_f, cells = selected_c)
dim(data.filt)
selected_mito <- WhichCells(data.filt, expression = percent_mito < 0.2)
#selected_ribo <- WhichCells(data.filt, expression = percent_ribo > 0.05)
# and subset the object to only keep those cells
data.filt <- subset(data.filt, cells = selected_mito)
#data.filt <- subset(data.filt, cells = selected_ribo)
dim(data.filt)
table(data.filt$orig.ident)
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito","percent_ribo")
png("QC_4_sc_PostQC_CRC_NonBRAF.png", units="in", width=20, height=10, res=300)
VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 4) +
  NoLegend()
dev.off()
data.filt <- data.filt[!grepl("MALAT1", rownames(data.filt)), ]
# Filter Mitocondrial
data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]
# Filter Ribossomal gene (optional if that is a problem on your data) data.filt
data.filt <- data.filt[ ! grepl('^RP[SL]', rownames(data.filt)), ]
# Filter Hemoglobin gene (optional if that is a problem on your data)
data.filt <- data.filt[!grepl("^HB[^(P)]", rownames(data.filt)), ]
dim(data.filt)
chrY.gene <- read.csv("/Users/merajam/Documents/CRC/1_Data/genes.table_Y.tsv", sep ="\t")
genes.table <- chrY.gene[chrY.gene$external_gene_name %in% rownames(data.filt),]
chrY.gene = genes.table$external_gene_name[genes.table$chromosome_name == "Y"]
data.filt$pct_chrY = colSums(data.filt@assays$RNA@counts[chrY.gene,
])/colSums(data.filt@assays$RNA@counts)
png("QC_4_sample_sex_NonBRAF.png", units="in", width=15, height=7, res=300)
VlnPlot(data.filt, group.by = "orig.ident" ,features = c("XIST", "pct_chrY"))
dev.off()
data.filt = NormalizeData(data.filt)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
data.filt <- CellCycleScoring(object = data.filt, s.features = s.genes, g2m.features = g2m.genes)
gc()
png("QC_5_cellcycle_score_NonBRAF.png", units="in", width=15, height=7, res=300)
VlnPlot(data.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
        ncol = 4, pt.size = 0.1)
dev.off()
dim(data.filt)
data.filt = FindVariableFeatures(data.filt, verbose = F)
data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "percent_mito"), verbose =
                        F)
gc()
data.filt = RunPCA(data.filt, verbose = F, npcs = 20)
gc()
data.filt = RunUMAP(data.filt, dims = 1:10, verbose = F)
gc()
sweep.res <- paramSweep_v3(data.filt)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn_data.filt <- find.pK(sweep.stats)
gc()
png("QC_6_paramSweep.bcmvn_data.filt_CRC_NonBRAF.png", units="in", width=10, height=7, res=300)
barplot(bcmvn_data.filt$BCmetric, names.arg = bcmvn_data.filt$pK, las=2)
dev.off()
pK = bcmvn_data.filt$pK[which.max(bcmvn_data.filt$BCmetric)]
pK = as.numeric(levels(pK))[pK]
# define the expected number of doublet cellscells.
nExp <- round(ncol(data.filt) * 0.04) # expect 4% doublets
gc()
data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = pK, nExp = nExp, PCs = 1:10)
gc()
# name of the DF prediction can htange, so extract the correct column name.
DF.name = colnames(data.filt@meta.data)[grepl("DF.classification",
                                              colnames(data.filt@meta.data))]
png("QC_7_doublet_data.filt_CRC_NonBRAF.png", units="in", width=15, height=7, res=300)
cowplot::plot_grid(ncol = 2,
                   DimPlot(data.filt, group.by = "orig.ident") + NoAxes(),
                   DimPlot(data.filt, group.by = DF.name) + NoAxes())
dev.off()
data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
gc()
png("QC_7a_removed_doublet_data.filt_CRC_NonBRAF.png", units="in", width=15, height=7, res=300)
cowplot::plot_grid(ncol = 2,
                   DimPlot(data.filt, group.by = "orig.ident") + NoAxes(),
                   DimPlot(data.filt, group.by = DF.name) + NoAxes())
dev.off()
dim(data.filt)
png("QC_9_most_expressed_after_QC_CRC_NonBRAF.png", units="in", width=20, height=10, res=300)
par(mar=c(4,8,2,1))
C <- data.filt@assays$RNA@counts
C <- Matrix::t( Matrix::t(C) / Matrix::colSums(C) ) * 100
most_expressed <- order(apply(C,1,median),decreasing = T)[20:1]
boxplot( as.matrix(t(C[most_expressed,])),cex=.1, las=1, xlab="% total count per
cell",col=scales::hue_pal()(20)[20:1],horizontal=TRUE)
dev.off()
gc()
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito","percent_ribo")
png("QC_10_sc_preQC_CRC_NonBRAF_sc_PostQC_CRC_NonBRAF.png", units="in", width=20, height=10,
    res=300)
VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 4) +
  NoLegend()
dev.off()
saveRDS(data.filt, "CRC_NonBRAF_afterQC_NonBRAF.rds")
save.image("CRC_NonBRAF_QC.RData")
#savehistory("CRC_NonBRAF_QC.R")
#### Dimensionality reduction
library(Seurat)
library(Matrix)
library(knitr)
library(cowplot)
library(ggplot2)
library(scran)
setwd("/Users/merajam/Documents/CRC/2_Analysis/2_DR/")
BRAF <-
  readRDS("/Users/merajam/Documents/CRC/2_Analysis/1_QC/BRAF/CRC_BRAF_afterQC_BRAF.rds")
BRAF$category <- "BRAF"
NonBRAF <-
  readRDS("/Users/merajam/Documents/CRC/2_Analysis/1_QC/NonBRAF/CRC_NonBRAF_afterQC_NonBRAF.rds")
NonBRAF$category <- "NonBRAF"
alldata <- merge(BRAF, y = NonBRAF)
saveRDS(alldata,"merged.rds")
## -------------------------------------------------------------------------------------------
print("2- Feature selection !")
suppressWarnings(suppressMessages(alldata <- FindVariableFeatures(alldata, selection.method =
                                                                    "vst", nfeatures = 2000 ,verbose = FALSE,assay = "RNA")))
top20 <- head(VariableFeatures(alldata), 20)
png("DR_1_high_variable_genes.png", units="in", width=10, height=10, res=300)
LabelPoints(plot = VariableFeaturePlot(alldata), points = top20, repel = TRUE)
dev.off()
## -------------------------------------------------------------------------------------------
print("3- Z-score transformation !")
alldata <- ScaleData(alldata, vars.to.regress = c("percent_mito", "nFeature_RNA"), assay = "RNA")
## -------------------------------------------------------------------------------------------
print("4- PCA !")
alldata <- RunPCA(alldata, npcs = 50, verbose = F)
## ---- fig.asp=.28---------------------------------------------------------------------------
png("DR_2_PCs.png", units="in", width=20, height=10, res=300)
plot_grid(ncol = 3,
          DimPlot(alldata, reduction = "pca", group.by = "orig.ident",dims = 1:2),
          DimPlot(alldata, reduction = "pca", group.by = "orig.ident",dims = 3:4),
          DimPlot(alldata, reduction = "pca", group.by = "orig.ident",dims = 5:6) )
dev.off()
## ----fig.asp=.5-----------------------------------------------------------------------------
png("DR_3_PCs.png", units="in", width=20, height=10, res=300)
VizDimLoadings(alldata, dims = 1:5, reduction = "pca",ncol = 5,balanced = T)
dev.off()
## ----fig.asp=.3-----------------------------------------------------------------------------
png("DR_4_variance_PCs.png", units="in", width=20, height=10, res=300)
ElbowPlot(alldata, reduction = "pca",ndims = 50)
dev.off()
## ----fig.asp=1------------------------------------------------------------------------------
print("4- tSNE !")
alldata <- RunTSNE(alldata, reduction = "pca", dims = 1:30,
                   perplexity=30,
                   max_iter=1000,
                   theta=0.5,
                   eta=200,
                   num_threads=0 )
#see ?Rtsne and ?RunTSNE for more info
## ----fig.asp=.28----------------------------------------------------------------------------
png("DR_5_tsne.png", units="in", width=10, height=10, res=300)
plot_grid(ncol = 1,DimPlot(alldata, reduction = "tsne", group.by = "orig.ident"))
dev.off()
## -------------------------------------------------------------------------------------------
print("5- RunUMAP !")
alldata <- RunUMAP(alldata, reduction = "pca", dims = 1:30,
                   n.components=2,
                   n.neighbors=30,
                   n.epochs=200,
                   min.dist=0.3,
                   learning.rate=1,
                   spread=1 )
#see ?RunUMAP for more info
## -------------------------------------------------------------------------------------------
# we can add in additional reductions, by defulat they are named "pca", "umap", "tsne" etc. But
we can specify alternative names with reduction.name
alldata <- RunUMAP(alldata, reduction.name = "UMAP10_on_PCA",
                   reduction = "pca",
                   dims = 1:30,
                   n.components=10,
                   n.neighbors=30,
                   n.epochs=200,
                   min.dist=0.3,
                   learning.rate=1,
                   spread=1 )
#see ?RunUMAP for more info
## ----fig.asp=.28----------------------------------------------------------------------------
png("DR_6_umap_colored_per_dataset.png", units="in", width=20, height=10, res=300)
plot_grid(ncol = 3,
          DimPlot(alldata, reduction = "umap", group.by = "orig.ident")+ ggplot2::ggtitle(label
                                                                                          ="UMAP_on_PCA"),
          DimPlot(alldata, reduction = "UMAP10_on_PCA", group.by = "orig.ident",dims = 1:2)+
            ggplot2::ggtitle(label ="UMAP10_on_PCA"),
          DimPlot(alldata, reduction = "UMAP10_on_PCA", group.by = "orig.ident",dims = 3:4)+
            ggplot2::ggtitle(label ="UMAP10_on_PCA")
)
dev.off()
## ----fig.asp=.28----------------------------------------------------------------------------
png("DR_7_umap_pca_tsne.png", units="in", width=20, height=10, res=300)
plot_grid(ncol = 3,
          DimPlot(alldata, reduction = "pca", group.by = "orig.ident"),
          DimPlot(alldata, reduction = "tsne", group.by = "orig.ident"),
          DimPlot(alldata, reduction = "umap", group.by = "orig.ident")
)
dev.off()
## -------------------------------------------------------------------------------------------
print("6- Run_UMAP_on_ScaleData !")
alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_ScaleData",
                   features = alldata@assays$RNA@var.features,
                   assay = "RNA",
                   n.components=2,
                   n.neighbors=30,
                   n.epochs=200,
                   min.dist=0.3,
                   learning.rate=1,
                   spread=1 )
## -------------------------------------------------------------------------------------------
#Build Graph
alldata <- FindNeighbors(alldata,
                         reduction = "pca",
                         graph.name = "SNN",
                         assay = "RNA",
                         k.param = 20,
                         features = alldata@assays$RNA@var.features)
save.image("CRC_merged_before_UMAP_on_Graph.RData")
#Run UMAP on a graph
alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_Graph",
                   graph = "SNN",
                   assay = "RNA" )
## ---- fig.asp=.28---------------------------------------------------------------------------
p1 <- DimPlot(alldata, reduction = "umap", group.by = "orig.ident")+ ggplot2::ggtitle(label
                                                                                      ="UMAP_on_PCA")
p2 <- DimPlot(alldata, reduction = "UMAP_on_ScaleData", group.by = "orig.ident")+
  ggplot2::ggtitle(label ="UMAP_on_ScaleData")
p3 <- DimPlot(alldata, reduction = "UMAP_on_Graph", group.by = "orig.ident")+
  ggplot2::ggtitle(label ="UMAP_on_Graph")
leg <- get_legend(p1)
png("DR_8_pca_umap_graph.png", units="in", width=40, height=20, res=300)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    p1 + NoLegend() + NoAxes(),
    p2 + NoLegend() + NoAxes(),
    p3 + NoLegend() + NoAxes(),
    leg,nrow=2),
  ncol=1,widths=c(1)
)
dev.off()
saveRDS(alldata, "CRC_merged_afterDI.rds")
save.image("CRC_merged_DI.RData")
3. Integration
library(Seurat)
library(Matrix)
library(knitr)
library(cowplot)
library(ggplot2)
library(scran)
library(harmony)
library(pheatmap)
setwd("/Users/merajam/Documents/CRC/2_Analysis/3_Integration/")

khaliq <- readRDS("/Users/merajam/Documents/CRC/2_Analysis/2_DR/khaliq/CRC_merged_afterDI.rds")
tan <- readRDS("/Users/merajam/Documents/CRC/2_Analysis/2_DR/tan/CRC_merged_afterDI.rds")

khaliq$dataset <- "khaliq"

tan$dataset <- "tan"

DefaultAssay(khaliq) <- "RNA"
DefaultAssay(tan) <- "RNA"
options(future.globals.maxSize = 200000000 * 1024^2)
alldata.list <- lapply(X = c(lee,pelka,khaliq,tan), FUN = function(x)
{
  x <- NormalizeData(x, verbose = FALSE)
  x <- SCTransform(object = x, verbose = FALSE, return.only.var.genes = FALSE)
})
alldata.features <- SelectIntegrationFeatures(object.list = alldata.list, nfeatures = 2000)
alldata.list <- PrepSCTIntegration(object.list = alldata.list, anchor.features =
                                     alldata.features)
rm(list=setdiff(ls(), c("alldata.features", "alldata.list")))
print("Finding anchors")
alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, normalization.method =
                                            "SCT", dims = 1:25, anchor.features = alldata.features)
gc()
saveRDS(alldata.anchors , "alldata.anchors.rds")
alldata.integrated <- IntegrateData( anchorset = alldata.anchors, normalization.method = "SCT",
                                     dims = 1:25)
gc()
rm(alldata.anchors)
alldata.integrated <- ScaleData(alldata.integrated, verbose = FALSE)
alldata.integrated <- RunPCA(alldata.integrated, pcs = 30, verbose = FALSE, reduction.name =
                               "PCA_on_CCA")
alldata.integrated <- RunUMAP(alldata.integrated, reduction = "PCA_on_CCA", dims = 1:25,
                              reduction.name = "UMAP_on_CCA")
alldata.integrated <- RunTSNE(alldata.integrated, reduction = "PCA_on_CCA", dims = 1:25,
                              reduction.name = "TSNE_on_CCA")
png("UMAP_SCT_on_CCA.png", units="in", width=15, height=15, res=300)
DimPlot(alldata.integrated, reduction = "UMAP_on_CCA", group.by = "orig.ident")
dev.off()
png("TSNE_SCT_on_CCA.png", units="in", width=15, height=15, res=300)
DimPlot(alldata.integrated, reduction = "TSNE_on_CCA", group.by = "orig.ident")
dev.off()
saveRDS(alldata.integrated , "alldata.integrated.rds")
save.image("alldata.integrated.RData")
##### Tcell Annotation
library(ProjecTILs)
cd4 <- readRDS("/Users/merajam/Documents/Tools/ProjecTILs/CD4T_human_ref_v1.rds")
cd8 <- readRDS("/Users/merajam/Documents/Tools/ProjecTILs/CD8T_human_ref_v1.rds")
getwd()
load("Tcell_clustering.RData")
ls()
rm(list=c("alldata", "Rtcells.0.8_markers","Rtcells.1.5", "Rtcells.1.5_markers", "Rtcells.2.0"))
ls()
rm("Rtcells.2.0_markers" )
ls()
Rtcells.cr
Rtcells.0.8
cd4.projected <- Run.ProjecTILs(Rtcells.0.8, ref=cd4)
cd4.projected
cd8.projected <- Run.ProjecTILs(Rtcells.0.8, ref=cd8)
cd8.projected
Tcells <- Rtcells.0.8
cd8.projected_no_filter <- Run.ProjecTILs(Rtcells.0.8, ref=cd8, filter.cells = FALSE)
cd4.projected_no_filter <- Run.ProjecTILs(Rtcells.0.8, ref=cd4, filter.cells = FALSE)
plot.statepred.composition(cd4, cd4.projected, metric = "Percent")
plot.statepred.composition(cd4, cd4.projected_no_filter, metric = "Percent")
plot.statepred.composition(cd4, cd4.projected, metric = "Percent")
plot.statepred.composition(cd4, cd4.projected_no_filter, metric = "Percent")
plot.statepred.composition(cd4, cd4.projected, metric = "Percent")
plot.statepred.composition(cd4, cd4.projected_no_filter, metric = "Percent")
cd4.projected
table(cd4.projected$functional.cluster)
table(cd4.projected_no_filter$functional.cluster)
table(cd8.projected_no_filter$functional.cluster)
table(cd8.projected$functional.cluster)
plot.states.radar(cd4, query = cd4.projected, min.cells = 30)
table(cd8.projected$functional.cluster, cd8.projected$seurat_clusters)
table(cd8.projected_no_filter$functional.cluster, cd8.projected_no_filter$seurat_clusters)
table(cd4.projected_no_filter$functional.cluster, cd4.projected_no_filter$seurat_clusters)
*Rtcells.0.8$seurat_clusters)
table(Rtcells.0.8$seurat_clusters)
barplot(cd4.projected_no_filter$functional.cluster, cd4.projected_no_filter$seurat_clusters)
plot(cd4.projected_no_filter$functional.cluster, cd4.projected_no_filter$seurat_clusters)
plot(cd4.projected_no_filter$functional.cluster, cd4.projected_no_filter$seurat_clusters,
     na.rm=T)
table(cd4.projected_no_filter$functional.cluster, cd4.projected_no_filter$seurat_clusters)
boxplot(cd4.projected_no_filter$functional.cluster, cd4.projected_no_filter$seurat_clusters)
table(cd4.projected_no_filter$functional.cluster, cd4.projected_no_filter$seurat_clusters)
barplot(table(cd4.projected_no_filter$functional.cluster,
              cd4.projected_no_filter$seurat_clusters),na.rm=T)
barplot(table(cd4.projected_no_filter$functional.cluster,
              cd4.projected_no_filter$seurat_clusters))
barplot(table(cd4.projected_no_filter$functional.cluster,
              cd4.projected_no_filter$seurat_clusters))
barplot(table(cd4.projected_no_filter$functional.cluster,
              cd4.projected_no_filter$seurat_clusters))
barplot(table(cd4.projected_no_filter$functional.cluster,
              cd4.projected_no_filter$seurat_clusters), col = c("red", "green", "purple", "blue",
                                                                "yellow","brown","pink")))
barplot(table(cd4.projected_no_filter$functional.cluster,
              cd4.projected_no_filter$seurat_clusters), col = c("red", "green", "purple", "blue",
                                                                "yellow","brown","pink"))
barplot(table(cd4.projected$functional.cluster, cd4.projected$seurat_clusters), col = c("red",
                                                                                        "green", "purple", "blue", "yellow","brown","pink"))
barplot(table(cd4.projected_no_filter$functional.cluster,
              cd4.projected_no_filter$seurat_clusters), col = c("red", "green", "purple", "blue",
                                                                "yellow","brown","pink"))
barplot(table(cd8.projected_no_filter$functional.cluster,
              cd8.projected_no_filter$seurat_clusters), col = c("red", "green", "purple", "blue",
                                                                "yellow","brown","pink"))
table(cd8.projected_no_filter$functional.cluster, cd8.projected_no_filter$seurat_clusters)
table(cd8.projected$functional.cluster, cd8.projected$seurat_clusters)
table(cd8.projected$functional.cluster, cd8.projected$seurat_clusters)
table(cd8.projected_no_filter$functional.cluster, cd8.projected_no_filter$seurat_clusters)
dev.off()
dev.off()
table(cd8.projected_no_filter$functional.cluster, cd8.projected_no_filter$seurat_clusters)
barplot(table(cd8.projected_no_filter$functional.cluster,
              cd8.projected_no_filter$seurat_clusters), col = c("red", "green", "purple", "blue",
                                                                "yellow","brown","pink"))
barplot(table(cd8.projected$functional.cluster, cd8.projected$seurat_clusters), col = c("red",
                                                                                        "green", "purple", "blue", "yellow","brown","pink"))
barplot(table(cd8.projected_no_filter$functional.cluster,
              cd8.projected_no_filter$seurat_clusters), col = c("red", "green", "purple", "blue",
                                                                "yellow","brown","pink"))
barplot(table(cd8.projected$functional.cluster, cd8.projected$seurat_clusters), col = c("red",
                                                                                        "green", "purple", "blue", "yellow","brown","pink"))
barplot(table(cd4.projected$functional.cluster, cd4.projected$seurat_clusters), col = c("red",
                                                                                        "green", "purple", "blue", "yellow","brown","pink"))
barplot(table(cd8.projected$functional.cluster, cd8.projected$seurat_clusters), col = c("red",
                                                                                        "green", "purple", "blue", "yellow","brown","pink"))
barplot(table(cd8.projected$functional.cluster, cd8.projected$seurat_clusters), col = c("red",
                                                                                        "green", "purple", "blue", "yellow","brown","pink"),
        legend=c("CD8.CM","CD8.EM","CD8.MAIT","CD8.NaiveLike","CD8.TEMRA","CD8.TEX","CD8.TPEX"))
barplot(table(cd8.projected$functional.cluster, cd8.projected$seurat_clusters), col = c("red",
                                                                                        "green", "purple", "blue", "yellow","brown","pink"))
plot.projection(cd4, cd4.projected, linesize = 0.5, pointsize = 0.5)
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000")
plot.projection(cd4, cd4.projected, linesize = 0.5, pointsize = 0.5, cols=refCols)
plot.projection(cd8, cd8.projected, linesize = 0.5, pointsize = 0.5, cols=refCols)
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#e812dd", "#FF0000")
plot.projection(cd8, cd8.projected, linesize = 0.5, pointsize = 0.5, cols=refCols)
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#87f6a5", "#00B6EB", "#e812dd", "#FF0000")
plot.projection(cd8, cd8.projected, linesize = 0.5, pointsize = 0.5, cols=refCols)
plot.projection(cd4, cd4.projected, linesize = 0.5, pointsize = 0.5, cols=refCols)
barplot(table(cd8.projected$functional.cluster, cd8.projected$seurat_clusters), col = c("red",
                                                                                        "green", "purple", "blue", "yellow","brown","pink"))
barplot(table(cd4.projected$functional.cluster, cd4.projected$seurat_clusters), col = c("red",
                                                                                        "green", "purple", "blue", "yellow","brown","pink"))
table_cd4_nofilter <- table(cd8.projected_no_filter$functional.cluster,
                            cd8.projected_no_filter$seurat_clusters)
table_cd4_nofilter
barplot( table_cd4_nofilter,col = c("red", "green", "purple", "blue", "yellow","brown","pink"),
         legend = rownames(table_cd4_nofilter))
barplot( table_cd4_nofilter,col = refcols, legend = rownames(table_cd4_nofilter))
barplot( table_cd4_nofilter,col = refCols, legend = rownames(table_cd4_nofilter))
table_cd8_nofilter <- table(cd8.projected_no_filter$functional.cluster,
                            cd8.projected$seurat_clusters)
table_cd8_nofilter <- table(cd8.projected_no_filter$functional.cluster,
                            cd8.projected_no_filter$seurat_clusters)
table_cd8 <- table(cd8.projected$functional.cluster, cd8.projected$seurat_clusters)
table_cd4 <- table(cd4.projected$functional.cluster, cd4.projected$seurat_clusters)
table_cd4_n <- table(cd4.projected$functional.cluster, cd4.projected$seurat_clusters)
barplot( table_cd4_nofilter,col = refCols, legend = rownames(table_cd4_nofilter))
barplot( table_cd4,col = refcols, legend = rownames(table))
barplot( table_cd4,col = refCols, legend = rownames(table))
barplot( table_cd4,col = refCols, legend = rownames(table_cd4))
