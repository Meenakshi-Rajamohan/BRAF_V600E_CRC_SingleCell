library(Seurat)
library(Matrix)
library(DoubletFinder)
library(knitr)
library(cowplot)
library(ggplot2)
library(scran)
setwd("/Users/merajam/Documents/CRC/2_Analysis/1_QC/BRAF/")
alldata <- readRDS("/Users/merajam/Documents/CRC/1_Data.rds")
total_counts_per_cell <- colSums(alldata@assays$RNA@counts)
mito_genes <- rownames(alldata)[grep("^MT-", rownames(alldata))]
alldata$percent_mito <- colSums(alldata@assays$RNA@counts[mito_genes, ])/total_counts_per_cell
head(mito_genes, 10)
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
saveRDS(data.filt, "CRC_BRAF_afterQC_BRAF.rds")
save.image("CRC_BRAF_QC.RData")
