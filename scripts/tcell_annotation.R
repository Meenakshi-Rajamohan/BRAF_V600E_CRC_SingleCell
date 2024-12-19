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
