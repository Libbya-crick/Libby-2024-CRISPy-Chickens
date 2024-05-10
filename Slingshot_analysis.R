##Install packages
install.packages("Seurat")
install.packages("viridis")
install.packages("pheatmap")
install.packages("slingshot")
install.packages("tidymodels")



library(Matrix)
library(Seurat)
library(dplyr)
library(patchwork)
library(scales)
library (ggplot2)
library(sctransform)
suppressPackageStartupMessages({
  library(scran)
  library(scater)
  library(igraph)
})
library(tidymodels)
library(slingshot)
library(viridis)
library(tradeSeq)
library(RColorBrewer)
library(data.table)
library(pheatmap)
library(dendextend)

## load wildtype dataset

PRENEURAL <- readRDS(file = "PRENEURAL_ES.500.rds")

## generate lineages through dataset
data <- FindVariableFeatures(PRENEURAL, nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data)
data <- FindNeighbors(data)
data <- FindClusters(data, resolution = 0.4)
data <- RunUMAP(data, n.neighbors = 10, dims = 1:50, spread = 2, min.dist = 0.3)

# Save the objects as separate matrices for input in slingshot
dimred <- PRENEURAL@reductions$umap@cell.embeddings
dimred

clustering <- data$RNA_snn_res.0.4
#clustering

counts <- as.matrix(data[["RNA"]]$counts)

# Run default Slingshot lineage identification
set.seed(1)
lineages <- getLineages(data = dimred, clusterLabels = clustering)
lineages


#Run default Slingshot
set.seed(1)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        #end.clus = c("2,10,4"), #define how many branches/lineages to consider
                        start.clus = c( "1")) #define where to start the trajectories

lineages

lin_sling<-SlingshotDataSet(lineages)
lin_sling


# Plot the lineages
options(repr.plot.width=15, repr.plot.height=10)
par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = cell_colors_clust, cex = 0.5, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = cell_colors_clust, cex = 0.5, pch = 16)
lines(lin_sling, lwd = 3, col = "black")

# generate curves
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves
curves_sling<-SlingshotDataSet(curves)

#plot curves
p<-plot(dimred, col = cell_colors_clust, asp = 1, pch = 16)
d<- lines(curves_sling, lwd = 3, col = "black")
d

#plot multiple lineages
nc <- 3
pt <- slingPseudotime(curves_sling)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(curves_sling), col = colors, pch = 16, cex = 0.5, main = i)
  lines(curves_sling, lwd = 3, col = "black",)
}

#fit gene counts to curves
sce <- fitGAM(counts = as.matrix(counts), sds = curves_sling)
options(repr.plot.width=15, repr.plot.height=10)
plotGeneCount(curves_sling, counts, clusters = clustering, models = sce)

#plot genes of interest along lineages

gene_list <- c("MSGN1","TBX6",
               "ID2","MSX2", "TBXT","MSX1","MLLT3","MNX1","SPRY2","CDX2",
               "CMTM8","CDX4","F2RL1","CLDN1","NKX1-2","PAX7",#"MAPK3",
               "SOX2","SOX21", "FGFR1", "FABP5","SOX1",
               "PAX6", "FGFR2", "OLIG2","FGFR3",
               "NKX6-2","FOXA2","SOX9")
yhatSmooth <- predictSmooth(sce, gene = gene_list, nPoints = 150, tidy = FALSE)
lin1<-yhatSmooth[,151:300]
heatSmooth <- pheatmap(t(scale(t(lin2))),
                       cluster_cols = FALSE,
                       cluster_rows = TRUE,
                       kmeans_k = NA,
                       #revC = TRUE,
                       #cutree_rows = 8,
                       viridis(100),
                       show_rownames = TRUE, 
                       show_colnames = FALSE)
heatSmooth