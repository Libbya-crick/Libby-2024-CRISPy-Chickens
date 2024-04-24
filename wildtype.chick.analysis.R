##Install packages
install.packages("Seurat")
install.packages("dplyr")
install.packages("rjson")
install.packages("viridis")
install.packages("pheatmap")
install.packages("rstatix")
install.packages("dendextend")


##Load libraries
library(Matrix)
library(Seurat)
library(dplyr)
library(patchwork)
library(scales)
library (ggplot2)
library(sctransform)
library(Matrix)
library(dplyr)
library(patchwork)
library(rjson)
library(viridis)
library(RColorBrewer)
library(data.table)
library(pheatmap)
library(rstatix)
library(dendextend)
library(stringr)
library(ggrepel)

## Load RDS files
## in this example we reanalyze the Rito et al. data set (accession number GSE223189), 
##the 10 somite embryo data (GSM6940809) and the 13 somite embryo data (GSM6940810) were processed first in the Cell Ranger count pipeline.
##This will begin by first uploading the two datasets and merging them together

data.dir_13<- 'directory to 13 somite 'filtered_feature_bc_matrix''
data.dir_10<- 'directory to 10 somite 'filtered_feature_bc_matrix''

chick_10S_raw.data <- Read10X(data.dir_10)
chick_13S_raw.data <- Read10X(data.dir_13)

chick_13S_raw <- CreateSeuratObject(counts = chick_13S_raw.data, project = "13S", min.cells = 3, min.features = 200)
chick_10S_raw <- CreateSeuratObject(counts = chick_10S_raw.data, project = "10S", min.cells = 3, min.features = 200)

chick_10_13S_raw <- merge(x = chick_10S_raw, y = chick_13S_raw, add.cell.ids = c("10S", "13S"), project = "combined")
chick_10_13S_raw

chick_10_13S_raw<- JoinLayers(object = chick_10_13S_raw, layer = "counts")
chick_10_13S_raw

##enter gene lists to be used to help filter raw data sets
mit_genes <- c("ATP6", "ATP8", "COII", "COX3", "CYTB", "MT-CO1", "MT-ND2","ND1", "ND3","ND4","ND4L","ND5","ND6")
W_genes <- c("ENSGALG00000014184","ENSGALG00000026991","ENSGALG00000027170","ENSGALG00000030382","ENSGALG00000030707","ENSGALG00000030997","ENSGALG00000031327","ENSGALG00000031355",
             "ENSGALG00000033705","ENSGALG00000033827","ENSGALG00000034125","ENSGALG00000034488","ENSGALG00000034905","ENSGALG00000035780","ENSGALG00000035785","ENSGALG00000038064",
             "ENSGALG00000039023","ENSGALG00000040263","ENSGALG00000041221","ENSGALG00000041500","ENSGALG00000043758","ENSGALG00000043838","ENSGALG00000044674","ENSGALG00000044734",
             "ENSGALG00000045131","ENSGALG00000045335","ENSGALG00000045984","ENSGALG00000046681","ENSGALG00000046757","ENSGALG00000046789","ENSGALG00000046864","ENSGALG00000046886",
             "ENSGALG00000046966","ENSGALG00000047164","ENSGALG00000047355","ENSGALG00000047426","ENSGALG00000047434","ENSGALG00000047588","ENSGALG00000047682","ENSGALG00000047839",
             "ENSGALG00000047904","ENSGALG00000048170","ENSGALG00000048286","ENSGALG00000048307","ENSGALG00000048432","ENSGALG00000048542","ENSGALG00000048633","ENSGALG00000048665",
             "ENSGALG00000048839","ENSGALG00000048969","ENSGALG00000049069","ENSGALG00000049073","ENSGALG00000049167","ENSGALG00000049210","ENSGALG00000049219","ENSGALG00000049305",
             "ENSGALG00000049351","ENSGALG00000049534","ENSGALG00000049732","ENSGALG00000049785","ENSGALG00000049821","ENSGALG00000049970","ENSGALG00000050141","ENSGALG00000050274",
             "ENSGALG00000050349","ENSGALG00000050356","ENSGALG00000050474","ENSGALG00000050520","ENSGALG00000050647","ENSGALG00000050817","ENSGALG00000050879","ENSGALG00000050881",
             "ENSGALG00000050905","ENSGALG00000050935","ENSGALG00000051029","ENSGALG00000051123","ENSGALG00000051197","ENSGALG00000051237","ENSGALG00000051283","ENSGALG00000051422",
             "ENSGALG00000051455","ENSGALG00000051508","ENSGALG00000051568","ENSGALG00000051624","ENSGALG00000051913","ENSGALG00000052041","ENSGALG00000052395","ENSGALG00000052405",
             "ENSGALG00000052623","ENSGALG00000052752","ENSGALG00000052920","ENSGALG00000052946","ENSGALG00000052997","ENSGALG00000053189","ENSGALG00000053504","ENSGALG00000053521",
             "ENSGALG00000053584","ENSGALG00000053798","ENSGALG00000053841","ENSGALG00000053933","ENSGALG00000054122","ENSGALG00000054257","ENSGALG00000054371","ENSGALG00000054501",
             "ENSGALG00000054546","ENSGALG00000054658","ENSGALG00000054809","ENSGALG00000054870","ENSGALG00000054907","ENSGALG00000054989","ENSGALG00000055012","gga-mir-122-2","gga-mir-7b",
             "HINTW","HNRNPKL","SMAD7B","SNORD58","SPIN1W", "U6", "UBAP2")

## first because these data sets contained multiple pooled embryos, 
## find genes associated with the W-sex chromosome to filter out so they don't skew clustering

W_genes_2 <- W_genes[W_genes %in% row.names(chick_10_13S_raw@assays$RNA)]
counts <- GetAssayData(chick_10_13S_raw, assay = "RNA")
counts <- counts[-(which(row.names(counts) %in% W_genes_2)),]
chick_1013_dW <- subset(chick_10_13S_raw, features = row.names(counts))

##then filter and visualize on mitochondrial genes

chick_1013_dW[["percent.mt"]] <- PercentageFeatureSet(chick_1013_dW, features = mit_genes)

plot1 <- VlnPlot(chick_1013_dW, features = c("nCount_RNA", "nFeature_RNA","percent.mt"), ncol = 3 )
plot2 <- FeatureScatter(chick_1013_dW, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
options(repr.plot.width=20, repr.plot.height=10)
plot1+plot2

chick_1013_dW <- subset(chick_1013_dW, subset = nCount_RNA > 2500 & nFeature_RNA >300 & nFeature_RNA < 3500 & percent.mt< 10 & percent.mt>.5)

##Lof normalize and use the standard Seurat workflow to find variable features in dataset

chick_1013_dW <- NormalizeData(chick_1013_dW, normalization.method = "LogNormalize", scale.factor = 10000)
chick_1013_dW <- FindVariableFeatures(chick_1013_dW,selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(chick_1013_dW), 10)
top20 <- head(VariableFeatures(chick_1013_dW), 20)
plot1 <- VariableFeaturePlot(chick_1013_dW, pt.size = 3)
plot2 <- LabelPoints(plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(chick_1013_dW)
chick_1013_dW <- ScaleData(chick_1013_dW, features = all.genes)
chick_1013_dW <- RunPCA(chick_1013_dW, features = VariableFeatures(object = chick_1013_dW))

VizDimLoadings(chick_1013_dW, dims = 1:2, reduction = "pca")
DimPlot(chick_1013_dW, pt.size = 1.5, reduction = "pca")
DimHeatmap(chick_1013_dW, dims = 1:3, cells = 500, balanced = TRUE)

chick_1013_dW <- JackStraw(chick_1013_dW, num.replicate = 100)
chick_1013_dW <- ScoreJackStraw(chick_1013_dW, dims = 1:20)
JackStrawPlot(chick_1013_dW, dims = 1:20)
ElbowPlot(chick_1013_dW, ndims = 50, reduction = "pca" )

## Find UMAP clusters based on Variable features (Figure 1B - left)
chick_1013_dW <- FindNeighbors(chick_1013_dW, dims = 1:40)
chick_1013_dW <- FindClusters(chick_1013_dW, resolution = 0.9)
chick_1013_dW <- RunUMAP(chick_1013_dW, dims = 1:20)
plot3 <- DimPlot(chick_1013S_dW, reduction = "umap",pt.size = 1.5, label = TRUE)
plot3

plot4 <- DimPlot(chick_1013_dW, reduction = 'umap', group.by = "orig.ident")
plot4

## save the seurat object
saveRDS(chick_1013_dW, file = "file_out_path.rds")

##plot key genes of interest to mark primitive streak and neural tube
p_favorite_genes <-FeaturePlot(chick_1013_dW, features = c( "NOTO","TBXT", "SOX2", "PAX6", "SOX10", "FOXA2"), pt.size = 1)
options(repr.plot.width=20, repr.plot.height=10)
p_favorite_genes

## To then isolate the lineages from paraxial mesoderm to neural tube subset the data based on clusters 
##that expressed the markers in p_favorite_genes, excluding notochord (high TBXT and NOTO)

all.genes <- rownames(chick_1013_dW)
PRENEURAL <- subset(chick_1013_dW, idents = c(1,4,7,8,15))
PRENEURAL <- ScaleData(PRENEURAL, features = all.genes)
PRENEURAL <- FindVariableFeatures(PRENEURAL,
                                  selection.method = "vst")

##Then load Entropy ranked gene list and use this to cluster the subsetting data

entropy.ranked.genes<- read.csv(file = "ES_Cluster_Ranked_Genes_Lists.csv" )
ES_500<- entropy.ranked.genes[c(1:500),-1]
ES_500<- unique(c(ES_500$X1,ES_500$X10,ES_500$X10,ES_500$X14,
                  ES_500$X4,ES_500$X7,ES_500$X8,ES_500$Mesoderm,ES_500$PS ))

PRENEURAL <- RunPCA(PRENEURAL, features =ES_500)
VizDimLoadings(PRENEURAL, dims = 1:2, reduction = "pca")
DimPlot(PRENEURAL, pt.size = 1.5, reduction = "pca")
DimHeatmap(PRENEURAL, dims = 1:3, cells = 500, balanced = TRUE)

PRENEURAL <- JackStraw(PRENEURAL, num.replicate = 100)
PRENEURAL <- ScoreJackStraw(PRENEURAL, dims = 1:20)
JackStrawPlot(PRENEURAL, dims = 1:20)
ElbowPlot(PRENEURAL, ndims = 50, reduction = "pca" )

PRENEURAL <- FindNeighbors(PRENEURAL, dims = 1:20)
PRENEURAL <- FindClusters(PRENEURAL, resolution = 0.9)
PRENEURAL <- RunUMAP(PRENEURAL, dims = 1:20)

## plot UMAP of subsetted dataset (Figure 1B - right)
Idents(PRENEURAL)<-PRENEURAL@meta.data$seurat_clusters
p<-DimPlot(PRENEURAL, reduction = 'umap', pt.size = 4,label = FALSE ) 
p

##plot genes of interest on this dataset (Figure 1D)

Idents(PRENEURAL)<-  PRENEURAL$seurat_clusters
my_levels <- c("0","1","2","6","5","8","7","4","3","9","10","11")
Idents(PRENEURAL)<-factor(Idents(PRENEURAL), levels = my_levels)
p<- DotPlot(PRENEURAL,  features = c("FOXC2","MSGN1","PARAXIS","TBX6", "TBXT","FGF8","ID4","CDX2","BMP4", 
                                     "NOTO","PAX7","SOX2","NKX1-2",
                                     "SOX21","SOX1",  "PAX6","NKX6-2","OLIG2",
                                     "SOX9","SOX10","FOXD3","WNT6","SHH","FOXA2",  
                                     "NEUROG1","ONECUT1"),
            cols = c("grey","blue"),dot.scale = 17)+
  theme(axis.text.x = element_text(angle = 90))

p

## We plot the expression of known lineage markers (Figure S1C)

blend <- TRUE
x<- "SOX10" ## TBX6, CDX2, PAX7, MSGN1, IRX3, FOXC2
y<- "FOXA2" ## SOX21, PAX6, NKX6-2, NKX1-2, OLIG2, NEUROG1
p <- FeaturePlot(PRENEURAL, reduction = "umap", features = c(x,y),  blend = TRUE, blend.threshold = 0.1,
                 cols = if (blend) {c("lightgrey", "blue", "red")} else {c("lightgrey", "blue")},  
                 pt.size = 2,  combine = FALSE)
p <- lapply(p, function(x) x + NoLegend())
p<- CombinePlots(p)
p

## We additionally plot the expression of SOX2, SOX3, and TBXT overlap to identify neuromesodermal progenitors (Figure S1C)

NMPs2.exp<- PRENEURAL@assays$RNA$data["SOX2",]
NMPs3.exp<- PRENEURAL@assays$RNA$data["SOX3",]
NMPst.exp<- PRENEURAL@assays$RNA$data["TBXT",]

TBXT <- which(NMPt.exp>0 & NMPs2.exp==0 & NMPs3.exp==0)
SOX2 <- which(NMPt.exp==0 & NMPs2.exp>0 & NMPs3.exp==0)
SOX3 <- which(NMPt.exp==0 & NMPs2.exp==0 & NMPs3.exp>0)
ALL <- which(NMPt.exp>0 & NMPs2.exp>0 & NMPs3.exp>0)
TBXTSOX2 <- which(NMPt.exp>0 & NMPs2.exp>0 & NMPs3.exp==0)
TBXTSOX3 <- which(NMPt.exp>0 & NMPs2.exp==0 & NMPs3.exp>0)
SOX2SOX3 <- which(NMPt.exp==0 & NMPs2.exp>0 & NMPs3.exp>0)
NEG <- which(NMPt.exp==0 & NMPs2.exp==0 & NMPs3.exp==0)

cell_ids <- rep(0,dim(PRENEURAL@meta.data)[1])

cell_ids[TBXT] <- 'TBXT+'
cell_ids[SOX2] <- 'SOX2+'
cell_ids[SOX3] <- 'SOX3+'
cell_ids[TBXTSOX2] <- 'TBXT+/SOX2+'
cell_ids[TBXTSOX3] <- 'TBXT+/SOX3+'
cell_ids[SOX2SOX3] <- 'SOX2+/SOX3+'
cell_ids[ALL] <- 'ALL+'
cell_ids[NEG] <- 'Negative'

PRENEURAL@meta.data['expression_type'] <- cell_ids

options(repr.plot.width=20, repr.plot.height=15)
p<- DimPlot(PRENEURAL, reduction = "umap", dims=c(1,2),
            group.by = 'expression_type', cols = c("purple","grey","red", "red","red","dodgerblue", "purple","purple" ), pt.size = 6)
p


## Later in paper we look into the Super Elongation Complex and plot the expression of associated genes (Figure 5D)

SEC<- c("AFF1","AFF2","AFF3","AFF4","TAT",
        "NELFA","NELFB","NELFCD","NELFE",
        "POLA1","MLLT1","ENSGALG00010002919","BRD4L",
        "MLLT3","ELL","ELL2","ELL3","CDK9","CCNT1",
        "GTF2F2","ELL",
        "BRD4","JMJD6", "EAR", "PAF1")
##CDK9 <- pTEFb
##MLLT1 <- ENL
##NELFa-E<- NELF complex
## BRD4 <- BRD4L = "ENSGALG00010002919"

Idents(PRENEURAL) <- PRENEURAL$seurat_clusters
my_levels <- c(1,2,3,7,6,9,8,5,4,10,11,12)

matrix<- AverageExpression(PRENEURAL, features = SEC)
matrix<- as.data.frame(matrix)
matrix<- matrix[,my_levels]
matrix<-scale(matrix)
matrix

plot<- pheatmap(matrix, color = viridis(25),
                cluster_cols = FALSE)
plot
