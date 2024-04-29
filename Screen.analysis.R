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

## in this analysis pipeline we first filter the raw screen single cell RNA sequencing result
## this will remove contaminating populations and map the resultant cells back onto the wild type dataset (Figure 1)
## analysis of wildtype dataset is described in "wildtype.chick.analysis.R" in theis repository

## first load four screen sequencing datasets and the wildtype dataset

chick_1013S_dW <- readRDS(file = "path to wildtype.rds")


data.dir_1 = 'directory to sample 1 'filtered_feature_bc_matrix''
data.dir_2 = 'directory to sample 2 'filtered_feature_bc_matrix''
data.dir_3 = 'directory to sample 3 'filtered_feature_bc_matrix''
data.dir_4 = 'directory to sample 4 'filtered_feature_bc_matrix''

sample1.data = Read10X(data.dir_1)
sample2.data = Read10X(data.dir_2)
sample3.data = Read10X(data.dir_3)
sample4.data = Read10X(data.dir_4)

##Create and merge seurat objects of each sample

sample1 <- CreateSeuratObject(counts = sample1.data$'Gene Expression', project = "sample 1", min.cells = 3, min.features = 200)
sample2 <- CreateSeuratObject(counts = sample2.data$'Gene Expression', project = "sample 2", min.cells = 3, min.features = 200)
sample3 <- CreateSeuratObject(counts = sample3.data$'Gene Expression', project = "sample 3", min.cells = 3, min.features = 200)
sample4 <- CreateSeuratObject(counts = sample4.data$'Gene Expression', project = "sample 4", min.cells = 3, min.features = 200)

screen_seq_ori <- merge(sample1, y = c(sample2,sample3,sample4), project = "screen_seq")
screen_seq_ori

## Now we filter out W (sex chromosome) genes, cells with high mitochondrial genes, cells with low RNA count or low feature count

W_genes <- c("ENSGALG00000014184","ENSGALG00000026991","ENSGALG00000027170","ENSGALG00000030382","ENSGALG00000030707","ENSGALG00000030997","ENSGALG00000031327","ENSGALG00000031355","ENSGALG00000033705","ENSGALG00000033827","ENSGALG00000034125","ENSGALG00000034488","ENSGALG00000034905","ENSGALG00000035780","ENSGALG00000035785","ENSGALG00000038064","ENSGALG00000039023","ENSGALG00000040263","ENSGALG00000041221","ENSGALG00000041500","ENSGALG00000043758","ENSGALG00000043838","ENSGALG00000044674","ENSGALG00000044734","ENSGALG00000045131","ENSGALG00000045335","ENSGALG00000045984","ENSGALG00000046681","ENSGALG00000046757","ENSGALG00000046789","ENSGALG00000046864","ENSGALG00000046886","ENSGALG00000046966","ENSGALG00000047164","ENSGALG00000047355","ENSGALG00000047426","ENSGALG00000047434","ENSGALG00000047588","ENSGALG00000047682","ENSGALG00000047839","ENSGALG00000047904","ENSGALG00000048170","ENSGALG00000048286","ENSGALG00000048307","ENSGALG00000048432","ENSGALG00000048542","ENSGALG00000048633","ENSGALG00000048665","ENSGALG00000048839","ENSGALG00000048969","ENSGALG00000049069","ENSGALG00000049073","ENSGALG00000049167","ENSGALG00000049210","ENSGALG00000049219","ENSGALG00000049305","ENSGALG00000049351","ENSGALG00000049534","ENSGALG00000049732","ENSGALG00000049785","ENSGALG00000049821","ENSGALG00000049970","ENSGALG00000050141","ENSGALG00000050274","ENSGALG00000050349","ENSGALG00000050356","ENSGALG00000050474","ENSGALG00000050520","ENSGALG00000050647","ENSGALG00000050817","ENSGALG00000050879","ENSGALG00000050881","ENSGALG00000050905","ENSGALG00000050935","ENSGALG00000051029","ENSGALG00000051123","ENSGALG00000051197","ENSGALG00000051237","ENSGALG00000051283","ENSGALG00000051422","ENSGALG00000051455","ENSGALG00000051508","ENSGALG00000051568","ENSGALG00000051624","ENSGALG00000051913","ENSGALG00000052041","ENSGALG00000052395","ENSGALG00000052405","ENSGALG00000052623","ENSGALG00000052752","ENSGALG00000052920","ENSGALG00000052946","ENSGALG00000052997","ENSGALG00000053189","ENSGALG00000053504","ENSGALG00000053521","ENSGALG00000053584","ENSGALG00000053798","ENSGALG00000053841","ENSGALG00000053933","ENSGALG00000054122","ENSGALG00000054257","ENSGALG00000054371","ENSGALG00000054501","ENSGALG00000054546","ENSGALG00000054658","ENSGALG00000054809","ENSGALG00000054870","ENSGALG00000054907","ENSGALG00000054989","ENSGALG00000055012","gga-mir-122-2","gga-mir-7b","HINTW","HNRNPKL","SMAD7B","SNORD58","SPIN1W", "U6", "UBAP2")
W_genes_2 <- W_genes[W_genes %in% row.names(screen_seq_ori@assays$RNA)]

counts <- GetAssayData(screen_seq_ori, assay = "RNA")
counts <- counts[-(which(row.names(counts) %in% W_genes_2)),]
screen_seq_ori <- subset(screen_seq_ori, features = row.names(counts))

mit_genes <- c("ATP6", "ATP8", "COII", "COX3", "CYTB", "MT-CO1", "MT-ND2","ND1", "ND3","ND4","ND4L","ND5","ND6")
screen_seq_ori[["percent.mt"]] <- PercentageFeatureSet(screen_seq_ori, features = mit_genes)

plot1 <- VlnPlot(screen_seq_ori, features = c("nCount_RNA", "nFeature_RNA","percent.mt"), ncol = 3 )
plot1
plot2 <- FeatureScatter(screen_seq_ori, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
options(repr.plot.width=20, repr.plot.height=10)
plot1+plot2

screen_seq <- subset(screen_seq_ori, subset = nCount_RNA > 2500 & nFeature_RNA >300 & nFeature_RNA < 5500 & percent.mt< 20 & percent.mt>.5)
screen_seq <- NormalizeData(screen_seq, normalization.method = "LogNormalize", scale.factor = 10000)

#now we follow the Seurat clustering pipeline from the Satija lab

screen_seq <- FindVariableFeatures(screen_seq,selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(screen_seq), 10)
top20 <- head(VariableFeatures(screen_seq), 20)
plot1 <- VariableFeaturePlot(screen_seq, pt.size = 3)
plot2 <- LabelPoints(plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(screen_seq)
screen_seq <- ScaleData(screen_seq, features = all.genes)
screen_seq <- RunPCA(screen_seq, features = VariableFeatures(object = screen_seq))

VizDimLoadings(screen_seq, dims = 1:2, reduction = "pca")
DimPlot(screen_seq, pt.size = 1.5, reduction = "pca")
DimHeatmap(screen_seq, dims = 1:3, cells = 500, balanced = TRUE)

screen_seq <- JackStraw(screen_seq, num.replicate = 100)
screen_seq <- ScoreJackStraw(screen_seq, dims = 1:20)
JackStrawPlot(screen_seq, dims = 1:20)
ElbowPlot(screen_seq, ndims = 70, reduction = "pca" )

#here we added in a step to determine clustering resolution by creating an elbow plot of inertias between clusters

screen_seq <- FindNeighbors(screen_seq, dims = 1:20)
res <- c(0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0)
for(x in res){screen_seq <- FindClusters(screen_seq, resolution = x)}

res_0.2 <- screen_seq@meta.data['RNA_snn_res.0.2']
res_0.4 <- screen_seq@meta.data['RNA_snn_res.0.4']
res_0.6 <- screen_seq@meta.data['RNA_snn_res.0.6']
res_0.8 <- screen_seq@meta.data['RNA_snn_res.0.8']
res_1.0 <- screen_seq@meta.data['RNA_snn_res.1']
res_1.2 <- screen_seq@meta.data['RNA_snn_res.1.2']
res_1.4 <- screen_seq@meta.data['RNA_snn_res.1.4']
res_1.6 <- screen_seq@meta.data['RNA_snn_res.1.6']
res_1.8 <- screen_seq@meta.data['RNA_snn_res.1.8']
res_2.0 <- screen_seq@meta.data['RNA_snn_res.2']
res_2.2 <- screen_seq@meta.data['RNA_snn_res.2.2']
res_2.4 <- screen_seq@meta.data['RNA_snn_res.2.4']
res_2.6 <- screen_seq@meta.data['RNA_snn_res.2.6']
res_2.8 <- screen_seq@meta.data['RNA_snn_res.2.8']
res_3.0 <- screen_seq@meta.data['RNA_snn_res.3']

mutli_res <- data.frame(c(res_0.2,res_0.4,res_0.6,res_0.8,res_1.0,res_1.2,res_1.4,res_1.6,res_1.8,res_2.0,res_2.2,res_2.4,res_2.6,res_2.8, res_3.0))
rownames(mutli_res)<- colnames(screen_seq)
mutli_res

#distances with clusters squared and then summed. then take average across all cluster sums and this is score,
#want a low score which means points are close to centriod
#then plot inertia score against resolution

min <- which( colnames(screen_seq@meta.data)=="RNA_snn_res.0.2" )
min
max <- which( colnames(screen_seq@meta.data)=="RNA_snn_res.3" )
max
no <- which( colnames(screen_seq@meta.data)=="seurat_clusters" )
no
res <- seq(min, max, 1)
res <- res[res != no ]
res

euclidean <- function(a,b)sqrt(rowSums((a-t(b))^2))

Df_res <- lapply(res, function(r){
  clusters <- c(0:last(levels(screen_seq@meta.data[,r])))
  df <- screen_seq@meta.data[r]
  
  inertias <- lapply(clusters, function(c){
    
    b <- GetAssayData(screen_seq[which(df == c)])
    tot <- ncol(b)
    a <- colSums(b)/tot
    inert <- euclidean(a,b)
    #print (inert)
    inertias <- sum(inert)/length(inert)
    inertias
    
  })
  
  inertias_final <- unlist(inertias)
  inertias_fin <- data.frame(avg_inertia=mean(inertias_final))
  inertias_fin$resolution <-  r
  inertias_fin
  
})

Res_table_inertias <- do.call(rbind,Df_res)

#label with resolutions for plotting
Res_table_inertias[1,2] <- 0.2
Res_table_inertias[2,2] <- 0.4
Res_table_inertias[3,2] <- 0.6
Res_table_inertias[4,2] <- 0.8
Res_table_inertias[5,2] <- 1.0
Res_table_inertias[6,2] <- 1.2
Res_table_inertias[7,2] <- 1.4
Res_table_inertias[8,2] <- 1.6
Res_table_inertias[9,2] <- 1.8
Res_table_inertias[10,2] <- 2.0
Res_table_inertias[11,2] <- 2.2
Res_table_inertias[12,2] <- 2.4
Res_table_inertias[13,2] <- 2.6
Res_table_inertias[14,2] <- 2.8
Res_table_inertias[15,2] <- 3.0

Res_table_inertias
ggplot(Res_table_inertias, aes(x=resolution, y=avg_inertia)) + geom_point(size = 5, col = "blue")
#select the resolution at the corner of the elbow plot
badres <- c('RNA_snn_res.0.2','RNA_snn_res.0.4','RNA_snn_res.0.6',
            'RNA_snn_res.0.8','RNA_snn_res.1.2',
            'RNA_snn_res.1.4','RNA_snn_res.1.6','RNA_snn_res.1.8',
            'RNA_snn_res.2','RNA_snn_res.2.2','RNA_snn_res.2.4',
            'RNA_snn_res.2.6','RNA_snn_res.2.8','RNA_snn_res.3')
good <- 'RNA_snn_res.1'
bres <- badres[badres != good ]
bres

screen_seq@meta.data <- screen_seq@meta.data[,!(names(screen_seq@meta.data) %in% bres)]
names(screen_seq@meta.data)
screen_seq@meta.data['seurat_clusters'] <- screen_seq@meta.data['RNA_snn_res.1'] 
screen_seq <- SetIdent(screen_seq, value = screen_seq@meta.data$seurat_clusters)

#now continue to UMAP plotting
screen_seq <- RunUMAP(screen_seq, dims = 1:20)
plot3 <- DimPlot(screen_seq, reduction = "umap",label = TRUE,label.col = "white", pt.size = 1.5)
plot3

## we then filtered out the contaminating blood from the dataset by plotting GATA2 and LMO2 expression and subsetting out the clusters with high dual expression

FeaturePlot(screen_seq, features = c( "LMO2", "GATA2"), pt.size = 1)
all.genes <- rownames(screen_seq)
screen_seq_no.blood <- subset(screen_seq, idents = c(14,13,12,20,18,23,15,3,4,16,22,7,11))

DimPlot(screen_seq_no.blood, reduction = "umap")

## Add guide RNA counts to the meta.data of the screen_seq_no.blood object
##this is to figure out where the guides are found in the barcodes for cells
## first load found guides and their corresponding cell id barcodes. this spreadsheet matches the cell ids with what barcodes were discovered in each cell

guides1 <- fromJSON(file = 'directory to sample 1 'cells_per_protospacer.json'')
guides2 <- fromJSON(file = 'directory to sample 2 'cells_per_protospacer.json'')
guides3 <- fromJSON(file = 'directory to sample 3 'cells_per_protospacer.json'')
guides4 <- fromJSON(file = 'directory to sample 4 'cells_per_protospacer.json'')

guides <- rbind(guides1, guides2, guides3,guides4)

## also want to figure out how many umis per feature and per cell
## because we amplify barcodes in sequencing prep, we want to eliminate possible ambiant amplification

samples<- c(1:4)
p.spacers<-lapply(samples, function(p){
  file<- paste0(directory to 'protospacer_calls_per_cell.csv')
  p.spacer<- read.csv(file = file)
  cells1<- p.spacer$cell_barcode 
  cell_samples<- lapply(cells1, function(c)
  {cell.sample<- paste0(c,"_",p)                                         
  cell.sample})
  cell_samples<-unlist(cell_samples)
  p.spacer$cell_barcode <- cell_samples
  p.spacer})

p.spacers<- do.call(rbind, p.spacers)
p.spacers

#make list of cells that repeats a cell if it have multiple protospacers
p.cells<- p.spacers$cell_barcode
tot_p.cells<- lapply(p.cells, function(p){x<- p.spacers$num_features[which(p.spacers$cell_barcode == p)]
p.cell<- c(rep(p,x))
p.cell
})
tot_p.cells<- unlist(tot_p.cells)

#make list of number of features for each cell
tot_p.feature<- lapply(p.cells, function(p){x<- p.spacers$num_features[which(p.spacers$cell_barcode == p)]
p.cell<- c(rep(x,x))
p.cell
})
tot_p.feature<- as.integer(unlist(tot_p.feature))


#make list of protospacers in same cell 
tot_p.spacers<-lapply(p.cells, function(p){x<- p.spacers$feature_call[which(p.spacers$cell_barcode == p)]
x<- strsplit(x, "\\|")
x})
tot_p.spacers<-unlist(tot_p.spacers)
#tot_p.spacers

#make list of protospacers counts in same cell
tot_p.counts<-lapply(p.cells, function(p){x<- p.spacers$num_umis[which(p.spacers$cell_barcode == p)]
x<- strsplit(x, "\\|")
x})
tot_p.counts<-as.integer(unlist(tot_p.counts))

# set a thrshold of what is considered a true count by using 2*standard deviation
log_p.counts<- log(tot_p.counts)
range(log_p.counts)
mean<- mean(log_p.counts)
sd<- sd(log_p.counts)
sd
threshold_low <- mean - (2*sd)
threshold_low
bad_low<- which(log_p.counts < threshold_low)

umis_log_thresh<- log_p.counts

umis_log_thresh[bad_low]<- NA


p.spacers_plot<- data.frame(cell_barcode = tot_p.cells,
                            num_features = tot_p.feature,
                            feature = tot_p.spacers,
                            umis = tot_p.counts,
                            umis_log = log_p.counts,
                            umis_log_thresh = umis_log_thresh)

#determine expression of fluorophore Citrine which marks CAS9 expression
#to check if low level fetures are getting Cas9
Cas9G.exp <- GetAssayData(object = screen_seq, assay = "RNA", slot = "data")["CITRINE",]
Cas9G.exp<- unlist(Cas9G.exp)
Cas9G.exp_df<- as.data.frame(Cas9G.exp)
Cas9G.exp_df$cell_barcode <- rownames(Cas9G.exp_df)

good_cells<- Cas9G.exp_df$cell_barcode

#make list of protospacers in same cell order
p.spacers_plot$Citrine <- NA
for(p in good_cells){x<- which(p.spacers_plot$cell_barcode == p)
y<- which(Cas9G.exp_df$cell_barcode == p)
p.spacers_plot$Citrine[x] <- Cas9G.exp_df$Cas9G.exp[y]}

p.spacers_plot

# Create a bar chart of total guide expression showing cut off (Figure S5C)
p<-ggplot(p.spacers_plot, aes(x = "", y = umis_log)) +
  geom_violin(fill = "blue") +
  geom_hline(yintercept = threshold_low, linetype = "dashed", color = "red")  +
  theme_classic()
p

# Create a violin plot without counts below the threshold
p2<-ggplot(p.spacers_plot, aes(x = "", y = umis_log_thresh)) +
  geom_violin(fill = "blue") +
  theme_classic()
p2

#generate a dataframe to add to metadata of screen_seq_no.blood
#figure out what guides were detected
gRNA_detected <- sort(unique(p.spacers_plot$feature))
#figure out cells detected with guides
cells_tot<- rownames(screen_seq@meta.data)
#cells_tot
gRNAs_df <- lapply(gRNA_detected, function(g){
  guide<- c(rep(0,length(cells_tot)))
  cells<- p.spacers_plot$cell_barcode[which(p.spacers_plot$feature == g & p.spacers_plot$umis_log_thresh != 'NA')]
  g_cells<- lapply(cells, function(c){
    c<- which(cells_tot == c)
    c})
  g_cells<- unlist(g_cells)
  g_cells
  guide[g_cells]<- 1
  guide})
gRNAs_df<- do.call(cbind, gRNAs_df)
colnames(gRNAs_df)<- gRNA_detected
rownames(gRNAs_df)<- cells_tot
gRNAs_df

#add sum of all guides per cell to dataframe
gRNAs_df_fin  <- cbind(gRNAs_df, rowSums(gRNAs_df))
gRNA_names_1 <- c(gRNA_detected, "SUM_all")
colnames(gRNAs_df_fin)<- gRNA_names_1
gRNAs_df_fin

#find sum of guides to same gene

names_all <- sapply(strsplit(colnames(gRNAs_df_fin),"_"), `[`, 1)
names_all <- paste("g_",names_all, sep = "")
#names_all 


names_uni <- names_all[!duplicated(names_all)]
#names_uni

length <- c(1:length(names_uni))
#length

gRNAs_tot_df <- lapply(length, function(i){n <- names_uni[i]
p <- grep(n, names_all)
f <- min(p)
e <- max(p)
if(f==e){v<-gRNAs_df_fin[,f]}
else{v <- rowSums(gRNAs_df_fin[,f:e])}
v})


gRNAs_tot_df_fin <- do.call(cbind,gRNAs_tot_df)
colnames(gRNAs_tot_df_fin)<- names_uni                    
gRNAs_tot_df_fin

#add to screen_seq_no.blood and include sample information for scramble guides
gRNAs_df <- cbind(gRNAs_df_fin,gRNAs_tot_df_fin)
screen_seq_no.blood<-AddMetaData(screen_seq_no.blood, metadata = gRNAs_df, col.name = colnames(gRNAs_df)) 


s1<- which(screen_seq_no.blood@meta.data$SCRAMBLE_g1 >0 & screen_seq_no.blood@meta.data$orig.ident == "sample 1")
s2<- which(screen_seq_no.blood@meta.data$SCRAMBLE_g1 >0 & screen_seq_no.blood@meta.data$orig.ident == "sample 2")
s3<- which(screen_seq_no.blood@meta.data$SCRAMBLE_g1 >0 & screen_seq_no.blood@meta.data$orig.ident == "sample 3")
s4<- which(screen_seq_no.blood@meta.data$SCRAMBLE_g1 >0 & screen_seq_no.blood@meta.data$orig.ident == "sample 4")

screen_seq_no.blood$SCRAMBLE_g1.1 <- c(rep(0, nrow(screen_seq_no.blood)))
screen_seq_no.blood$SCRAMBLE_g1.2 <- c(rep(0, nrow(screen_seq_no.blood)))
screen_seq_no.blood$SCRAMBLE_g1.3 <- c(rep(0, nrow(screen_seq_no.blood)))
screen_seq_no.blood$SCRAMBLE_g1.4 <- c(rep(0, nrow(screen_seq_no.blood)))

screen_seq_no.blood$SCRAMBLE_g1.1[s1] <- screen_seq_no.blood$SCRAMBLE_g1[s1]
screen_seq_no.blood$SCRAMBLE_g1.2[s2] <- screen_seq_no.blood$SCRAMBLE_g1[s2]
screen_seq_no.blood$SCRAMBLE_g1.3[s3] <- screen_seq_no.blood$SCRAMBLE_g1[s3]
screen_seq_no.blood$SCRAMBLE_g1.4[s4] <- screen_seq_no.blood$SCRAMBLE_g1[s4]

## plot number of guides per cell in dataset
##determine how many guides per cell in dataset

Sum_guides<- screen_seq_no.blood@meta.data$SUM_all
#Sum_guides
ncount_RNA<- screen_seq_no.blood@meta.data$nCount_RNA

cluster<- screen_seq_no.blood@meta.data$seurat_clusters
value_counts <- table(Sum_guides)
value_counts_df <- as.data.frame(value_counts)
#value_counts_df

#determine expression of fluorophore Citrine which marks CAS9 expression
Cas9G.exp <- GetAssayData(object = screen_seq_no.blood, assay = "RNA", slot = "data")["CITRINE",]
Cas9G.exp<- unlist(Cas9G.exp)

g2GFP <- data.frame(guide_sum = Sum_guides,
                    n_countRNA = ncount_RNA,
                    Citrine = Cas9G.exp)

# Create a bar chart of total guide expression (Figure 3E)
ggplot(value_counts_df, aes(x = Sum_guides, y = Freq)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_cartesian(ylim = c(0, 1500))+
  theme_classic()


# Create a scatter plot (Figure S5B)
ggplot(g2GFP, aes(x = guide_sum, y = n_countRNA)) +
  geom_point( fill = "blue") +
  theme_classic()

## generate a umap transform of the screen dataset onto the umap of the wildtype dataset to determine cells that sit within the tissues of interest (ie. primitive streak, CLE, and neural tube)

chick_1013S_dW <- ScaleData(chick_1013S_dW, verbose = FALSE)
chick_1013S_dW <- RunPCA(chick_1013S_dW, npcs = 30, verbose = FALSE, return.model=TRUE)
chick_1013S_dW <- RunUMAP(chick_1013S_dW, reduction = "pca", dims = 1:30, verbose = FALSE, return.model=TRUE)

chick.query <- screen_seq_no.blood
chick.anchors <- FindTransferAnchors(reference = chick_1013S_dW, query = screen_seq_no.blood,
                                     dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = chick.anchors, refdata = chick_1013S_dW$seurat_clusters,
                            dims = 1:30)
chick.query <- AddMetaData(chick.query, metadata = predictions)

chick_1013S_dW.umap <- RunUMAP(chick_1013S_dW, dims = 1:30, reduction = "pca", return.model = TRUE)
chick.query <- MapQuery(anchorset = chick.anchors, reference = chick_1013S_dW.umap, query = chick.query,
                        refdata = list(celltype = "seurat_clusters"), reference.reduction = "pca", reduction.model = "umap")

#use wildtype cluster references in screen_seq_no.blood dataset
screen_seq_no.blood@meta.data$seurat_cluster <- chick.query$predicted.id

## Plot resultant transform of screen onto wildtype dataset UMAP (Figure 3B)

p1 <- DimPlot(chick_1013S_dW, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 6,
              repel = TRUE, pt.size = 2, )  + ggtitle("Reference annotations")
p2 <- DimPlot(chick.query, reduction = "ref.umap", group.by = "predicted.id",label = TRUE,
              cols = c( '#F8766D','#E7851E','#00BCD6','#00B3F2','#29A3FF','#9C8DFF','#D277FF','#F166E8','#FF61C7','#FF689E','#D09400','#B2A100','#89AC00','#45B500','#00BC51','#00C087','#00C0B2'),
              label.size = 6, repel = TRUE,pt.size = 1)  + ggtitle("Query transferred labels")
p1 + p2
p2

##Similar to wildtype dataset, subset and recluster just the primitive streak, CLE, and neural tube
all.genes <- rownames(chick.query)
Idents(chick.query) <- chick.query$predicted.id

PRENEURAL <- subset(chick.query, idents = c(1,4,7,8,15))

##determine how many guides per cell in dataset and plot against detected Citrine

Sum_guides<- PRENEURAL@meta.data$SUM_all
#Sum_guides
ncount_RNA<- PRENEURAL@meta.data$nCount_RNA
nfeature_RNA<- PRENEURAL@meta.data$nFeature_RNA

cluster<- PRENEURAL@meta.data$seurat_clusters
value_counts <- table(Sum_guides)
value_counts_df <- as.data.frame(value_counts)

value_counts_df

#determine expression of fluorophore Citrine which marks CAS9 expression
Cas9G.exp <- GetAssayData(object = PRENEURAL, assay = "RNA", slot = "data")["CITRINE",]
Cas9G.exp<- unlist(Cas9G.exp)

g2GFP <- data.frame(guide_sum = Sum_guides,
                    nfeature_RNA = nfeature_RNA,
                    nCount_RNA = ncount_RNA,
                    Citrine = Cas9G.exp)

# Create a bar chart of total guide expression
p1<- ggplot(value_counts_df, aes(x = Sum_guides, y = Freq)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Value Counts", x = "Value", y = "Count")+
  #coord_cartesian(ylim = c(0, 1000))+
  theme_classic()

p1

#Remove Citrine from dataset so that it doesn't effect clustering
counts <- GetAssayData(PRENEURAL, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c("CITRINE"))),]
PRENEURAL <- subset(PRENEURAL, features = rownames(counts))

## now we recluster follow the Seurat clustering pipeline from the Satija lab
PRENEURAL <- NormalizeData(PRENEURAL, normalization.method = "LogNormalize", scale.factor = 10000)
PRENEURAL <- FindVariableFeatures(PRENEURAL,selection.method = "vst", nfeatures = 2000)

#because we are dealing with cycling progenitors we regressed out hte cell cycle

##regress out cell cycle
exp.mat <- read.table(file = "nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
PRENEURAL <- CellCycleScoring(PRENEURAL, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
PRENEURAL <- CellCycleScoring(PRENEURAL, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Visualize the distribution of cell cycle markers across
RidgePlot(PRENEURAL, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
PRENEURAL <- RunPCA(PRENEURAL, features = c(s.genes, g2m.genes))
DimPlot(PRENEURAL, reduction = "umap")
PRENEURAL <- ScaleData(PRENEURAL, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(PRENEURAL))

PRENEURAL <- RunPCA(PRENEURAL, features =VariableFeatures(object = PRENEURAL))
PRENEURALURALNEURALNEURAL <- RunPCA(PRENEURAL, features = c(s.genes, g2m.genes))
DimPlot(PRENEURAL)

VizDimLoadings(PRENEURAL, dims = 1:2, reduction = "pca")
DimPlot(PRENEURAL, pt.size = 1.5, reduction = "pca")
DimHeatmap(PRENEURAL, dims = 1:3, cells = 500, balanced = TRUE)

PRENEURAL <- JackStraw(PRENEURAL, num.replicate = 100)
PRENEURAL <- ScoreJackStraw(PRENEURAL, dims = 1:20)
JackStrawPlot(PRENEURAL, dims = 1:20)
ElbowPlot(PRENEURAL, ndims = 50, reduction = "pca" )

## again use an elbow plot of resoultion inertias to determine effective resolution to cluster the dataset
PRENEURAL <- FindNeighbors(PRENEURAL, dims = 1:20)
res <- c(0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0)
for(x in res){PRENEURAL <- FindClusters(PRENEURAL, resolution = x)}

res_0.2 <- PRENEURAL@meta.data['RNA_snn_res.0.2']
res_0.4 <- PRENEURAL@meta.data['RNA_snn_res.0.4']
res_0.6 <- PRENEURAL@meta.data['RNA_snn_res.0.6']
res_0.8 <- PRENEURAL@meta.data['RNA_snn_res.0.8']
res_1.0 <- PRENEURAL@meta.data['RNA_snn_res.1']
res_1.2 <- PRENEURAL@meta.data['RNA_snn_res.1.2']
res_1.4 <- PRENEURAL@meta.data['RNA_snn_res.1.4']
res_1.6 <- PRENEURAL@meta.data['RNA_snn_res.1.6']
res_1.8 <- PRENEURAL@meta.data['RNA_snn_res.1.8']
res_2.0 <- PRENEURAL@meta.data['RNA_snn_res.2']
res_2.2 <- PRENEURAL@meta.data['RNA_snn_res.2.2']
res_2.4 <- PRENEURAL@meta.data['RNA_snn_res.2.4']
res_2.6 <- PRENEURAL@meta.data['RNA_snn_res.2.6']
res_2.8 <- PRENEURAL@meta.data['RNA_snn_res.2.8']
res_3.0 <- PRENEURAL@meta.data['RNA_snn_res.3']

mutli_res <- data.frame(c(res_0.2,res_0.4,res_0.6,res_0.8,res_1.0,res_1.2,res_1.4,res_1.6,res_1.8,res_2.0,res_2.2,res_2.4,res_2.6,res_2.8, res_3.0))
rownames(mutli_res)<- colnames(PRENEURAL)
mutli_res

#distances with clusters squared and then summed. then take average across all cluster sums and this is score,
#want a low score which means point are close to centriod
#then plot inertia score against resolution

min <- which( colnames(PRENEURAL@meta.data)=="RNA_snn_res.0.2" )
min
max <- which( colnames(PRENEURAL@meta.data)=="RNA_snn_res.3" )
max
no <- which( colnames(PRENEURAL@meta.data)=="seurat_clusters" )
no
res <- seq(min, max, 1)
res <- res[res != no ]
res

euclidean <- function(a,b)sqrt(rowSums((a-t(b))^2))

Df_res <- lapply(res, function(r){
  clusters <- c(0:last(levels(PRENEURAL@meta.data[,r])))
  df <- PRENEURAL@meta.data[r]
  
  inertias <- lapply(clusters, function(c){
    
    b <- GetAssayData(PRENEURAL[which(df == c)])
    tot <- ncol(b)
    a <- colSums(b)/tot
    inert <- euclidean(a,b)
    #print (inert)
    inertias <- sum(inert)/length(inert)
    inertias
    
  })
  
  inertias_final <- unlist(inertias)
  inertias_fin <- data.frame(avg_inertia=mean(inertias_final))
  inertias_fin$resolution <-  r
  inertias_fin
  
})

Res_table_inertias <- do.call(rbind,Df_res)
Res_table_inertias

Res_table_inertias[1,2] <- 0.2
Res_table_inertias[2,2] <- 0.4
Res_table_inertias[3,2] <- 0.6
Res_table_inertias[4,2] <- 0.8
Res_table_inertias[5,2] <- 1.0
Res_table_inertias[6,2] <- 1.2
Res_table_inertias[7,2] <- 1.4
Res_table_inertias[8,2] <- 1.6
Res_table_inertias[9,2] <- 1.8
Res_table_inertias[10,2] <- 2.0
Res_table_inertias[11,2] <- 2.2
Res_table_inertias[12,2] <- 2.4
Res_table_inertias[13,2] <- 2.6
Res_table_inertias[14,2] <- 2.8
Res_table_inertias[15,2] <- 3.0

ggplot(Res_table_inertias, aes(x=resolution, y=avg_inertia)) + geom_point(size = 5, col = "blue")

#select resolution to continue analysis with
badres <- c('RNA_snn_res.0.2','RNA_snn_res.0.4','RNA_snn_res.0.6',
            'RNA_snn_res.0.8','RNA_snn_res.1.2','RNA_snn_res.1',
            'RNA_snn_res.1.4','RNA_snn_res.1.6','RNA_snn_res.1.8',
            'RNA_snn_res.2','RNA_snn_res.2.2','RNA_snn_res.2.4',
            'RNA_snn_res.2.6','RNA_snn_res.2.8','RNA_snn_res.3')
good <- 'RNA_snn_res.1'
bres <- badres[badres != good ]
bres

PRENEURAL@meta.data <- PRENEURAL@meta.data[,!(names(PRENEURAL@meta.data) %in% bres)]
names(PRENEURAL@meta.data)
PRENEURAL@meta.data$seurat_clusters<- PRENEURAL@meta.data$RNA_snn_res.1

#run umap
PRENEURAL <- RunUMAP(PRENEURAL, dims = 1:20)

#cluster 13 has markers of blood (LMO2) and surface ectoderm (GATA2) contamination so remove from dataset and rerun clustering
PRENEURAL <- subset(PRENEURAL, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,14))


Idents(PRENEURAL)<- PRENEURAL@meta.data$seurat_clusters

k<- DimPlot(PRENEURAL, reduction = 'umap', pt.size = 2,label = FALSE
            , label.size = 12 )
k

##plot genes used in Figure 3D

my_levels <- c("5", "6","10","2","8","0","7","3","1","4","12","9","11","13")
Idents(PRENEURAL)<-factor(Idents(PRENEURAL), levels = my_levels)
p<- DotPlot(PRENEURAL,  features = c("FOXC2","MSGN1","PARAXIS","TBX6", "TBXT","FGF8","ID4","MSX1","CDX2","BMP4", "NOTO","PAX7","SOX2","NKX1-2",
                                     "SOX21","SOX1",  "PAX6","NKX6-2","OLIG2",
                                     "SOX9","SOX10","FOXD3","WNT6","SHH","FOXA2",  
                                     "NEUROG1","ONECUT1"),
            cols = c("grey","blue"),dot.scale = 17)+
  theme(axis.text.x = element_text(angle = 90))

p

## plot overlaping genes from Figure S5A

blend <- TRUE
x<- "PAX7"
y<- "NKX6-2"
p <- FeaturePlot(PRENEURAL, reduction = "umap", features = c(x,y),  blend = TRUE, blend.threshold = 0.1,
                 cols = if (blend) {c("lightgrey", "blue", "red")} else {c("lightgrey", "blue")},  
                 pt.size = 1, combine = FALSE)
p <- lapply(p, function(x) x + NoLegend())
p<- CombinePlots(p)
p

## plot NMPs marked by overlap of SOX2/3 expression and TBXT expression
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
p<-DimPlot(PRENEURAL, reduction = "umap", dims=c(1,2),
           group.by = 'expression_type', cols = c("purple","grey","red", "red","red","deepskyblue", "purple","purple" ), pt.size = 6) 
p


## to generate the distribution of scramble guides across clusters in Figure 3F

clusters <- PRENEURAL@meta.data[,'seurat_clusters']

scramble_dist <- cbind(PRENEURAL@meta.data[,c("SCRAMBLE_g1","SCRAMBLE_g2")], clusters)

p<- which(scramble_dist[1] == 1)
scramble_dist[p,1]<- "SCRAMBLE"
p<- which(scramble_dist[2] == 1)
scramble_dist[p,2]<- "SCRAMBLE"
#scramble_dist


clusters <- c(0:13)
dist1<- lapply(clusters, function(x){i<- length(which(scramble_dist[,'SCRAMBLE_g1'] == 'SCRAMBLE' & 
                                                        scramble_dist[,'clusters'] == x))
tot <- length(which(scramble_dist[,'clusters'] == x))
value<- i/tot
value})

dist2<- lapply(clusters, function(x){i<- length(which(scramble_dist[,'SCRAMBLE_g2'] == 'SCRAMBLE' & 
                                                        scramble_dist[,'clusters'] == x))
tot <- length(which(scramble_dist[,'clusters'] == x))
value<- i/tot
value})

dist1<-unlist(dist1)
dist2<- unlist(dist2)
all<- c(dist1, dist2)
clusters <- c(clusters, clusters)
all
clusters

guide <- c(rep("scramble g1",14), rep('scramble g2',14) )
guide
final <- cbind(all, clusters, guide)
final <- as.data.frame(final)
final$clusters <- as.character(final$clusters)
final$all <- as.numeric(final$all)
final

p<-ggplot(final, aes(x=factor(clusters, levels = c(0:13)), y = all, fill = guide))+ 
  geom_bar(position="dodge", stat="identity")+
  coord_cartesian( ylim = c(0, 1))+
  theme_classic()
p

## To determine what guides had effective knockouts to generate Figure 3G we first calculated the gene expression of each gene for each guide
guides <- gRNA_detected

#exclude scramble and get just name of gene
x<-length(guides)
sa<- (1:x)
genes<- lapply(sa, function(s){
  guides <- str_split(guides,"_")
  g<- guides[[s]][[1]] 
  g})
genes <- unlist(genes)
genes<- unique(genes)
scr<- which(genes=="SCRAMBLE")
genes<- genes[-scr]


Gene.exp<- lapply(genes, function(gene){
  if(gene == "F2RL1"){guides <- 4} else{if(gene == "PAX7" |gene == "MAPK3"){guides <- c(1:3)} 
    else{ guides<- c(1:4)}}
  
  gene.exps<- lapply(guides, function(g){
    guide<- paste0(gene,"_g",g)
    Idents(PRENEURAL)<- PRENEURAL@meta.data[,which(colnames(PRENEURAL@meta.data) == guide)]
    #if low detection count as not significant or unable to make claim
    if(sum(PRENEURAL@meta.data[,which(colnames(PRENEURAL@meta.data) == guide)]) <3){
      gene.exp <- data.frame(p_val = 1,
                             avg_log2FC = 0,
                             pct.1 = 0,
                             pct.2 = 0,
                             p_val_adj = 1)
      gene.exp$gene<- gene
      gene.exp$guide<- guide
      gene.exp$g_num<- as.character(g)
      gene.exp$cluster <- c
    } else {
      guide_sum<- paste0("g_",gene)
      gene.exp <- FindMarkers(object = PRENEURAL, features = gene,
                              ident.1 = 1)
      
      if(nrow(gene.exp) == 0){
        gene.exp <- data.frame(p_val = 1,
                               avg_log2FC = 0,
                               pct.1 = 0,
                               pct.2 = 0,
                               p_val_adj = 1)
        gene.exp$gene<- gene
        gene.exp$guide<- guide
        gene.exp$g_num<- as.character(g)
        gene.exp$cluster <- c
      }else{
        gene.exp$gene<- gene
        gene.exp$guide<- guide
        gene.exp$g_num<- as.character(g)
        gene.exp$cluster <- c}}
    gene.exp})
  gene.exps<- do.call(rbind, gene.exps)
  gene.exps})
Gene.exp<- do.call(rbind,Gene.exp)
Gene.exp

# do the same thing for scramble guides in each sample
Gene.exp.scram<- lapply(genes, function(gene){
  guides<- c("SCRAMBLE_g1","SCRAMBLE_g2")
  gene.exps<- lapply(guides, function(g){
    #guide<- paste0(gene,"_g",g)
    gene_sum <- paste0("g_", gene)
    Idents(screen_seq)<- PRENEURAL@meta.data[,which(colnames(PRENEURAL@meta.data) == gene_sum)]
    exp<- subset(PRENEURAL, idents = "0")
    Idents(exp)<- PRENEURAL@meta.data[,which(colnames(PRENEURAL@meta.data) == g)]
    #if low detection count as not significant
    if(sum(exp@meta.data[,which(colnames(exp@meta.data) == g)]) <3){
      gene.exp <- data.frame(p_val = 1,
                             avg_log2FC = 0,
                             pct.1 = 0,
                             pct.2 = 0,
                             p_val_adj = 1)
      gene.exp$gene<- gene
      gene.exp$guide<- g
      gene.exp$g_num<- g
    } else {
      
      gene.exp <- FindMarkers(object = exp, features = gene,
                              ident.1 = 1)
      
      if(nrow(gene.exp) == 0){
        gene.exp <- data.frame(p_val = 1,
                               avg_log2FC = 0,
                               pct.1 = 0,
                               pct.2 = 0,
                               p_val_adj = 1)
        gene.exp$gene<- gene
        gene.exp$guide<- g
        gene.exp$g_num<- g  
      }else{
        gene.exp$gene<- gene
        gene.exp$guide<- g
        gene.exp$g_num<- g}}
    gene.exp})
  gene.exps<- do.call(rbind, gene.exps)
  gene.exps})
Gene.exp.scram<- do.call(rbind,Gene.exp.scram)
Gene.exp.scram

## find average expression of scrambles

gene<- unique(Gene.exp.scram$gene)

scram.avg.exp<- lapply(gene, function(g){
  df<- Gene.exp.scram[which(Gene.exp.scram$gene == g),]
  mean2FC<- mean(df$avg_log2FC)
  df2<- data.frame(gene = g,
                   avg_log2FC = mean2FC)
  df2
})
scram.avg.exp<- do.call(rbind, scram.avg.exp)
scram.avg.exp

#correct change in gene expression by subtracting changes also found in cells with scramble 
Gene.exp.tot.corrected<- lapply(genes, function(gene){
  guides<- Gene.exp$guide[which(Gene.exp$gene == gene)]
  guides
  scr.avg<- scram.avg.exp$avg_log2FC[which(scram.avg.exp$gene ==gene )]
  scr.avg
  gene.cor<- lapply(guides, function(g){
    guide.exp<- Gene.exp$avg_log2FC[which(Gene.exp$guide == g)]
    guide.exp<- guide.exp - scr.avg
    df<- data.frame(avg_log2FC.cor = guide.exp,
                    guide.cor= g)
    df})
  gene.cor<- do.call(rbind,gene.cor)
  gene.cor
})
Gene.exp.tot.corrected<- do.call(rbind, Gene.exp.tot.corrected)
Gene.exp.tot.corrected

#now finally add a column that quantifies whether the guide shows larger negative change in the query gene 
#compared to the scramble control

Gene.exp.tot<-cbind(Gene.exp,Gene.exp.tot.corrected)

#list of genes
genes<- Gene.exp.tot$gene
genes<-unique(genes)
genes
sd<- sd(Gene.exp.tot$avg_log2FC)
avg<- mean(Gene.exp.tot$avg_log2FC)
avg
n <- length(unique(genes))*2
error <- qnorm(0.975)*sd/sqrt(n)
left <- 0-error
left
right <- 0+error


index<- nrow(Gene.exp.tot)
index<- c(1:index)
less_than<- lapply(index, function(i){
  p<- Gene.exp.tot$avg_log2FC[i]
  if(p<(left)){l<- "down"} else {if(p> right){l<- "up"}else{l<-"noise"}}
  #t<- Gene.exp.tot$p_val_adj[i]
  #if(t<0.05){l<- "noise"}
  guide<- Gene.exp.tot$guide[i]
  p<- data.frame(successful = l)
  p
})
less_than<- do.call(rbind, less_than)
Gene.exp.tot<- cbind(Gene.exp.tot, less_than)
Gene.exp.tot

# plot
p<- ggplot(Gene.exp.tot, aes(x=avg_log2FC, y= gene)) + 
  geom_point(aes( color=successful, shape = g_num ), size = 15,position=position_dodge(0.05))+
  scale_shape_manual(values = c( 15, 16, 17, 18))+ #0,1,2,5))+
  scale_color_manual(values = c( "dodgerblue", "grey", "red"))+
  #geom_errorbar(aes(ymin=avg_log2FC-error, ymax=avg_log2FC+error), width=0.2, position= position_dodge(0.7)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  #geom_hline(yintercept = (avg-sd), linetype = "dashed", color = "red")+
  geom_vline(xintercept = (left), linetype = "dashed", color = "grey")+
  geom_vline(xintercept = (right), linetype = "dashed", color = "grey")


options(repr.plot.width=25, repr.plot.height=15)
p<- p+theme(
  panel.background = element_rect(fill = "transparent",
                                  colour = NA_character_), # necessary to avoid drawing panel outline
  panel.grid.major = element_blank(), # get rid of major grid
  panel.grid.minor = element_blank(), # get rid of minor grid
  plot.background = element_rect(fill = "transparent",
                                 colour = NA_character_), # necessary to avoid drawing plot outline
  legend.background = element_rect(fill = "transparent"),
  legend.box.background = element_rect(fill = "transparent"),
  legend.key = element_rect(fill = "transparent")
)

p

##for remaining analysis identify "good guides" and save in csv file
bad_g<- unique(Gene.exp.tot$guide[which(Gene.exp.tot$successful == "up" |
                                          Gene.exp.tot$successful == "noise" )])
bad_g
good_g<- unique(Gene.exp.tot$guide[which(Gene.exp.tot$successful == "down")])
good_g<- unique(c(good_g, "SCRAMBLE_g1","SCRAMBLE_g2"))
sort(good_g)

write.csv(good_g,file = "good_guides.csv")

## next we generated tests by randomly subsetting each sample to use in confidence testing
samples <- c("sample 1","sample 2", "sample 3","sample 4")
s1 <- which(PRENEURAL@meta.data$orig.ident == 'sample 1')
sample_size = floor(0.8*length(s1))
s1t1<- sample(s1, sample_size)
s1t2 <- sample(s1, sample_size)
s1t3<- sample(s1,sample_size)
s1t4<- sample(s1,sample_size)
s1t5<- sample(s1, sample_size)

s2 <- which(PRENEURAL@meta.data$orig.ident == 'sample 2')
sample_size = floor(0.8*length(s2))
s2t1<- sample(s2, sample_size)
s2t2 <- sample(s2, sample_size)
s2t3<- sample(s2, sample_size)
s2t4<- sample(s2, sample_size)
s2t5<- sample(s2, sample_size)

s3 <- which(PRENEURAL@meta.data$orig.ident == 'sample 3')
sample_size = floor(0.8*length(s3))
s3t1<- sample(s3, sample_size)
s3t2 <- sample(s3, sample_size)
s3t3<- sample(s3, sample_size)
s3t4<- sample(s3, sample_size)
s3t5<- sample(s3, sample_size)

s4 <- which(PRENEURAL@meta.data$orig.ident == 'sample 4')
sample_size = floor(0.8*length(s4))
s4t1<- sample(s4, sample_size)
s4t2 <- sample(s4, sample_size)
s4t3<- sample(s4, sample_size)
s4t4<- sample(s4, sample_size)
s4t5<- sample(s4, sample_size)

test1 <- rep(0,dim(PRENEURAL@meta.data)[1])
test2 <- rep(0,dim(PRENEURAL@meta.data)[1])
test3 <- rep(0,dim(PRENEURAL@meta.data)[1])
test4 <- rep(0,dim(PRENEURAL@meta.data)[1])
test5 <- rep(0,dim(PRENEURAL@meta.data)[1])

test1[c(s1t1,s2t1,s3t1,s4t1)] <- 't1'
test2[c(s1t2,s2t2,s3t2,s4t2)] <- 't2'
test3[c(s1t3,s2t3,s3t3,s4t3)] <- 't3'
test4[c(s1t4,s2t4,s3t4,s4t4)] <- 't4'
test5[c(s1t5,s2t5,s3t5,s4t5)] <- 't5'


PRENEURAL@meta.data['test1'] <- test1
PRENEURAL@meta.data['test2'] <- test2
PRENEURAL@meta.data['test3'] <- test3
PRENEURAL@meta.data['test4'] <- test4
PRENEURAL@meta.data['test5'] <- test5

###make dataframe of all theguides and what clusters they fall into, then take random sampling to get 5 test samples within each guide to test for outliers

guides4<- good_g
guides4



ALL <- lapply(guides4, function(gu){TBXT<- rownames(PRENEURAL@meta.data)
TBXT<- data.frame(TBXT)    
TBXT$count<- PRENEURAL@meta.data[,gu]
c <- which( colnames(PRENEURAL@meta.data)=="RNA_snn_res.1" )
test <- which(colnames(PRENEURAL@meta.data)=="tests" )
TBXT$seurat_clusters <- PRENEURAL@meta.data[,c]
TBXT$test1 <- PRENEURAL@meta.data[,'test1']
TBXT$test2 <- PRENEURAL@meta.data[,'test2']
TBXT$test3 <- PRENEURAL@meta.data[,'test3']
TBXT$test4 <- PRENEURAL@meta.data[,'test4']
TBXT$test5 <- PRENEURAL@meta.data[,'test5']
TBXT$orig.ident <- PRENEURAL@meta.data$orig.ident                                    
guide <- c(rep(gu,nrow(TBXT)))
TBXT$guide <- guide
g<- strsplit(gu,"_")
g<- unlist(g)
g<- g[[1]]
TBXT$gene <- g
colnames(TBXT)<- c("X","count","seurat_clusters","test1","test2","test3","test4","test5","orig.ident",
                   "guide","gene")

TBXT})

ALL<- do.call(rbind,ALL)
ALL

#do the same for the scramble guides
s1 <- which( colnames(PRENEURAL@meta.data)=="SCRAMBLE_g1" )
s2 <- which( colnames(PRENEURAL@meta.data)=="SCRAMBLE_g2" )
c <- which( colnames(PRENEURAL@meta.data)=="RNA_snn_res.1" )
t<- which( colnames(PRENEURAL@meta.data)=="tests" )

SCR <- PRENEURAL@meta.data[,c(s1:s2)]
SCR$avg <- rowSums(SCR)/2
SCR$seurat_clusters <- PRENEURAL@meta.data[,c]
SCR$test1 <- PRENEURAL@meta.data[,'test1']
SCR$test2 <- PRENEURAL@meta.data[,'test2']
SCR$test3 <- PRENEURAL@meta.data[,'test3']
SCR$test4 <- PRENEURAL@meta.data[,'test4']
SCR$test5 <- PRENEURAL@meta.data[,'test5']
SCR$orig.ident <- PRENEURAL@meta.data[,'orig.ident']

SCR

## next we calculate the p - proportion of cells with gene knockout per test set in certain cluster, q - proportion of cells with scramble per test set in certain cluster
## p_all - proportion of all cells with knockout in certain cluster, q - proportion of all cells with scramble in certain cluster

ALL_done <- ALL
tests <- c('t1','t2','t3','t4','t5')
max <- max(as.numeric(ALL_done$seurat_clusters), na.rm = TRUE)
max
sc <- c(0:(max-1))

df_seurat<- lapply(sc, function(x){
  guide <- which(ALL_done$seurat_clusters == x &
                   ALL_done$gene != "SCRAMBLE")
  
  g <- ALL_done[guide, -10]
  scr <- ALL_done[which(ALL_done$gene == "SCRAMBLE"),-11]
  colnames(scr)<-     c("X","count","seurat_clusters","test1","test2","test3","test4",
                        "test5","orig.ident","gene")
  g<- rbind(g,scr)
  g_uni<-unique(g)
  good_q<- unique(unlist(g_uni$gene))
  target_p<- lapply(good_q, function(q){
    samples<- unique(g$orig.ident[which(g_uni$gene == q)])
    samples<- unlist(samples)
    good_p<- lapply(samples, function(s){num<- which(g_uni$orig.ident == s&
                                                       g_uni$gene == q)
    num})
    good_p<- unlist(good_p)
    good_p<- g_uni[good_p,]
    test<- lapply(tests, function (t){
      n <- which(good_p$seurat_clusters == x & (good_p$test1 == t |
                                                  good_p$test2 == t |
                                                  good_p$test3 == t |
                                                  good_p$test4 == t |
                                                  good_p$test5 == t ) & 
                   good_p$gene == q &
                   good_p$count > 0)
      l <- good_p[n,]
      l_fin<- nrow(l)
      b<- PRENEURAL@meta.data[which(PRENEURAL@meta.data$RNA_snn_res.1 == x & 
                                      (PRENEURAL@meta.data$test1 == t |
                                         PRENEURAL@meta.data$test2 == t |
                                         PRENEURAL@meta.data$test3 == t |
                                         PRENEURAL@meta.data$test4 == t |
                                         PRENEURAL@meta.data$test5 == t )),]
      indices <- which(b$orig.ident %in% samples)
      b<- b[indices, ] 
      len <- nrow(b)
      if(len == 0){len <- 0.0000001}
      p <- l_fin/len
      indices <- which(SCR$orig.ident %in% samples)
      SCR1<- SCR[indices, ]
      n1 <- which(SCR1$seurat_clusters == x & (SCR1$test1 == t |
                                                 SCR1$test2 == t |
                                                 SCR1$test3 == t |
                                                 SCR1$test4 == t |
                                                 SCR1$test5 == t ) & 
                    (SCR1["SCRAMBLE_g1"] != 0 ))
      n2 <- which(SCR1$seurat_clusters == x & (SCR1$test1 == t |
                                                 SCR1$test2 == t |
                                                 SCR1$test3 == t |
                                                 SCR1$test4 == t |
                                                 SCR1$test5 == t ) & 
                    (SCR1["SCRAMBLE_g2"] != 0))
      l1 <- (length(n1))
      l2 <- (length(n2))
      l3<- (l1+l2)/2
      if(l3 == 0){fin <- 0.0000001}
      q <- l3/len
      fin<- data.frame(p,q)
      fin})
    test<-do.call(rbind,test)
    n_all <- which(good_p$seurat_clusters == x &
                     good_p$gene == q &
                     good_p$count > 0)
    l_all <- good_p[n_all,]
    lall_fin<- nrow(l_all)
    if(lall_fin == 0){lall_fin <- 0.0000001}
    b_all<- PRENEURAL@meta.data[which(PRENEURAL@meta.data$RNA_snn_res.1 == x),]
    indices <- which(b_all$orig.ident %in% samples)
    b_all<- b_all[indices, ]
    len_a <- nrow(b_all)
    if(len_a == 0){len_a <- 0.0000001}
    test$p_all<- lall_fin/len_a
    indices <- which(SCR$orig.ident %in% samples)
    SCR1<- SCR[indices, ]
    n_all1 <- which(SCR1$seurat_clusters == x & SCR1["SCRAMBLE_g1"] != 0 )
    n_all2 <- which(SCR1$seurat_clusters == x &  SCR1["SCRAMBLE_g2"] != 0)
    l_all_fin <- (length(n_all1)+length(n_all2))/2
    fin_a <- l_all_fin/len_a
    test$q_all <- fin_a
    test$tests <- tests
    test$guide <- q
    test$seurat_cluster <-x
    test$n <- len_a
    test})
  target_p<- do.call(rbind, target_p)
  target_p})
df_seurat_all<- do.call(rbind, df_seurat)
df_seurat_all_good_sum<-df_seurat_all
df_seurat_all_good_sum

## to determine which gene knockouts were significantly enriched or depleted in certain clusters across the dataset, we performed both chi-squared tests and Kullback-Leibler divergence tests

# for Chi-squared test reorient data to grid of genes vers clusters where p_all is the value
genes<- unique(df_seurat_all_sum$guide)
#o<- which(genes=="OLIG2")
#genes<- genes[-o]
clusters<- unique(df_seurat_all_sum$seurat_cluster)
clusters

chi<- lapply(genes, function(x){
  g.row<- lapply(clusters, function(c){
    g<- df_seurat_all_sum$p_all[which(df_seurat_all_sum$seurat_cluster == c & 
                                        df_seurat_all_sum$guide == x &
                                        df_seurat_all_sum$tests == "t1")]
    g
  })
  g.row<- do.call(cbind, g.row)
  colnames(g.row)<- clusters
  g.row
})
chi<- do.call(rbind, chi)
rownames(chi)<- genes
chi<- chi*100
chi

# Perform a chi-square test
test_result <- chisq.test(chi)
test_result

chisq.test(chi)$statistic

# Extract the residuals
residuals <- test_result$residuals

# Plot the heatmap
heatmap(residuals, 
        Rowv = NA, 
        Colv = NA,
        col = cm.colors(256),  # Color palette
        scale = "none",        # Don't scale the values
        margins = c(5, 10))    # Add extra margins

# to see if certain gene knockouts are different from scramble create contingency tables of the query gene distribution across clusters, the average of all gene knockout distributions, and the scramble distributions
genes<- rownames(chi)
genes
length(genes)
genes<- genes[-c(24)]
genes



contingency_table <- lapply(genes, function(g){
  avg_clus <- colSums(chi)/nrow(chi)
  table <- data.frame(g = chi[g,],
                      SCRAMBLE_g1 = chi[c(25),],
                      SCRAMBLE_g2 = chi[c(24),],
                      average = avg_clus
  )
  p.val<- chisq.test(table)$p.value
  X.sq<- chisq.test(table)$statistic
  table.chi<- data.frame(X.sq,p.val)
  table.chi
})
contingency_table<- do.call(rbind, contingency_table)
rownames(contingency_table) <- genes
colnames(contingency_table) <- c("chi.squared", "p.value")
contingency_table$guide <- genes
contingency_table

#plot results see in Figure 4A
plot2<- ggplot(contingency_table, aes(x = -log2(chi.squared), y = -log2(p.value), fill = guide, label = guide)) +
  geom_point(aes(colour = guide, label = guide), size = 15) +
  scale_shape_manual(values= c(rep(c(15, 16,17,18),10)))+
  geom_hline(yintercept = -log2(0.001), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log2(0.0001), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log2(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(box.padding = 0.5, max.overlaps = 15, size = 12) +
  theme_classic()
plot2

# remove unsignificant guides 

bad<-which(contingency_table$p.value>0.05)
#bad <- bad(-c(which(bad == 'SCRAMBLE_g1'| bad =='SCRAMBLE_g2')))
ns.guides<- contingency_table$guide[bad]
#ns.guides<- c(ns.guides, "SCRAMBLE_g1")
scrams <- which(ns.guides == 'SCRAMBLE_g1'| ns.guides =='SCRAMBLE_g2')
ns.guides<- ns.guides[-scrams]
##remove unsignificant gudies

chi.bad<- lapply(ns.guides, function(g){
  bad<- which(rownames(chi) == g)
  bad
})
chi.bad<- unlist(chi.bad)


chi.sig<- chi[-(chi.bad),]
chi.sig

## cluster using k.means clustering to see what gene knockouts behave similarly in terms of enrichment scores (Figure 4B)
kmeans_result_all <- kmeans(chi.sig, centers = 6,nstart = 10)
kmeans_result <- kmeans(chi.sig, centers = 5,nstart = 10)

p<-fviz_cluster(kmeans_result, chi.sig, axes = c(1,2),pointsize = 5, line = 10, shape = 16, labelsize = 20,repel = TRUE)+
  geom_point(shape = 16, size = 6, color = "black")+
  theme_bw()
p

#plot PCA plots from Figure 4D

res.pca <- prcomp(chi.sig,  scale = TRUE)
p<- fviz_pca_var(res.pca, axes = c(1,2), repel = TRUE)
p

# to make heatmap of chisquared enrichments in Figure 4C

##make fancy heatmap

my_hclust_gene <- hclust(dist(chi.sig), method = "complete")

# install if necessary
#install.packages("dendextend")

# load package
library(dendextend)

as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)

my_gene_col <- data.frame(cutree(tree = as.dendrogram(my_hclust_gene), k = 5))
colnames(my_gene_col)<- "k.means"
my_gene_col
p<- rownames(residuals)
rownames(my_gene_col)<- p
sub_samp_ordered <- order(my_gene_col$k.means)
my_sample_col <- data.frame(cluster = c(0:13))
row.names(my_sample_col) <- colnames(residuals)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=14)
color_list

min<- min(residuals)
max<- max(residuals)
range(residuals)
values1 <- seq(-2, min, by = -0.1)
values2 <- seq( 2, max, by = 0.1)# Adjust the range and step as needed
colors <- c(seq(min, -2.01, by = 0.1),seq(-2,2, by = 0.1),seq(2.01,max, by = 0.1))

my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(values1)-3),
                rep( "lightgrey", length(seq(-2,2, by = 0.1))),
                c(colorRampPalette(colors = c("lightpink","darkred"))(n = length(values2)-3)))
my_colour = list(
  kmeans = c('1' = "#D8F3DC", '2' = "#95D5B2", "3" = "#52B788", '4' = "#2D6A4F", '5' ='#1B4332', '6' = "#081C15"),
  cluster = c('0' =color_list[1], '1' =color_list[2],'2' =color_list[3],'3' =color_list[4],'4' =color_list[5],'5' = color_list[6],
              '6' =color_list[7], '7' =color_list[8],'8' =color_list[9],'9' =color_list[10],'10' =color_list[11],
              '11' =color_list[12], '12' =color_list[13],'13' =color_list[14])
  
)

my_sample_col$cluster <- factor(my_sample_col$cluster, levels = c(0:13))


P<- pheatmap(n.sig[-c(23,24),],  annotation_col = my_sample_col, annotation_row = my_gene_col, 
             annotation_colors = my_colour, 
             border_color=NA,,
             color = my_palette , 
             fontsize_col = 6,
             fontsize_row = 15,
             cluster_cols = TRUE,
             cluster_rows= TRUE,
             gaps_row=c(6,9,12,15))
P

## to perform Kullback Leiber divergence tests used Philentropy R package
library (philentropy)

max <- max(as.numeric(ALL_done$seurat_clusters), na.rm = TRUE)
sc <- c(0:(max))
guides<- unique(df_seurat_all_sum$guide)
KLD <- lapply(guides, function(g){
  tests<- c("t1","t2","t3","t4","t5")
  test <- lapply(tests, function(t){
    klds <- lapply(sc, function(c){
      p <- df_seurat_all_sum$p[which(df_seurat_all_sum$test == t &
                                       df_seurat_all_sum$seurat_cluster == c &
                                       df_seurat_all_sum$guide == g)]
      if(length(p)==0){p <- 0}
      q <- df_seurat_all_sum$q[which(df_seurat_all_sum$test == t & 
                                       df_seurat_all_sum$seurat_cluster == c &
                                       df_seurat_all_sum$guide == g)]
      if(length(q)==0){q <- 0}
      n <- df_seurat_all_sum$n[which(df_seurat_all_sum$test == t & 
                                       df_seurat_all_sum$seurat_cluster == c&
                                       df_seurat_all_sum$guide == g )]
      if(length(n)==0){n <- 0}
      #p_all<- df_seurat_all_sum$p_all[which(df_seurat_all_sum$test == t &
      # df_seurat_all_sum$seurat_cluster == c & 
      # df_seurat_all_sum$guide == g)]
      #q_all <- df_seurat_all_sum$q_all[which(df_seurat_all_sum$test == t & 
      #df_seurat_all_sum$seurat_cluster == c &
      # df_seurat_all_sum$guide == g)]
      #if(q_all==0){ratio_all <- 1} else {ratio_all <- (p_all/q_all)}
      cluster<- c
      #m <- data.frame( p, p_all, q, q_all, n, ratio_all, cluster,t,g)
      m <- data.frame( p,  q,  n,  cluster,t,g)
      m})
    klds
    klds<- do.call(rbind,klds)
    #KLs<- klds[,c(1,3)]
    KLs<- klds[,c(1,2)]
    sum_p<- colSums(KLs)
    if(sum_p[1] == 0){p1 <- KLs[,1]} else {p1<- KLs[,1]/sum_p[1]}
    if(sum_p[2] == 0){KLs_fin<- 0} else {
      p1<- KLs[,1]/sum_p[1]
      q1<- KLs[,2]/sum_p[2]
      KLs2<- rbind(p1,q1)
      KLs_fin<- KL(KLs2, test.na = TRUE, unit = "log2", est.prob = NULL, epsilon = 1e-7)
      KLs_fin}
    klds$KLD_all <- KLs_fin
    klds})
  test<- do.call(rbind, test)
  test})
KLD_all<-do.call(rbind, KLD)
KLD_all

##perform anova between guide KLDs
max <- max(as.numeric(ALL_done$seurat_clusters), na.rm = TRUE)
max
sc <- c(0:(max))
two.way <- aov(KLD_all ~ g, data = KLD2[which(KLD2$cluster == 0),])
two.way
tukey.two.way<-TukeyHSD(two.way)
tukey.two.way
result_table <- tukey.two.way$g %>% as_tibble(rownames = "comparitors")
result_table
colnames(result_table)<- c("comparitors", "diff","lwr", "upr", "p_adj")
rt_new2 <- result_table[grepl("SCRAMBLE_g2", result_table$comparitors), ] 
rt_new1 <- result_table[grepl("SCRAMBLE_g1", result_table$comparitors), ] 
rt_new<- rbind(rt_new1,rt_new2)
rt_new

##clean up table

two_way<-rt_new[-(which(rt_new$comparitors=="SCRAMBLE_g1-SCRAMBLE_g2"|rt_new$comparitors=="SCRAMBLE_g2-SCRAMBLE_g1")),] 
p<- str_split(two_way$comparitors,"-") 
p <- unlist(p)
p<- p[which(grepl("SCRAMBLE_g1", p) == FALSE&
              grepl("SCRAMBLE_g2", p) == FALSE)]
p
two_way$gene <- p
scr<- rt_new[(which(rt_new$comparitors=="SCRAMBLE_g1-SCRAMBLE_g2"|rt_new$comparitors=="SCRAMBLE_g2-SCRAMBLE_g1")),] 
scr$gene <- c("SCRAMBLE 1", "SCRAMBLE 2")

two_way<- rbind(two_way,scr)
two_way

#calculate p.value comparing scramble 1 or scramble 2
KLD2<-KLD_all
KLD2$P.val1 <- c(rep(1,nrow(KLD2)))
KLD2$P.val2 <- c(rep(1,nrow(KLD2)))

genes<- two_way$gene
genes
x<-'TBXT'
for(x in genes){
  
  p.val<- two_way$p_adj[which(two_way$gene == x)]
  p.val
  KLD2$P.val1[which(KLD2$g == x)]<- p.val[1]
  KLD2$P.val2[which(KLD2$g == x)]<- p.val[2]}

#plot seen in figure S6A
plot2<- ggplot(kld_pval, aes(x = kld, y = -log2(p.val1), fill = guide)) +
  geom_point(aes(shape = guide,colour = guide), size = 4) +
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "grey") +
  scale_shape_manual(values= c(rep(c(15, 16,17,18),10)))+
  geom_hline(yintercept = -log2(0.001), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log2(0.0001), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log2(0.05), linetype = "dashed", color = "red") +
  theme_classic()
plot2

## save seurat objects
saveRDS(PRENEURAL, file = "PRENEURAL.rds")
saveRDS(screen_seq_no.blood, file = "screen_seq_no.blood.rds")