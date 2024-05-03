##Install packages
install.packages("Seurat")
install.packages("dplyr")
install.packages("rjson")
install.packages("viridis")
install.packages("pheatmap")
install.packages("rstatix")
install.packages("dendextend")
install.packages("clustree")
install.packages("ggdendro")
install.packages("Hmisc")
install.packages("ggpubr")
install.packages("AICcmodavg")
install.packages("enrichR")
install.packages("ggvenn")


library(Matrix)
library(Seurat)
library(dplyr)
library(patchwork)
library(scales)
library (ggplot2)
library(sctransform)
library(rjson)
library(clustree)
library(viridis)
library(RColorBrewer)
library(data.table)
library(pheatmap)
library(ggdendro)
library(Hmisc)
library(ggpubr)
library(broom)
library(AICcmodavg)
library(rstatix)
library(dendextend)
library(cowplot)
library(stringr)
library(ggrepel)
library(httr)
library(igraph)
library(enrichR)
library(igraph)
library(tidygraph)

## This analysis pipeline uses the scMAGeCK R pipeline to generate differential expression of pooled CRISPR screens with single cell RNA sequencing data
##EXAMPLES ON THIS PAGE https://bitbucket.org/weililab/scmageck/src/master/tools/R/virtual_facs_functions.R

### first load BARCODE file contains cell identity information, generated from the cell identity collection step in this case from the Cell Ranger analysis

GUIDES1 <- read.csv(file = 'sample1_protospacer_calls_per_cell.csv')
x <- paste(GUIDES1$cell_barcode, "_1", sep='')
GUIDES1$cell_barcode <- x
GUIDES2 <- read.csv(file = 'sample2_protospacer_calls_per_cell.csv')
x <- paste(GUIDES2$cell_barcode, "_2", sep='')
GUIDES2$cell_barcode <- x
GUIDES3 <- read.csv(file = 'sample3_protospacer_calls_per_cell.csv')
x <- paste(GUIDES3$cell_barcode, "_3", sep='')
GUIDES3$cell_barcode <- x
GUIDES4 <- read.csv(file = 'sample4_protospacer_calls_per_cell.csv')
x <- paste(GUIDES4$cell_barcode, "_4", sep='')
GUIDES4$cell_barcode <- x
GUIDES <- rbind(GUIDES1,GUIDES2,GUIDES3,GUIDES4)

SGRNAS1 <- read.csv(file ='sample1_feature_reference.csv')
SGRNAS2 <- read.csv(file ='sample2_feature_reference.csv')
SGRNAS3 <- read.csv(file ='sample3_feature_reference.csv')
SGRNAS4 <- read.csv(file ='sample4_feature_reference.csv')
SGRNAS <- rbind(SGRNAS1, SGRNAS2,SGRNAS3,SGRNAS4)

##generate a dataframe of the guide, gene target, umi_count, and cell

barcode2 <- SGRNAS1[1]
#barcode2
sgrnas <- SGRNAS1[5]
df2 <- cbind(barcode2,sgrnas)

length <- c(1:nrow(GUIDES))
#length

df2<- lapply(length, function(x){
  barcode <- unlist(strsplit(as.character(GUIDES[x,3]),'|',fixed = T))
  barcode
  gene <- sapply(strsplit(barcode, "_"), `[`, 1)
  #gene
  umi_count <- unlist(strsplit(as.character(GUIDES[x,4]),'|',fixed = T))
  #umi_count
  cell <- rep(list(GUIDES[x,1]), length(barcode))
  #cell
  cell <- unlist(cell)
  sgrna<- lapply( barcode, function(b){sgrna<- df2$sequence[which(df2$id == b)]
  sgrna})
  sgrna<- unlist(sgrna)
  
  df <- cbind(cell,barcode,sgrna,gene,umi_count)
  df})

gRNAs_df_fin  <- do.call(rbind,df2)
gRNAs_df_fin

gRNAs_df_fin<-as.data.frame(gRNAs_df_fin)
gRNAs_df_fin$umi_count<-as.numeric(gRNAs_df_fin$umi_count)

## filter out bad calls by removing guide calls below 2 standard deviations of the log(umi count)

log_p.counts<- log(as.numeric(gRNAs_df_fin$umi_count))
gRNAs_df_fin$log_umi_count<- log_p.counts
range(log_p.counts)
mean<- mean(log_p.counts)
sd<- sd(log_p.counts)
sd
threshold_low <- mean - (2*sd)
threshold_low

gRNAs_df_fin<- gRNAs_df_fin[-which(gRNAs_df_fin$log_umi_count< threshold_low),]
gRNAs_df_fin

## load significant guides (good_g from Screen.analysis.R) and filter dataframe to only include good guides
barcode.sig<- read.csv(file = "good_g.csv")
gRNAs_df_fin<-data.frame(gRNAs_df_fin)
gRNAs_df_fin_sig<- lapply(barcode.sig, function(g){ p<- which(gRNAs_df_fin$barcode == g)
cells<- gRNAs_df_fin[p,]
cells})
gRNAs_df_fin_sig<-do.call(rbind,gRNAs_df_fin_sig)
gRNAs_df_fin_sig

#save dataframe as table for future reference
write.table(gRNAs_df_fin_sig, file = 'scMAGeCK_barcode_sig.txt', append = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

file<-'scMAGeCK_barcode_sig.txt'
BARCODE<- read.delim2(file, header = TRUE, sep = "\t")
BARCODE<- gRNAs_df_fin_sig
BARCODE<- BARCODE[,(-6)]

### RDS can be a Seurat object or local RDS file path that contains the scRNA-seq dataset
#RDS <- our PRENEURAL dataset
RDS <- readRDS(file = "PRENEURAL.rds")
RDS<-JoinLayers(object = RDS, layer = 'scale.data')
RDS
Idents(RDS)<- RDS@meta.data$seurat_clusters
DimPlot(RDS, reduction = "umap", label = TRUE, label.size = 10)

## load Entropy Ranked genes and subset RDS so that only entropy sorted genes are included
entropy.ranked.genes<- read.csv(file = "ES_Cluster_Ranked_Genes_Lists.csv" )
ES_500<- entropy.ranked.genes[c(1:500),-1]

ES_500<- unique(c(ES_500$X1,ES_500$X10,ES_500$X10,ES_500$X14,
                  ES_500$X4,ES_500$X7,ES_500$X8,ES_500$Mesoderm,ES_500$PS ))

RDS_ES.500 <- subset(RDS, features = ES_500)

## load scMAGeCK functions
scmageck_lr <- function(BARCODE, RDS, NEGCTRL, SELECT_GENE = NULL, LABEL = NULL, PERMUTATION = NULL,
                        SIGNATURE = NULL, SAVEPATH = "./", LAMBDA = 0.01, GENE_FRAC = 0.01) {
  if (!is.null(LABEL)) {
    data_label = LABEL
  } else {
    data_label = "sample1"
  }
  
  if (!is.null(PERMUTATION)) {
    n_permutation = as.integer(PERMUTATION)
  } else {
    n_permutation = 10000
  }
  
  if (!is.null(SIGNATURE)) {
    data_signature=SIGNATURE
    message(paste("run_signature: TRUE"))
  } else {
    data_signature = NULL
    message(paste("run_signature: FALSE"))
  }
  
  # read cell assignment and libray file ####
  bc_dox = read.table(BARCODE, header = TRUE, as.is = TRUE)
  
  if (sum(colnames(bc_dox) %in% c("cell", "barcode", "gene")) != 3) {
    stop("cell, barcode, or gene column names not found in barcode file.")
  }
  
  guide_count = table(bc_dox$cell)
  ncnt = table(table(bc_dox$cell))
  message(paste("Total barcode records:", nrow(bc_dox)))
  
  # load neg control guides ####
  ngctrlgenelist = strsplit(NEGCTRL, ",")[[1]]
  message(paste("Neg Ctrl guide:", paste(ngctrlgenelist, collapse = ";")))
  
  # read Seurat RDS file ####
  if (is.character(RDS)) {
    message(paste("Reading RDS file:", RDS))
    targetobj = readRDS(RDS)
  } else {
    targetobj = RDS
  }
  # check if names are consistent 
  nmatch = sum(bc_dox[, 1] %in% colnames(x = targetobj))
  if (nmatch == 0) {
    message("Cell names in expression matrix and barcode file do not match. Try to remove possible trailing \"-1\"s...")
    if (length(grep("-\\d$", bc_dox[, 1])) > 0) {
      bc_dox[, 1] = sub("-\\d$", "", bc_dox[, 1])
    }
    nmatch = sum(bc_dox[, 1] %in% colnames(x = targetobj))
    if (nmatch == 0) {
      stop("No cell names match in expression matrix and barcode file.")
    }
  }
  # bc_dox[,1]=sub('-\\d$','',bc_dox[,1])
  
  # convert to ind_matrix ####
  ind_matrix <- frame2indmatrix(bc_dox, targetobj)
  message(paste("Index matrix dimension:", nrow(ind_matrix), ",", ncol(ind_matrix)))
  
  # try to perform matrix regresson on single genes ####
  mat_for_single_reg = single_gene_matrix_regression(targetobj, selected_genes_list = SELECT_GENE, 
                                                     ngctrlgene = ngctrlgenelist, indmatrix = ind_matrix, high_gene_frac = GENE_FRAC)
  Xmat = mat_for_single_reg[[1]]
  
  # Xmat[,which(colnames(Xmat)%in%ngctrlgenelist)[1]]=1 # already integrated into function
  Ymat = mat_for_single_reg[[2]]
  
  # Optional function
  # Get the results based on gmt file
  if(!is.null(data_signature)){
    x <- scan(data_signature, what = "", sep = "\n")
    x <- strsplit(x, "\t") # split string by white space
    max.col <- max(sapply(x, length))
    cn <- paste("V", 1:max.col, sep = "")
    gmt <- read.table(data_signature, fill = TRUE, col.names = cn)
    # gmt <- read.delim(data_signature, header = FALSE)
    gmt <- t(as.matrix(gmt))
    colnames(gmt) <- gmt[1, ]
    gmt <- gmt[-1:-2, ]
    message(paste("Total signature records:", ncol(gmt)))
    sig_mat <- getsigmat(Ymat, gmt_file = gmt)
    if (ncol(sig_mat) > 0) {
      Amat_sig_lst = getsolvedmatrix_with_permutation_cell_label(Xmat, sig_mat, lambda = LAMBDA, npermutation = n_permutation)
      sig_score = Amat_sig_lst[[1]]
      sig_pval = Amat_sig_lst[[2]]
      sig_re <- getsigresult(signature_score = sig_score, signature_pval = sig_pval)
      sig_re$Fdr <- p.adjust(sig_re$p_value, method = "fdr")
      write.table(data.frame(sig_re), file = file.path(SAVEPATH, paste(data_label, "_signature.txt", sep = "")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      return(list(data.frame(sig_re)))
    }
  } else {
    # remove values in Y mat
    Amat_pm_lst = getsolvedmatrix_with_permutation_cell_label(Xmat, Ymat, lambda = LAMBDA, npermutation = n_permutation)
    Amat = Amat_pm_lst[[1]]
    Amat_pval = Amat_pm_lst[[2]]
    if(!is.null(SAVEPATH)){
      write.table(data.frame(Perturbedgene = rownames(Amat), Amat), file = file.path(SAVEPATH, paste(data_label,
                                                                                                     "_score.txt", sep = "")), sep = "\t", quote = FALSE, row.names = FALSE)
      write.table(data.frame(Perturbedgene = rownames(Amat), Amat_pval), file = file.path(SAVEPATH, paste(data_label,
                                                                                                          "_score_pval.txt", sep = "")), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    return(list(data.frame(Perturbedgene = rownames(Amat), Amat), data.frame(Perturbedgene = rownames(Amat),
                                                                             Amat_pval)))
  }
}
TRUE

# function definitions ##### version: 02-22-2019 should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R

frame2indmatrix <- function(bc_d, targetobj) {
  
  rnm = unique(bc_d$cell)
  cnm = unique(bc_d$gene)
  scalef = getscaledata(targetobj)
  message(paste(length(rnm), "..."))
  message(paste(ncol(scalef), "..."))
  # if(sum(rnm%in%colnames(scalef))==0){ message('Cell names in expression matrix and barcode file do
  # not match. Try to remove possible trailing '-1's...')
  # if(length(grep('-\\d$',colnames(scalef)))>0){ colnames(scalef) =
  # sub('-\\d$','',colnames(scalef)) } if(length(grep('-\\d$',rnm))>0){ rnm =
  # sub('-\\d$','',rnm) } }
  rnm = rnm[!is.na(rnm)]
  rnm = rnm[rnm %in% colnames(scalef)]
  if (length(rnm) == 0) {
    stop("Cell names do not match in expression matrix and barcode.")
  }
  cnm = cnm[!is.na(cnm)]
  ind_matrix = matrix(rep(FALSE, length(rnm) * length(cnm)), nrow = length(rnm))
  rownames(ind_matrix) = rnm
  colnames(ind_matrix) = cnm
  row <- bc_d[, 'cell']
  col <- bc_d[, 'gene']
  test <- (row %in% rnm) & (col %in% cnm)
  idx <- cbind(row[test], col[test])
  ind_matrix[idx]  <- TRUE
  return(ind_matrix)
}
TRUE

# function definitions ##### version: 02-22-2019 should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R


getscaledata <- function(targetobj, scaled = TRUE) {
  # if scaled=FALSE, return raw.data
  if ("scale.data" %in% names(attributes(targetobj))) {
    if (scaled) {
      scalef = targetobj@scale.data  # for version 2
    } else {
      scalef = targetobj@raw.data  # for version 2
    }
  } else {
    if (scaled) {
      scalef = GetAssayData(object = targetobj, slot = "scale.data")
    } else {
      scalef = GetAssayData(object = targetobj, slot = "counts")
    }
  }
  return(scalef)
}
TRUE
# Get the gene signature expression matrix from Ymat
getsigmat <- function(Ymat, gmt_file) {
  colgmt <- colnames(gmt_file)
  sig_mat <- data.frame(row.names = rownames(Ymat))
  for (num in (1:ncol(gmt_file))) {
    genes <- as.character(gmt_file[, num])
    if (any(genes %in% colnames(Ymat))) {
      genes <- genes[genes %in% colnames(Ymat)]
      Ymat_sig <- as.data.frame(Ymat[, genes])
      Ymat_sig$m <- rowMeans(Ymat_sig)
      sig_mat <- cbind(sig_mat, as.data.frame(Ymat_sig$m)) # identify whether the genome is mouse or human
    } else {
      genes <- capitalize(tolower(genes))
      if(any(genes %in% colnames(Ymat))) {
        genes <- genes[genes %in% colnames(Ymat)]
        Ymat_sig <- as.data.frame(Ymat[, genes])
        Ymat_sig$m <- rowMeans(Ymat_sig)
        sig_mat <- cbind(sig_mat, as.data.frame(Ymat_sig$m))
      } else {
        message(paste(colnames(gmt_file)[num], "can not found in this dataset"))
        colgmt <- subset(colgmt, colgmt != colnames(gmt_file)[num])
        next
      }
    }
  }
  if (ncol(sig_mat) > 0) {
    colnames(sig_mat) <- colgmt
    sig_mat <- as.matrix(sig_mat)
  } else {
    message("No signatures can be found in this dataset")
  }
  return(sig_mat)
}
TRUE
#convert signature score/pvalue matrix to the format of final results
getsigresult <- function(signature_score, signature_pval) {
  p <- rownames(signature_score)
  output <- data.frame(rep(p, ncol(signature_score)))
  sig_name <- NULL
  for (i in colnames(signature_score)) {
    sig_name <- rbind(sig_name, data.frame(rep(i, nrow(signature_score))))
  }
  sig_sc <- NULL
  for (n in 1:ncol(signature_score)) {
    sig_sc <- rbind(sig_sc, data.frame(signature_score[, n]))
  }
  sig_p <- NULL
  for (n in 1:ncol(signature_pval)) {
    sig_p <- rbind(sig_p, data.frame(signature_pval[, n]))
  }
  output <- cbind(output, sig_name)
  output <- cbind(output, sig_sc)
  output <- cbind(output, sig_p)
  colnames(output) <- paste(c("sgrna", "gene_signature", "LR_score", "p_value"))
  rownames(output) <- (1:nrow(output))
  return(output)
}
TRUE

# function definitions ##### version: 02-22-2019 should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R

# construct a matrix of Y=XA, Y= (cells*expressed genes), X=(cells* KO genes), A=(KO genes *
# expressed genes)

single_gene_matrix_regression <- function(targetobj, ngctrlgene = c("NonTargetingControlGuideForHuman"), 
                                          indmatrix = NULL, high_gene_frac = 0.01, selected_genes_list = NULL) {
  # return X matrix and Y matrix for regression note that all the ngctrlgene are merged into one
  # column, 'NegCtrl' if indmatrix is provided, the Xmat will be constructed from indmatrix
  outlier_threshold = 0.95
  rawf = getscaledata(targetobj, scaled = FALSE)
  scalef = getscaledata(targetobj)
  if(nrow(rawf)>0){
    select_genes = rownames(rawf)[which(rowSums(as.matrix(rawf) != 0) >= ncol(rawf) * high_gene_frac)]
    message(paste('Filter genes whose expression is greater than 0 in raw read count in less than',high_gene_frac,'single-cell populations.' ))
  } else {
    message(paste('Cannot find raw read count in Seurat object. Use scaled data instead, and filter genes whose expression is greater than 0 in less than',high_gene_frac,'single-cell populations.' ))
    select_genes = rownames(scalef)[which(rowSums(as.matrix(scalef) >= 0) >= ncol(scalef) * high_gene_frac)]
  }
  if (is.null(selected_genes_list) == FALSE) {
    select_genes = select_genes[select_genes %in% selected_genes_list]
  }
  select_genes = select_genes [!is.na(select_genes)& select_genes %in% rownames(scalef)]
  
  if (length(select_genes) == 0) {
    stop("No genes left for regression. Check your selected gene list.")
  }
  message(paste("Selected genes:", length(select_genes)))
  # browser()
  
  
  if (is.null(indmatrix)) {
    select_cells = rownames(targetobj@meta.data)[which(!is.na(targetobj@meta.data$geneID))]
  } else {
    select_cells = rownames(indmatrix)
    select_cells = select_cells[select_cells %in% colnames(scalef)]
  }
  select_cells = select_cells[!is.na(select_cells) & select_cells %in% colnames(scalef)]
  message(paste("Selected cells:", length(select_cells)))
  YmatT = scalef[select_genes, select_cells]
  
  Ymat = as.matrix(t(YmatT))  # (cells * expressed genes)
  if (is.null(indmatrix)) {
    tgf = targetobj@meta.data[select_cells, "geneID"]
    tgf[tgf %in% ngctrlgene] = "NegCtrl"
    tgphenotype = as.factor(tgf)
    Xmat = matrix(rep(0, length(select_cells) * length(unique(tgphenotype))), nrow = length(select_cells))
    rownames(Xmat) = select_cells
    colnames(Xmat) = levels(tgphenotype)
    Xmat[as.matrix(cbind(1:nrow(Xmat), as.numeric(tgphenotype)))] = 1
    Xmat[, "NegCtrl"] = 1  # set up base line
  } else {
    tgf = colnames(indmatrix)
    tgf[tgf %in% ngctrlgene] = "NegCtrl"
    tgphenotype = as.factor(tgf)
    
    Xmat = matrix(rep(0, length(select_cells) * length(unique(tgphenotype))), nrow = length(select_cells))
    rownames(Xmat) = select_cells
    colnames(Xmat) = levels(tgphenotype)
    for (cnl in colnames(indmatrix)) {
      cellns = which(indmatrix[, cnl] == TRUE)  #make sure indmatrix 
      if (cnl %in% ngctrlgene) {
        Xmat[cellns, "NegCtrl"] = 1
      } else {
        Xmat[cellns, cnl] = 1
      }
    }
    Xmat[, "NegCtrl"] = 1
    
  }  # end if
  
  # remove outliers
  Ymat_outlier = apply(Ymat, 2, function(X) {
    return(quantile(X, probs = outlier_threshold))
  })
  outlier_mat = t(matrix(rep(Ymat_outlier, nrow(Ymat)), ncol = nrow(Ymat)))
  Ymat_corrected = ifelse(Ymat > outlier_mat, outlier_mat, Ymat)
  Ymat = Ymat_corrected
  
  return(list(Xmat, Ymat))
}
TRUE

# function definitions ##### version: 02-22-2019 should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R


getsolvedmatrix_with_permutation_cell_label <- function(Xm, Ym, lambda = 0.01, npermutation = 1000) {
  Amat_ret = getsolvedmatrix(Xm, Ym, lambda = lambda)
  Amat_ret_higher = matrix(rep(0, ncol(Amat_ret) * nrow(Amat_ret)), nrow = nrow(Amat_ret))
  rownames(Amat_ret_higher) = rownames(Amat_ret)
  colnames(Amat_ret_higher) = colnames(Amat_ret)
  # permute N times randomly shuffle cell labels
  for (npm in 1:npermutation) {
    if (npm%%100 == 0) {
      message(paste("Permutation:", npm, "/", npermutation, "..."))
    }
    cells_shu = sample(rownames(Ym), nrow(Ym))
    Xm_s = Xm[cells_shu, ]
    Ym_s = Ym  # [cells_shu,]
    rownames(Ym_s) = cells_shu
    Amat_random = getsolvedmatrix(Xm_s, Ym_s, lambda = lambda)
    
    Amat_ret_higher = Amat_ret_higher + (abs(Amat_random) > abs(Amat_ret)) * 1
    # browser()
  }
  Amat_ret_higher = Amat_ret_higher/npermutation
  return(list(Amat_ret, Amat_ret_higher))
}
TRUE
# function definitions ##### version: 02-22-2019 should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R

# perform matrix decomposition get A for Y=XA A= (X^TX)^-1X^TY Y: (cells * expressed genes) X:
# design matrix, (cells * KO genes), can also be (cells * double KO genes) A: (KO genes * expressed
# genes)
getsolvedmatrix <- function(Xm, Ym, lambda = 0.01) {
  # Amat=solve(Xmat,Ymat) # solve AX=B, or Xmat * A =Ymat
  TMmat_g = (t(Xm) %*% Xm) + lambda * diag(ncol(Xm))
  
  Amat_g = solve(TMmat_g) %*% t(Xm) %*% Ym
  return(Amat_g)
}
TRUE

## lr_scMAGeCK
lr_result <- scmageck_lr(BARCODE=file, RDS=RDS_ES.500, LABEL='Perturbseq_3.0_scmageck_lr', 
                         NEGCTRL = 'SCRAMBLE', PERMUTATION = 1000)
lr_score <- lr_result[[1]]
lr_score
lr_score_pval <- lr_result[[2]]
lr_score_pval

#calculate the significantly changed genes based on the lr_score_pval and write as csv files
guides <- lr_score_pval[,1]
guides
lapply(guides, function(g){
  t.genes <- lr_score[g,]
  max<- ncol(t.genes)
  t.genes<- t.genes[,2:max]
  t.genes <- t(t.genes)
  len<- nrow(t.genes)
  t.pval<- lr_score_pval[g,]
  t.pval <- t.pval[,2:max]
  t.pval <- t(t.pval)
  t.genes<- data.frame(t.genes,t.pval)
  colnames(t.genes)<- c("lr_score", "p_val")
  t.genes$lr_score<- as.numeric(t.genes$lr_score)
  t.genes$p_val<- as.numeric(t.genes$p_val)
  t.genes$low.p_val <- NA
  low<- which(t.genes$p_val < 0.001)
  t.genes$low.p_val[low]<- rownames(t.genes[low,])
  NO <- c(rep("NO", len))
  t.genes$diffexpressed <- NO
  #t.genes
  ZERO<-which( t.genes$p_val == 0)
  t.genes$p_val[c(ZERO)]<- 0.00001
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
  UP<-which(t.genes$lr_score >= 0.1 & t.genes$p_val < 0.05)
  #UP
  DOWN<-which( t.genes$lr_score <= -0.1 & t.genes$p_val < 0.05)
  DOWN
  t.genes$diffexpressed[c(UP)] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  t.genes$diffexpressed[c(DOWN)] <- "DOWN"
  b <- which(t.genes$diffexpressed != "NO")
  t.genes$t.geneslabel <- NA
  t.genes$t.geneslabel[which(t.genes$diffexpressed != "NO")] <- rownames(t.genes[b,])
  t.genes <- t.genes[order( t.genes$lr_score,t.genes$p_val),]
  t.genes
  
  t.genes.down <- t.genes[which(t.genes$diffexpressed == 'DOWN'),]
  t.genes.up <- t.genes[which(t.genes$diffexpressed == 'UP'),]
  
  
  pUP<- paste0("Neural_scMAGeCK_",g, "_UP_sig.csv")
  pDOWN<- paste0("Neural_scMAGeCK_",g, "_DOWN_sig.csv")
  pALL<-paste0("Neural_scMAGeCK_",g, "_ALL_sig.csv")
  
  write.csv(t.genes.up, file =  pUP)
  write.csv(t.genes.down, file = pDOWN)
  write.csv(t.genes, file = pALL)})

##Calculate jaccard similarities between genes that have been misregualted (Figure 4E)
guides <- row.names(lr_score)
guides<-guides[-c(which(guides =="NegCtrl"))]

jaccards<- lapply(guides, function(g){
  pALL1<-paste0("Neural_scMAGeCK_",g, "_ALL_sig.csv")
  list1 <-  read.csv(file = pALL1)
  list1<- list1[which(list1$diffexpressed != "NO"),]
  list1<- list1$X
  if(length(list1)== 0){ guide = 0} else {
    guide<- lapply(guides, function(x){
      pALL2<-paste0("Neural_scMAGeCK_",x, "_ALL_sig.csv")
      pALL2
      list2 <-  read.csv(file = pALL2)
      list2<- list2[which(list2$diffexpressed != "NO"),]
      list2<- list2$X
      if(length(list2)== 0){ p<-data.frame(0,x)
      colnames(p)<-c(g, "gene2")} else {
        #calculate_jaccard_similarity 
        # Convert lists to sets
        set1 <- unique(list1)
        set2 <- unique(list2)
        
        # Calculate intersection and union
        intersection <- intersect(set1, set2)
        union <- union(set1, set2)
        
        # Handle case to avoid division by zero
        if (length(union) == 0) {
          return(0)
        }
        
        # Calculate Jaccard similarity
        jaccard_similarity <- length(intersection) / length(union)
        p<-data.frame(jaccard_similarity,x)
        colnames(p)<-c(g, "gene2")}
      p})
    guide<- do.call(rbind,guide)
    rownames(guide)<- guide$gene2
    guide<- guide[,1]}
  
  guide})

jaccards<- do.call(cbind, jaccards)
rownames(jaccards)<- guides
colnames(jaccards)<- guides
jaccards

p<- pheatmap(jaccards, 
             scale = "none", 
             border_color=NA,
             fontsize_row = 8,
             cluster_cols = TRUE, #as.hclust(row_dend),
             cluster_rows = TRUE) #as.hclust(row_dend))
p


## find kmeans clusters of enrichment distributions (Figure S6C)

#Assuming 'data' is your combined dataset
k<-c(1:10)
#k <- 3  # Number of clusters

set.seed(123)  # for reproducibility
type
k.means<- lapply(k, function(x){
  kmeans_result <- kmeans(jaccards, centers = x,nstart = 10)
  # Access cluster assignments
  clusters <- kmeans_result$cluster
  
  # Access cluster centers
  centroids <- kmeans_result$centers
  
  # Access within-cluster sum of squares
  wcss <- kmeans_result$tot.withinss
  k.mean<- data.frame(x,wcss)
  k.mean
  
})
k.means<-do.call(rbind,k.means)
k.means

p<- ggplot(k.means, aes(x=x, y=wcss)) + geom_point(size = 6, col = "blue")
p

##Calculate Gene Ontology and KEGG PAthway terms via Enrichr and generate jaccard distances

library(enrichR)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}else{print("F")}

guides <- row.names(lr_score)


#guides<- c(PC1.enrich, "NegCtrl")
#guides<- guides[-which(guides =="SCRAMBLE")]


websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015","KEGG_2021_Human")

#Calculate KEGG terms for a list
KEGG.terms<- function(g){
  if (websiteLive) {
    pALL1<-paste0("Neural_scMAGeCK_",g, "_ALL_sig.csv")
    list1 <-  read.csv(file = pALL1)
    list1<- list1[which(list1$diffexpressed != "NO"),]
    if(nrow(list1) == 0){
      file<-paste0("Neural_scMAGeCK_",g, "_ALL_sig.csv")
      p2<- "null"
      write.csv(p2, file =  file)}else {
        list1<- list1$X
        enriched <- enrichr(list1, dbs)
        
        p<-enriched[["KEGG_2021_Human"]]
        
        p2 <- p[which(p$P.value < 0.05),]
        if(length(p2) == 0){p3<- "null"}else {p3<- p2$Term}
        file<-paste0("Neural_scMAGeCK_",g, "_KEGGpaths_ALL_sig.csv")
        write.csv(p2, file =  file)
        p3
      } }}

# Calculate GO terms for a list
GO.terms<- function(g){
  if (websiteLive) {
    pALL1 <- paste0("Neural_scMAGeCK_",g, "_ALL_sig.csv")
    list1 <-  read.csv(file = pALL1)
    list1<- list1[which(list1$diffexpressed != "NO"),]
    if(nrow(list1) == 0){
      file<-paste0("Neural_scMAGeCK_",g, "_ALL_sig.csv")
      p2<- "null"
      write.csv(p2, file =  file)}else {
        list1<- list1$X
        enriched <- enrichr(list1, dbs)
        
        p<-enriched[["GO_Biological_Process_2015"]]
        
        p2 <- p[which(p$P.value < 0.05),]
        if(length(p2) == 0){p3<- "null"}else {p3<- p2$Term}
        file<-paste0("Neural_scMAGeCK_",g, "_GOterms_ALL_sig.csv")
        write.csv(p2, file =  file)
        p3
      } }}

lapply(guides, function(g){
  KEGG.terms(g)
  GO.terms(g)})


##generate jaccard distances used in Figure 4F and Figure S6B
test<- "KEGGpaths"
#test<- "GOterms"


jaccards<- lapply(guides, function(g){
  pALL1<-paste0("Neural_scMAGeCK_",g, "_GOterms_ALL_sig.csv")
  list1 <-  read.csv(file = pALL1)
  list1
  list1 <- list1[which(list1$P.value < 0.05),]
  list1<- list1$Term
  list1
  if(length(list1)== 0){ guide = rep(0, length(guides))} else {
    guide<- lapply(guides, function(x){
      pALL2<-paste0("Neural_scMAGeCK_",x, "_GOterms_ALL_sig.csv")
      pALL2
      list2 <-  read.csv(file = pALL2)
      list2
      list2 <- list2[which(list2$P.value < 0.05),]
      list2<- list2$Term
      #list2
      
      #calculate_jaccard_similarity 
      # Convert lists to sets
      set1 <- unique(list1)
      set2 <- unique(list2)
      
      # Calculate intersection and union
      intersection <- intersect(set1, set2)
      union <- union(set1, set2)
      
      # Handle case to avoid division by zero
      if (length(union) == 0) {
        return(0)
      }
      
      # Calculate Jaccard similarity
      jaccard_similarity <- length(intersection) / length(union)
      p<-data.frame(jaccard_similarity,x)
      colnames(p)<-c(g, "gene2")
      p})
    guide<- do.call(rbind,guide)
    rownames(guide)<- guide$gene2
    guide<- guide[,1]}
  
  guide})

jaccards<- do.call(cbind, jaccards)
rownames(jaccards)<- guides
colnames(jaccards)<- guides
jaccards

neg<- which(rownames(jaccards) == "NegCtrl")
jaccards2<-jaccards[-(neg), -neg]
range(jaccards2)

heat1<- pheatmap(jaccards2, 
                 scale = "none", 
                 border_color=NA,
                 fontsize_row = 8,
                 cluster_cols = TRUE,
                 cluster_rows = TRUE)
heat1


## load guides of interest (CMTM8, MLLT3, SPRY2, GREB1) to generate vendiagram in Figure 5A
a<- "CMTM8"
b<- "MLLT3"
c<- "GREB1"
d<- "SPRY2"

test<- "GOterms"

file_ALL<- paste0("Neural_scMAGeCK_",a, "_",test,"_ALL_sig.csv")
list1 <-  read.csv(file = file_ALL)
list1 <- list1[which(list1$P.value < 0.05),]
list1$KO<- a

file_ALL2<- paste0("Neural_scMAGeCK_",b, "_",test,"_ALL_sig.csv")
list2 <-  read.csv(file = file_ALL2)
list2 <- list2[which(list2$P.value < 0.05),]
list2$KO<- b

file_ALL3<- paste0("Neural_scMAGeCK_",c, "_",test,"_ALL_sig.csv")
list3 <-  read.csv(file = file_ALL3)
list3 <- list3[which(list3$P.value < 0.05),]
list3$KO<- c

file_ALL4<- paste0("Neural_scMAGeCK_",d, "_",test,"_ALL_sig.csv")
list4 <-  read.csv(file = file_ALL4)
list4 <- list4[which(list4$P.value < 0.05),]
list4$KO<- d

ALL.df<- rbind(list1, list2, list3, list4)

#now generate list for venndiagram and plot

x.ALL<- list ( CMTM8 = ALL.df$Term[which(ALL.df$KO == " CMTM8")], 
              MLLT3 = ALL.df$Term[which(ALL.df$KO == "MLLT3")],
              SPRY2 = ALL.df$Term[which(ALL.df$KO == "SPRY2")],
              GREB1 = ALL.df$Term[which(ALL.df$KO == "GREB1")])

library(ggvenn)
# Default plot
gvenn<- ggvenn(x.ALL,
               fill_color = c("red", "pink",  "orange", "yellow", "magenta","rose"),)

gvenn


## to generate Figure 5B we mannually assigned "parent terms" to the GO terms in ALL.df in a new column
## plot GO.terms

plot<- ggplot(data=ALL.df, aes(x=Odds.Ratio, y=-log10(Adjusted.P.value), color = parent.term, shape = KO)) + 
  geom_point(size = 6) + 
  geom_hline(yintercept=-log10(0.05), col="grey")+
  scale_shape_manual(values = c( 0,1,2,5))+
  scale_color_manual(values=c("black", "blue","green","grey", "dodgerblue","orange","pink", "purple","red","yellow")) + 
  ggtitle(g)
plot


#generate differentially expressed SEC members in MLLT3 KO Figure 5E

a<-"MLLT3"
file_ALL<- paste0("Neural_scMAGeCK_", a, "_ALL_sig_no_cell_cycle.csv")
a.genes <-  read.csv(file = file_ALL)


SEC<- c("AFF1","AFF2","AFF3","AFF4","TAT","NELFA","NELFB","NELFCD","NELFE","POLA1","MLLT1",
        "MLLT3","ELL","ELL2","ELL3","CDK9","ELL",
        "SUPT5H","SUPT4H1",
        "BRD4","JMJD6", "EAR", "PAF1")

GOI<- c(SEC)

low_lr<- which(a.genes$lr_score > -0.01 & a.genes$lr_score < 0.01)
a.genes_high <- a.genes[-c(low_lr),]
#t.genes_high

a.genes_high_ordered <- a.genes_high[order(a.genes_high$p_val), ]

SEC.genes<- lapply(GOI, function(g){
  p<- a.genes_high_ordered[which(a.genes_high_ordered$X == g),]
  p})
SEC.genes<- do.call(rbind,SEC.genes)
SEC.genes

plot<- ggplot(data=SEC.genes, aes(x=lr_score, y=-log10(p_val), col=diffexpressed)) + 
  geom_point(size = 6) + #xlim(-3, 3) + ylim(-0.01, 6.5)+
  geom_vline(xintercept=c(-0.1, 0.1), col="grey") +
  geom_hline(yintercept=-log10(0.05), col="grey")+
  geom_label_repel(data = subset(SEC.genes,  
                                 X =="RORA" |X=="AFF2"|X=="JMJD6"|
                                   X=="WNT5A"|X=="WNT5B"|X=="HES1"|X=="DLL1"),
                   aes(label=X, segment.color="red"), 
                   segment.color = 'grey50', max.overlaps = 100,  colour= "black") +
  #geom_label_repel(label=a.genes$X,force = 1, max.overlaps = 50, colour = "black")+
  scale_color_manual(values=c("blue", "darkgrey", "red")) + ggtitle(a)
plot

## To generate the heatmap of overlapping genes in Figure S7A first create a dataframe of the differentially expressed 
#genes of each gene of interest

file_ALL<- paste0("Neural_scMAGeCK_",a,"_ALL_sig.csv")
gene1 <-  read.csv(file = file_ALL)
gene1 <- gene1[which(gene1$p_val < 0.05),]
gene1$KO<- a

file_ALL2<- paste0("Neural_scMAGeCK_",b,"_ALL_sig.csv")
gene2 <-  read.csv(file = file_ALL2)
gene2 <- gene2[which(gene2$p_val < 0.05),]
gene2$KO<- b

file_ALL3<- paste0("Neural_scMAGeCK_",c,"_ALL_sig.csv")
gene3 <-  read.csv(file = file_ALL3)
gene3 <- gene3[which(gene3$p_val < 0.05),]
gene3$KO<- c

file_ALL4<- paste0("Neural_scMAGeCK_",d,"_ALL_sig.csv")
gene4 <-  read.csv(file = file_ALL4)
gene4 <- gene4[which(gene4$p_val < 0.05),]
gene4$KO<- d

gene.df<- rbind(gene1, gene2, gene3, gene4)

#find genes of overlap between KO
m2g<- intersect(gene1,gene2)
m2s<- intersect(gene1,gene3)
m2c<-intersect(gene1,gene4)
g2c<-intersect(gene2,gene4)
g2s<-intersect(gene2,gene3)
s2c<-intersect(gene3,gene4)
genes.sig<- unique(c(m2g,m2s,m2c,g2c,g2s,s2c))
genes.sig

# generate dataframe of the lr_scores of significant genes

genes.of.interest<- c("MLLT3","GREB1","CMTM8","SPRY2")

tot.genes<-genes.sig
x<- "FGFR3"
g<- "MLLT3"

lr_scores.goi<-lapply(genes.of.interest,function(g){
  KO<- lapply(tot.genes, function(x){
    p<-genes$p_val[which(genes$X == x& genes$KO == g)]
    if(p<0.05){
      q<- genes$lr_score[which(genes$X == x& genes$KO == g)]} else {q<- 0}
    df<- data.frame(lr_score = q)
    df})
  KO<- do.call(rbind, KO)
  KO
})
lr_scores.goi<- do.call(cbind, lr_scores.goi)
colnames(lr_scores.goi)<- genes.of.interest
rownames(lr_scores.goi)<- tot.genes

lr_scores.goi

#plot

plot<- pheatmap(lr_scores.goi,
                color = viridis(n = 256, alpha = 1, 
                                begin = 0, end = 1, option = "viridis"))
plot
