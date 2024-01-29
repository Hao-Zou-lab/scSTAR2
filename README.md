![scSTAR2图片](https://github.com/Hao-Zou-lab/scSTAR2/assets/115384996/301ac18e-ef77-4c39-a549-5fede0a008c2)
# scSTAR2

1. scSTAR2,R version 4.1.2.
2. scSTAR2 can map prognostic information to the spatial transcriptome, scSTAN-seq, ascATAC-seq data

# Step 1: intallation
1. install 'Seurat','pls', 'MASS','ggplot2','readxl','tcltk','wrMisc' R package from Cran.
2. installing STutility from github:
   dowaload the file 'STutility-1.1.1tar.gz'/ 'STutility-1.1.1.zip'and install the package from local path.
3. Installing scSTAR from github:
   1. We can install scSCTAR2 by downloading the installation file (scSTAR2_0.1.1.tar.gz) from the Assets page of Releases and installing it:
   ```
   install.packages("scSTAR2_0.1.1.tar.gz", repos = NULL, type="source")
   ```
   2. It should be noted that if you use this method, you must ensure that there is the file "scSTAR2_0.1.1.0.tar.gz" under setwd().
 4. Another method is to download the file "scSTAR2_0.1.1.0.tar.gz" and install the package from a local path.
 5.After making sure that your R package is downloaded (you can check whether there is scSTAR2 package in Packags), open the run_scSTAR2_demo_ST script and run it directly.
# Step 2: run scSTAR2
```
rm(list=ls())
graphics.off()
library(Seurat)#Seurat version 4.1.1
library(pls)#pls version 2.8.1
library(MASS)#MASS version 7.3-58.1
library(ggplot2)#ggplot2 version 3.3.6
library(readxl)#readxl version 1.4.1
library(tcltk)#tcltk version 4.1.2
library(STutility)#STutility version 1.1.1
library(wrMisc)#wrMisc version 1.13.0
library(scSTAR2)
setwd("G:/algorithms/scSTAR2/R package")
```
## load bulk RNA-seq data with prognosis information
```
extdir<-system.file("extdata",package="scSTAR2")
BRCA_bulk_1<-paste(extdir,"/BRCA_bulk.txt",sep="")
data_bulk = read.table(BRCA_bulk_1, sep = "\t", stringsAsFactors = F, row.names = NULL)
geneList_bulk = data_bulk[, 1]
data_matrix_bulk = as.matrix(data_bulk[, 2:ncol(data_bulk)])
data_matrix_bulk <- log2(data_matrix_bulk + 1)

extdir<-system.file("extdata",package="scSTAR2")
BRCA_bulk_prognosis_1<-paste(extdir,"/BRCA_bulk_prognosis.txt",sep="")
meta_bulk = read.table(BRCA_bulk_prognosis_1, sep = "\t", stringsAsFactors = F, row.names = NULL)
# patient_ID = table2cell(meta_bulk(:,1));
resp_pattern = meta_bulk[, 2]


idx_Good = which(resp_pattern == 'Good', arr.ind=T)
resp_pattern_binary = rep(0, length(resp_pattern))
resp_pattern_binary[idx_Good] <- 1#1 Good, 0 Poor
```

## load ST data with prognosis information
```
extdir<-system.file("extdata",package="scSTAR2")
BRCA_ST_H1<-paste(extdir,"/BRCA ST H1.txt",sep="")
data_sc = read.table(BRCA_ST_H1 ,sep = "\t", stringsAsFactors = F, row.names = NULL)
geneList_sc <- data_sc[, 1]
data_matrix_sc = as.matrix(data_sc[, 2:ncol(data_sc)])
data_matrix_sc <- log2(data_matrix_sc + 1)
```

## keep the shared genes between datasets
```
geneList <- intersect(geneList_bulk, geneList_sc)
geneList <- sort(geneList)
ia <- match(geneList, geneList_bulk)
ib <- match(geneList, geneList_sc)
data_matrix_bulk = data_matrix_bulk[ia, ]
data_matrix_sc = data_matrix_sc[ib, ]
data_matrix_sc = log2(data_matrix_sc + 1)
```

## gene filtering
This step can reduce the interferences of random noise by removing noise corrupted genes
```
par(mar=c(1,1,1,1))
idx_OGFSC = OGFSC(data_matrix_sc, nBins = 10, plot_option = 1)$OGFSC_idx; # OGFSC (Hao, et al. 2019)
idx_vector <- unlist(idx_OGFSC)
data_matrix_sc <- data_matrix_sc[idx_vector, , drop = FALSE]
idx_vectorr <- unlist(idx_OGFSC)
data_matrix_bulk <- data_matrix_bulk[idx_vectorr, , drop = FALSE]
idx_vector <- unlist(idx_OGFSC)
geneList_OGFSC_share <- geneList[idx_vector]
```
## PLS1 projection
In this step, the PLS1 method is used to remove batch effect bewteen ST and RNA-seq data
```
NCV <- 5 # Suggest to keep this unchanged，default
minNC <- 2 # Suggest to keep this unchanged，default
PLScomp <- 2 # Suggest to keep this unchanged
MODEL <- PLSconstruct(t(data_matrix_sc), t(data_matrix_bulk), 'mc', NCV, PLScomp, minNC)
```
## reconstruct ST
```
temp <- t(data_matrix_sc)
temp = temp-matrix(1,dim(temp)[1],1)%*%MODEL$mu;
temp = temp-temp%*%MODEL$XL%*%ginv(MODEL$XL);
temp = temp+matrix(1,dim(temp)[1],1)%*%MODEL$mu;
S_sc <- t(temp)
```
## reconstruct bulk
```
temp = t(data_matrix_bulk)
temp = temp-matrix(1,dim(temp)[1],1)%*%MODEL$mu;
temp = temp-temp%*%MODEL$XL%*%ginv(MODEL$XL);
temp = temp+matrix(1,dim(temp)[1],1)%*%MODEL$mu;
S_bulk = t(temp)
```

## meta prediction model training
```
PLScomp2 = 10 # by default, 10. Might be slightly adjusted to 8 or 9
FCV = 7 #Do not worry about this parameter for now
data_1 = S_bulk[, resp_pattern_binary == 1]
data_2 = S_bulk[, resp_pattern_binary == 0]
```
## train discriminatory model on bulk data
```
MODEL <- PLSconstruct(t(data_1), t(data_2), 'mc', NCV, PLScomp2, minNC)
```
## apply the model on ST data
```
temp <- t(S_sc)
temp = temp-matrix(1,dim(temp)[1],1)%*%MODEL$mu;
temp = temp%*%MODEL$XL%*%ginv(MODEL$XL);
temp = temp+matrix(1,dim(temp)[1],1)%*%MODEL$mu;
SS_sc <- t(temp)
```
## save results
```
geneList_OGFSC_share <- as.list(geneList_OGFSC_share)
resp_pattern <- as.list(resp_pattern)
scSTAR2_H1 = list(data_matrix_sc=data_matrix_sc, data_matrix_bulk=data_matrix_bulk, S_bulk=S_bulk, S_sc=S_sc, SS_sc=SS_sc, geneList_OGFSC_share=geneList_OGFSC_share, resp_pattern=resp_pattern)
names(scSTAR2_H1) = c('data_matrix_sc','data_matrix_bulk','S_bulk','S_sc','SS_sc','geneList_OGFSC_share','resp_pattern')
saveRDS(scSTAR2_H1,file = 'scSTAR2_H1.rds')
```

## load data
this part first clusters the spots based on the reconstructed data, then associates meta informaton with each spot cluster
```
data = readRDS("scSTAR2_H1.rds")
SS_sc <- data$SS_sc
SS_sc <- matrix(as.numeric(SS_sc),nrow=nrow(SS_sc),dimnames=dimnames(SS_sc))
genelist <- unlist(data$geneList_OGFSC_share)
rownames(SS_sc) <- genelist
colnames(SS_sc) <- paste0("Cell",c(1:dim(SS_sc)[2]))

data_matrix_sc <- data$data_matrix_sc
data_matrix_sc <- matrix(as.numeric(data_matrix_sc),nrow=nrow(data_matrix_sc),dimnames=dimnames(data_matrix_sc))
rownames(data_matrix_sc) <- genelist
colnames(data_matrix_sc) <- paste0("Cell",c(1:dim(data_matrix_sc)[2]))

S_sc <- data$S_sc
S_sc <- matrix(as.numeric(S_sc),nrow=nrow(S_sc),dimnames=dimnames(S_sc))
genelist <- unlist(data$geneList_OGFSC_share)
rownames(S_sc) <- genelist
colnames(S_sc) <- paste0("Cell",c(1:dim(S_sc)[2]))
```
## clustering and identify marker genes of each cluster
```
Integrated <- CreateSeuratObject(counts = SS_sc, project = "Integrated", min.cells = 0, min.features = 0)
Integrated <- NormalizeData(Integrated, verbose = FALSE)
Integrated <- FindVariableFeatures(Integrated, selection.method = "vst", nfeatures = 2000)
Integrated <- ScaleData(Integrated, verbose = FALSE)
Integrated <- RunPCA(Integrated, npcs = 20, verbose = FALSE)
Integrated <- RunTSNE(Integrated, reduction = "pca", dims = 1:20) 
Integrated <- FindNeighbors(Integrated, reduction = "pca", dims = 1:20)
Integrated <- FindClusters(Integrated, resolution = 0.5)
scS_cluster <- paste0('scS_', Integrated$seurat_clusters)

col_scSTAR2 <- c("#00A08A", "#E2D200", "#9986A5",  "#E58601","#46ACC8", "#B40F20")


tplot <- DimPlot(Integrated, reduction = "tsne", label=TRUE, pt.size=3, label.size = 5)+
  scale_color_manual(values=col_scSTAR2)+
  theme(legend.position="none")
tplot[[1]]$layers[[1]]$aes_params$alpha = .5

tiff('scSTAR2 scatter ST H1.tiff', units="in", width=5, height=5, res=300, compression = 'lzw')
tplot
dev.off()


Integrated.markers <- FindAllMarkers(Integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
write.table(Integrated.markers , "MarkerGenesByCluster_scSTAR2_H1.combined.xls",sep="\t")
```


## load bulk data
```
data_matrix_bulk <- data$data_matrix_bulk
data_matrix_bulk <- matrix(as.numeric(data_matrix_bulk),nrow=nrow(data_matrix_bulk),dimnames=dimnames(data_matrix_bulk))
genelist <- unlist(data$geneList_OGFSC_share)
rownames(data_matrix_bulk) <- genelist
colnames(data_matrix_bulk) <- paste0("Patient",c(1:dim(data_matrix_bulk)[2]))
```
## meta information
```
resp_pattern <- unlist(data$resp_pattern)
```
## load gene markers
```
markers <- read.table("MarkerGenesByCluster_scSTAR2_H1.combined.xls",sep="\t")
```
## Processing of gene marker data may involve extracting related genes by cluster. This part of the code involves selecting specific genes from the clustering information for subsequent analysis.
```
clusters = markers$cluster
clusters_uni = unique(clusters)
FC <- markers$avg_log2FC
selectedGenes = list()
for (i in 1:length(clusters_uni)){
  idx = which(clusters==clusters_uni[i])
  x <- sort(FC[idx], decreasing = TRUE, index.return = TRUE)
  selectedGenes[[i]] = as.vector(markers$gene[idx[x$ix]])
}
```


## This part evaluates genetic similarity between each cluster and bulk data with different meta information
```
buffer <- data$S_bulk
buffer <- matrix(as.numeric(buffer),nrow=nrow(buffer),dimnames=dimnames(buffer))
genelist <- unlist(data$geneList_OGFSC_share)
rownames(buffer) <- genelist

ngenes =200
n <- 5#This number can be changed
```
## For each cluster, extract the specified number of genes from the bulk data and the single-cell data of the corresponding cluster, then calculate the Pearson correlation between these genes, and save the results in the list corrMat_H1
```
corrMat_H1 <- list()

for (i in 0:n) {
  cluster_name <- paste0("scS_", i)
  X <- SS_sc[, which(scS_cluster == cluster_name)]
  selected_genes <- selectedGenes[[i + 1]][1:ngenes]
  Cor <- cor(buffer[selected_genes, ], X[selected_genes, ], method = "pearson")
  corrMat_H1[[paste0("Cor", i)]] <- Cor
}
corrMat_H1[["resp_pattern"]] <- as.list(resp_pattern)
names(corrMat_H1) <- c(paste0("Cor", 0:n), "resp_pattern")
saveRDS(corrMat_H1, file = 'corrMat_H1.rds')
```


## This code mainly implements the OPLS-DA (Partial Least Squares Discriminant Analysis) process
The OPLSDA function, encompassing model training, visualization, and significance testing based on specified parameters and data, is a crucial aspect of the scSTAR2R package. It's essential to note that the OPLSDA para.txt file resides in the inst directory of the package. Parameters within this file can be modified as needed. Specifically, scSTAR2 categorizes clusters 0, 4, and 5 as related to a good prognosis, while cluster 2 is identified as associated with a poor prognosis.
```
data = readRDS("corrMat_H1.rds")
pattern <- do.call(rbind, data$resp_pattern)
idx_R <- which(pattern == "Good", arr.ind = TRUE)
idx_NR <- which(pattern == "Poor", arr.ind = TRUE)

for (i in 0:n) {
  data1 = data[[paste0("Cor", i)]][idx_R[, 1], ]
  data2 = data[[paste0("Cor", i)]][idx_NR[, 1], ]
  cellindex = matrix(seq(1, ncol(data[[paste0("Cor", i)]])), nrow = 1, ncol = ncol(data[[paste0("Cor", i)]]))
  model <- OPLSDA(data1, data2,cellindex )
  ggsave(file = paste0('C', i, '.tiff'), dpi = 300, compression = 'lzw', width = 6, height = 5, units = "in")
  dev.off()
}
```

## predicted phenotype category plot
The main purpose of this code is to create, visualize and save the scatter plot of t-SNE analysis for the sample grouping (phenotype) information
of the Seurat object Integrated, which can be grouped according to the visualization results obtained by OPLSDA above.
```
pheno <- rep('Backgroud', length(scS_cluster))
pheno[c(which(scS_cluster=='scS_1'), which(scS_cluster=='scS_3'))] <- 'Backgroud'
pheno[c(which(scS_cluster=='scS_0'),
        which(scS_cluster=='scS_4'), which(scS_cluster=='scS_5'))] <- 'Good'
pheno[c(which(scS_cluster=='scS_2'))] <- 'Poor'
Integrated$pheno <- pheno

tplot <- DimPlot(Integrated, reduction = "tsne", group.by = "pheno", pt.size=3, cols = c("#9986A5", "#00A08A","#B40F20")) #wes_palette("Moonrise3")
tplot[[1]]$layers[[1]]$aes_params$alpha = .5

tiff('scSTAR2 scatter ST H1 phenotype.tiff', units="in", width=6, height=5, res=300, compression = 'lzw')
tplot
dev.off()
```



## spatial plot to highlight associated spots
```
extdir <- system.file("extdata", package = "scSTAR2")
H1_seurat_1 <- file.path(extdir,  "H1_seurat.rds")
STObj <- readRDS(H1_seurat_1)
STObj <- SCTransform(STObj, assay = "Spatial", ncells = 3000, verbose = FALSE)
STObj <- RunPCA(STObj)
STObj <- RunUMAP(STObj, dims = 1:20)
STObj <- FindNeighbors(STObj, dims = 1:20)
STObj <- FindClusters(STObj, resolution = 1, verbose = FALSE)


STObj$scS_cluster <- scS_cluster
col_scSTAR2 <- c("#E2D200", "#00A08A", "#9986A5", "#46ACC8",  "#E58601","#B40F20")

tiff('scSTAR2 H1 ST.tiff', units="in", width=6, height=6, res=300, compression = 'lzw')
FeatureOverlay(STObj,features = "scS_cluster", pt.size = 4,type ="raw", cols=col_scSTAR2)
dev.off()
```
