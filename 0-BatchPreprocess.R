#Powered by RainShalder
systemInformation <- Sys.info()
operationSystem <- systemInformation[1]
setwd("***")

library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(harmony)
library(clustree)
library(dplyr)
library(tidyverse)
library(future)
library(glmGamPoi)
library(presto)

options(future.globals.maxSize = 50 * 1000 * 1024^2)

DataFiltering <- function(scdat, nGSM){
  HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  HB.genes <- CaseMatch(HB.genes, rownames(scdat))
  scdat[["percent.mt"]] <- PercentageFeatureSet(scdat, pattern = "^MT-")
  scdat[["percent.HB"]] <- PercentageFeatureSet(scdat, features = HB.genes)
  setwd("***/QC_Pic")
    picpath1 <- paste0(nGSM, "Before.png")
    png(picpath1)
    a <- VlnPlot(scdat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"))
    print(a)
    dev.off()
  minGene = 200; maxGene = quantile(scdat$nFeature_RNA, 0.99)
  minCount = 500; maxCount = quantile(scdat$nCount_RNA, 0.99)
  pctMT = 25; pctHB = 1
  scdat <- subset(scdat, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
                  nCount_RNA > minCount & nCount_RNA < maxCount &
                  percent.mt < pctMT & percent.HB < pctHB)
    picpath2 <- paste0(nGSM, "QC.png")
    png(picpath2)
    a <- VlnPlot(scdat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"))
    print(a)
    dev.off()
  setwd("***")
  return(scdat)
}

scRNA_3Steps_PCA <- function(scdat, nGSM){
  #### Three Steps ####
  scdat <- NormalizeData(scdat, normalization.method = "LogNormalize",
                            scale.factor = 10000)
  scdat <- FindVariableFeatures(scdat, selection.method = "vst", nfeatures = 2000)
  scdat <- ScaleData(scdat, features = rownames(scdat))
  
  #### PCA Analysis ####
  scdat <- RunPCA(scdat, features = VariableFeatures(object = scdat))
  
  scdat <- FindNeighbors(scdat, dims = pcs)
  scdat <- FindClusters(scdat, resolution = 0.5)
  head(Idents(scdat), 5)
  scdat <- RunUMAP(scdat, dims = pcs)
  DimPlot(scdat, reduction = "umap")
  setwd("***/QC_Pic")
    picpath <- paste0(nGSM, "_4Features.png")
    png(picpath)
    a <- FeaturePlot(scdat, reduction = "umap", features = c("EIF2AK2", "TNFSF10", "CASP1", "EEF2"))
    print(a)
    dev.off()
  setwd("***")
  return(scdat)
}

FindDoublet <- function(scdat, nGSM){
  sweep.res.list <- paramSweep(scdat, PCs = pcs)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  sweep.stats[order(sweep.stats$BCreal), ]
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  #Doublet Rate Calculation, or Search in the Table
  homotypic.prop <- modelHomotypic(scdat$seurat_clusters)
  DoubletRate <- ncol(scdat) * 8*1e-6
  nExp_poi <- round(DoubletRate * nrow(scdat@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  #Doublet Find
  scdat <- doubletFinder(scdat, PCs = pcs, pN = 0.25, pK = pK_bcmvn, 
                               nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
  DF.c_num <- grep("DF.classifications", colnames(scdat@meta.data))
  DF.c <- colnames(scdat@meta.data)[DF.c_num]
  setwd("***/QC_Pic")
    picpath <- paste0(nGSM, "_Doublet.png")
    png(picpath)
    a <- DimPlot(scdat, reduction = "umap", group.by = DF.c)
    print(a)
    dev.off()
  setwd("***")
  Idents(scdat) <- DF.c
  scdat <- subset(scdat, idents = "Singlet")
}

main{
  Metadata <- read.csv("***/Metadata.csv")  #A csv file including: ID, condition, and other needed info
  GSMlist <- Metadata$GSM
  pcs <- 1:40
  for(i in 1:length(GSMlist)){
    datapath <- paste0("***/", GSMlist[i])   #CellRanger files named "ID"
    datapath
    Disease.data <- Read10X(datapath)
    Disease <- CreateSeuratObject(counts = Disease.data, project = "Disease", 
                              min.cells = 3, 
                              min.features = 200)
    Disease <- DataFiltering(Disease, GSMlist[i])
    Disease <- scRNA_3Steps_PCA(Disease, GSMlist[i])
    Disease <- FindDoublet(Disease, GSMlist[i])
    #### Post Process ####
    Disease$condition <- Metadata$Condition[i]
    Disease$orig.ident <- factor(GSMlist[i])
    Disease$age <- Metadata$Age[i]
    Disease$subject <- Metadata$Subject[i]
    Disease$tissue <- "PBMC"
    Idents(Disease) <- "orig.ident"
    outpath <- paste0("***", GSMlist[i], "_QC.rds")
    write_rds(Disease, outpath)
  }
}
