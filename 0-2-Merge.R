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

options(future.globals.maxSize = 500 * 1000 * 1024^2)

PostProcess <- function(s, orig, proj){    #s - SeuratObject, orig - String for orig.ident, proj - String for project.name
  greplist <- c("percent.mt", "percent.HB", "RNA_snn", "pANN", "DF.class", "clusters")
  maxL <- length(greplist)
  for(i in 1:maxL){
    erase.num <- grep(greplist[i], colnames(s@meta.data))
    erase.ident <- colnames(s@meta.data)[erase.num]
    s[[erase.ident]] <- NULL
  }
  s@graphs[] <- NULL
  s@reductions[] <- NULL
  s@commands[] <- NULL
  s$orig.ident <- factor(orig)
  Idents(s) <- "subject"
  s <- merge(s, add.cell.ids = orig, project = proj)
  return(s)
}

#### Post Process ####
dat1 <- readRDS("***/dat1.rds")
dat2 <- readRDS("***/dat2.rds")
#dat n

dat1 <- PostProcess(dat1, "ID1", "Subject1")
dat2 <- PostProcess(dat2, "ID2", "Subject2")
#dat n

#### Merge ####
Merged <- merge(dat1, y = c(dat2, #to dat n
                            ), project = "ProjectID")
table(Merged$orig.ident)
table(Merged$condition)
table(Merged$subject)
Merged <- JoinLayers(Merged)

#### Erase Scale.data ####
idents.df <- data.frame("orig.ident" = Merged$orig.ident, 
                        "nCount_RNA" = Merged$nCount_RNA, 
                        "nFeature_RNA" = Merged$nFeature_RNA, 
                        "subject" = Merged$subject, 
                        "condition" = Merged$condition,
                        "age" = Merged$age,
                        "tissue" = Merged$tissue) #Add other metadata
export <- CreateSeuratObject(counts = LayerData(Merged, assay = "RNA", layer = "counts"),
                             project = "Merged", 
                             meta.data = idents.df)
export@reductions$pca <- Merged@reductions$pca
VariableFeatures(export) = VariableFeatures(Merged)
Idents(export) = "orig.ident"
write_rds(export, "***/Merged.rds")

#Caution: Output Seurat object does not have "data" layer, should run Normalization.
