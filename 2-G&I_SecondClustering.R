#Powered by RainShalder
#Excecute this using an annotated Seurat object.
#### Environment Setup ####
systemInformation <- Sys.info()
operationSystem <- systemInformation[1]
setwd("***")

library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(harmony)
library(clustree)
library(dplyr)
library(tidyverse)
library(future)
library(glmGamPoi)
library(presto)

options(future.globals.maxSize = 50 * 1000 * 1024^2)

#### Original Data Input ####
renal_fc <- readRDS("***/renal_anno.rds")
Idents(renal_fc) <- "celltype"
sc_active <- renal_fc[, Idents(renal_fc) == "Glo & Inter"]
DefaultAssay(sc_active) <- "RNA"
sc_active[['SCT']] <- NULL
sc_active[['pca']] <- NULL
sc_active[['harmony']] <- NULL
sc_active[['umap']] <- NULL
sc_active <- JoinLayers(sc_active)

#### SCTransfrom, PCA & Harmony ####
#SCTransform
sc_SCT <- SCTransform(sc_active,
                      vars.to.regress = c('subject', 'percent.mt', 'percent.HB'))
DefaultAssay(sc_SCT) <- "SCT"
Idents(sc_SCT) <- "subject"
top10 <- head(VariableFeatures(sc_SCT), 10)
LabelPoints(plot = VariableFeaturePlot(object = sc_SCT), points = top10, 
            repel = T)

#PCA
sc_SCT <- RunPCA(sc_SCT, verbose = F)
DimPlot(sc_SCT, reduction = "pca")
ElbowPlot(sc_SCT, ndims = 50)
pct <- sc_SCT[["pca"]]@stdev / sum(sc_SCT[["pca"]]@stdev) * 100
cumsum(pct)
pcs = 1:40

#Harmony
sc_integ <- sc_SCT
sc_integ[["RNA"]] <- split(sc_integ[["RNA"]], f = sc_integ$subject)
sc_integ <- RunHarmony(sc_integ, group.by.vars = "subject", 
                       assay.use = "SCT")
table(sc_integ$subject)
DimPlot(sc_integ, reduction = "harmony")
sc_integ <- FindNeighbors(sc_integ, reduction = "harmony", dims = pcs)

#### Second-Level Clustering ####
sc_clutest <- sc_integ
res_seq <- seq(0.1, 1, by = 0.1)
for(res in res_seq){
  sc_clutest <- FindClusters(sc_clutest, resolution = res)
}
clustree(sc_clutest, prefix = 'SCT_snn_res.') + coord_flip()
sc_integ <- FindClusters(sc_integ, resolution = 0.5)
sc_integ <- RunUMAP(sc_integ, reduction = "harmony", dims = pcs)
DimPlot(sc_integ, reduction = "umap", label = T)
DimPlot(sc_integ, reduction = "umap", label = T, group.by = "subject")

saveRDS(sc_integ, "***/sc_G&I.rds")

#### Secondary-Cluster Annotation ####
#sc_integ <- readRDS("***/sc_G&I.rds")
#FindAllMarkers
markers <- FindAllMarkers(sc_integ, min.pct = 0.1, logfc.threshold = 0.5, 
                          only.pos = T, test.use = "wilcox")
markers_df <- markers %>% 
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)
write.table(markers_df, "***",
            quote = F, sep = "\t", row.names = F, col.names = T)

#Second-level cellmarkers
sc_anno <- sc_integ
Anno_Markers <- c("PODXL",  #ALL
                  "PECAM1",  #Endo
                  "GJA4",  #LA, AEA, DVR
                  "SEMA3G", "SERPINE2", "FBLN2", "VWF",  #LA
                  "SLC6A6",  #AEA
                  "SCIN",  #DVR
                  "RGCC",  #PTC1/2
                  "NRP2",  #AVR
                  "NTN4", "EHD3",  #GEC
                  "PDGFRB", "GATA3",  #Mesangial Cells
                  "CLDN1", "PAX8",  #PEC
                  "WT1", "NPHS1", "NPHS2",  #Podo
                  "MUSTN1", "ACTA2",  #VSMC
                  "LUM", "DCN", "PDGFRA",  #Fibro/Myofibro
                  "CSPG4",  #Peri
                  "MIOX",  #Tub(Pollution)
                  "PTPRC"  #Immuno(Pollution)
)
DotPlot(sc_anno, features = unique(Anno_Markers), group.by = "seurat_clusters") + 
  RotatedAxis() + 
  scale_x_discrete("") +
  scale_y_discrete("")
FeaturePlot(sc_anno, reduction = "umap", 
            features = Anno_Markers)
DimPlot(sc_anno, reduction = "umap", split.by = "condition")

#Second-level annotation
sc_anno$celltype <- recode(sc_anno$SCT_snn_res.0.5,
                           "0" = "PTC2",
                           "1" = "PTC1") # Example
Idents(sc_anno) <- "celltype"
sc_anno <- sc_anno[, Idents(sc_anno) != "NA"]
sc_anno <- sc_anno[, Idents(sc_anno) != "Tub"]
sc_anno <- sc_anno[, Idents(sc_anno) != "Immuno"]
table(sc_anno$celltype)
DimPlot(sc_anno, reduction = "umap", label = T)
DotPlot(sc_anno, features = unique(Anno_Markers), group.by = "celltype") + 
  RotatedAxis() + 
  scale_x_discrete("") +
  scale_y_discrete("")
FeaturePlot(sc_anno, reduction = "umap", 
            features = Anno_Markers, label = T, label.size = 2)

#Data postprocess
sc_anno@meta.data[1:5,]
sc_anno[['percent.mt']] <- NULL
#sc_anno[['percent.hsp']] <- NULL
sc_anno[['percent.HB']] <- NULL
sc_anno[['SCT_snn_res.0.2']] <- NULL
sc_anno[['SCT_snn_res.0.5']] <- NULL
#sc_anno[['SCT_snn_res.1']] <- NULL
#sc_anno[['DS.tri_class']] <- NULL
#sc_anno[['predicted_doublets']] <- NULL
DefaultAssay(sc_anno) <- "RNA"
sc_anno[['SCT']] <- NULL
sc_anno <- JoinLayers(sc_anno)
sc_anno@meta.data[1:5,]

#Construct groups
sc_anno@meta.data$group_celltype <- paste(sc_anno$condition, sc_anno$celltype,
                                          sep = "_")
sc_anno@meta.data[1:5,]

saveRDS(sc_anno, "***/sc_G&I_anno.rds")
sc_anno <- readRDS("***/sc_G&I_anno.rds")

