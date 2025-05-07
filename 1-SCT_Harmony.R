#Powered by RainShalder
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
dat <- readRDS("***")
table(dat$subject)
renal <- dat
Idents(renal) <- "condition"
renal <- NormalizeData(renal, normalization.method = "LogNormalize",
                       scale.factor = 10000)
table(renal$subject)

#### Quality Control ####
renal[["percent.mt"]] <- PercentageFeatureSet(renal, pattern = "^MT-")
HB.genes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(renal))
renal[["percent.HB"]] <- PercentageFeatureSet(renal, features = HB.genes)
#FeatureScatter(renal, "nCount_RNA", "percent.mt", group.by = "subject")
#FeatureScatter(renal, "nCount_RNA", "percent.HB", group.by = "subject")
#FeatureScatter(renal, "nCount_RNA", "nFeature_RNA", group.by = "subject")
#Violin
theme.set2 <- theme(axis.title.x = element_blank())
plot.features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
group <- "subject"
plots <- list()
for (i in seq_along(plot.features)){
  plots[[i]] <- VlnPlot(renal, group.by = group, pt.size = 0, 
                        features = plot.features[i]) + theme.set2 + NoLegend()
}
wrap_plots(plots = plots, nrow = 2)

#### Normalization - SCTransform ####
renal_SCT <- SCTransform(renal,
                         vars.to.regress = c('subject', 'percent.mt', 'percent.HB'))
top10 <- head(VariableFeatures(renal_SCT), 10)
LabelPoints(plot = VariableFeaturePlot(object = renal_SCT), points = top10, 
            repel = T)

#### PCA ####
DefaultAssay(renal_SCT) <- "SCT"
Idents(renal_SCT) <- "subject"
renal_SCT <- RunPCA(renal_SCT, verbose = F)
DimPlot(renal_SCT, reduction = "pca")
ElbowPlot(renal_SCT, ndims = 50)
pct <- renal_SCT[["pca"]]@stdev / sum(renal_SCT[["pca"]]@stdev) * 100
cumsum(pct)
pcs = 1:40

#### Harmony ####
renal_integ <- renal_SCT
renal_integ[["RNA"]] <- split(renal_integ[["RNA"]], f = renal_integ$subject)
renal_integ <- RunHarmony(renal_integ, group.by.vars = "subject", 
                          assay.use = "SCT")
table(renal_integ$subject)
DimPlot(renal_integ, reduction = "harmony")
renal_integ <- FindNeighbors(renal_integ, reduction = "harmony", dims = pcs)

#### Cell Clustering ####
renal_clutest <- renal_integ
res_seq <- seq(0.1, 1, by = 0.1)
for(res in res_seq){
  renal_clutest <- FindClusters(renal_clutest, resolution = res)
}
cp1 <- clustree(renal_clutest, prefix = 'SCT_snn_res.') + coord_flip()
cp1
renal_integ <- FindClusters(renal_integ, resolution = 0.2)
renal_integ <- RunUMAP(renal_integ, reduction = "harmony", dims = pcs)
DimPlot(renal_integ, reduction = "umap", label = T, raster = FALSE)
DimPlot(renal_integ, reduction = "umap", label = T, split.by = "subject")

saveRDS(renal_integ, "***/renal_integ.rds")

#### Cluster Annotation ####
#FindAllMarkers
#renal_integ <- readRDS("***/renal_integ.rds")
markers <- FindAllMarkers(renal_integ, min.pct = 0.1, logfc.threshold = 0.5, 
                          only.pos = T, test.use = "wilcox")
markers_df <- markers %>% 
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)
write.table(markers_df, "***",
            quote = F, sep = "\t", row.names = F, col.names = T)
renal_anno <- renal_integ

#First-level cellmarkers
Immuno <- c("PTPRC")
#Myeloid
Mye <- c("CD14", "CD16", "CD68",  #Monocyte
         "FCGR2B", "CD1C",  #cDC
         "CD163", "MRC1",  #M2 Macrophage (& CD68)
         "CD86", "TREM2",  #Lipid Macrophage
         "TPSB2",  "CPA3",  #Mast
         "FCGR3B", "MNDA", "CXCL8", "SELL"  #Neutrophil
)
#T cells
TNK <- c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "CD25", "FOXP3",  #T cells
         "KLRD1", "NKG7")  #NK cells
Bc <- c("CD19", "CD20", "CD27", "CD138")  #B cells
#Tubular
Tub <- c("MIOX", "SLC22A8", "VCAM1",  #PT
         "SLC12A1",  #LOH
         "SLC8A1", "SLC12A3",  #DT
         "ATP6V0D2", "AQP2", "CA2")  #CD
#Glomerulus & Interstitium
GEndo <- c("PODXL", "PECAM1", "KDR", "EHD3", "SEMA3G")
Inter <- c("GATA3", "PDGFRB",  #Mes
           "CLDN1", "PAX8",  #PECs
           "CSPG4", "NG2",  #Peri
           "WT1", "NPHS1", "NPHS2",  #Podo
           "LUM", "DCN", "PDGFRA",  #Fibro/Myofibro
           "MUSTN1", "ACTA2")  #VSMC
Anno_Markers <- c(Immuno, Mye, TNK, Bc, Tub, GEndo, Inter)
DotPlot(renal_anno, features = unique(Anno_Markers), group.by = "seurat_clusters") + 
  RotatedAxis() + 
  scale_x_discrete("") +
  scale_y_discrete("")
FeaturePlot(renal_anno, reduction = "umap",
            features = c(Immuno, Mye, TNK, Bc))
FeaturePlot(renal_anno, reduction = "umap",
            features = c(Tub, GEndo, Inter))

#First-Level Annotation
renal_anno$celltype <- renal_anno$SCT_snn_res.0.2  # Change with the value of resolution
table(renal_anno$celltype)
renal_anno$celltype <- recode(renal_anno$celltype,
                              "0" = "Tubule",
                              "1" = "T & NK",
                              "2" = "Glo & Inter",
                              "3" = "Myeloid",
                              "4" = "B")  # Example
Idents(renal_anno) <- "celltype"

#### Visualization ####
#Dimplot
DimPlot(renal_anno, reduction = "umap", group.by = "celltype", 
        label = T, label.size = 6, raster = FALSE)
renal_anno <- renal_anno[, Idents(renal_anno) != "Ribosome"]
DimPlot(renal_anno, reduction = "umap", group.by = "celltype", 
        label = T, label.size = 6, raster = FALSE, cols = color5) # 8*8 inches
DotPlot(renal_anno, features = unique(Anno_Markers), group.by = "celltype") + 
  RotatedAxis() + 
  scale_x_discrete("") +
  scale_y_discrete("")
saveRDS(renal_anno, "***/renal_anno.rds")
#renal_anno <- readRDS("***/renal_anno.rds")

