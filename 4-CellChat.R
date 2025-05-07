#Powered by RainShalder
#Excecute this using an annotated Seurat object.
#### Environment Setup ####
systemInformation <- Sys.info()
operationSystem <- systemInformation[1]
setwd("***")

library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(future)
library(parallel)
library(ComplexHeatmap)
library(uwot)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 40 * 1000 * 1024^2)

CCCompute <- function(cellchat){
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  #cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  #cellchat <- filterCommunication(cellchat, min.cells = 5)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  return(cellchat)
}


#### Data Preprocess ####
sc_anno <- readRDS("***")
sc_chat <- sc_anno
Idents(sc_chat) <- "celltype"
sc_chat <- subset(sc_chat, idents = c("1", "2"))
sc_chat$celltype <- as.factor(as.character(sc_chat$celltype))
Idents(sc_chat) <- "condition"
sc_Disease <- sc_chat[, Idents(sc_chat) == "Disease"]
sc_Control <- sc_chat[, Idents(sc_chat) == "Control"]

cc_Disease <- createCellChat(sc_Disease@assays$RNA$data, meta = sc_Disease@meta.data, group.by = "celltype")
cc_Control <- createCellChat(sc_Control@assays$RNA$data, meta = sc_Control@meta.data, group.by = "celltype")
#save(cc_Disease, cc_Control, file = "***/cc_Disease.Ctrl.rda")

cc_Disease <- CCCompute(cc_Disease)
cc_Control <- CCCompute(cc_Control)
#saveRDS(cc_Disease, "***/cc_Disease.rds")
#saveRDS(cc_Control, "***/cc_Control.rds")
cc_Disease <- readRDS("***/cc_Disease.rds")
cc_Control <- readRDS("***/cc_Control.rds")

cc_list <- list(Control = cc_Control, Disease = cc_Disease)
cellchat <- mergeCellChat(cc_list, add.names = names(cc_list), cell.prefix = FALSE)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, umap.method = "uwot", type = "functional")
cellchat <- netClustering(cellchat, type = "functional", do.parallel = FALSE)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, umap.method = "uwot", type = "structural")
cellchat <- netClustering(cellchat, type = "structural", do.parallel = FALSE)
p <- rankSimilarity(cellchat, type = "structural") + ggtitle("Structural similarity of pathway")
ggsave("Pathway_Similarity.pdf", p, width = 8, height = 5)

pathway.union <- union(cc_list[[1]]@netP$pathways, cc_list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cc_list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cc_list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cc_list[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cc_list)[2], width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

levels(cellchat@idents$joint)
netVisual_bubble(cellchat, sources.use = c(3, 10), targets.use = c(3, 10),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(3, 10), targets.use = c(3, 10), comparison = c(1, 2), 
                 max.dataset = 2, title.name = "Increased signaling in Disease", angle.x = 45, remove.isolate = T)
netVisual_bubble(cellchat, sources.use = c(3, 10), targets.use = c(3, 10), comparison = c(1, 2), 
                 max.dataset = 1, title.name = "Decreased signaling in Disease", angle.x = 45, remove.isolate = T)

netVisual_bubble(cellchat, sources.use = c(3, 10), targets.use = c(3, 10), comparison = c(1, 2), 
                 signaling = "***",
                 max.dataset = 1, title.name = "Decreased signaling in Disease", angle.x = 45, remove.isolate = T)
plotGeneExpression(cellchat, signaling = "***", angle.x = 45, idents = c("1", "2"), 
                   slot = "counts", split.by = "condition", color.use = c("lightblue", "red"))
