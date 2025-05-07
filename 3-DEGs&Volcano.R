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
library(MAST)
library(DESeq2)
library(ggplot2)
library(ggrepel)

options(future.globals.maxSize = 50 * 1000 * 1024^2)

#### Data Preprocess ####
sc_anno <- readRDS("***")

#### Differential Expression Genes Analysis ####
Idents(sc_anno) <- "group_celltype"

#FindMarkers
#For once execution, here
table(sc_anno$group_celltype)
mks1 <- FindMarkers(sc_anno, ident.1 = "1", ident.2 = "2",
                    test.use = "MAST", slot = "counts")
write.csv(mks1, "***")

#For loop execution, here
c1 <- "Disease"
c2 <- "Control"
clist <- c("1", "2", "3")

for(i in clist){
  id1 <- paste0(c1, "_", i)
  id2 <- paste0(c2, "_", i)
  writepath <- paste0("***", i, "_", c1, "vs", c2, ".csv")
  mks <- FindMarkers(sc_anno, ident.1 = id1, ident.2 = id2,
                     test.use = "MAST", slot = "counts")
  write.csv(mks, writepath)
}

#### Volcano Plot ####
dat <- read.csv("***")
length <- nrow(dat)
for (i in 1:length){
  dat[i, "FDR"] <- signif(dat[i, "p_val"] * length / i, 4)
}
cutoff_fdr <- 0.05
cutoff_log2fc <- 0.5
dat <- dat[!grepl("ENSG", dat$X), ]
dat <- dat[!grepl("RPL", dat$X), ]
dat$change <- ifelse(dat$FDR < cutoff_fdr & abs(dat$avg_log2FC) > cutoff_log2fc,
                     ifelse(dat$avg_log2FC > cutoff_log2fc, 'Up', 'Down'),
                     'Stable')
table(dat$change)
dat_DEGs <- dat[dat$change == "Up" | dat$change == "Down", ]
write.csv(dat_DEGs, "***",
          row.names = F)

p <- ggplot(dat, aes(x = avg_log2FC, y = -log10(FDR), colour = change)) +
  geom_point(size = 1.0) +
  scale_color_manual(values = c("#546de5", "#d2dae2", "#ff4757")) +
  geom_vline(xintercept = c(-1, 1), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(cutoff_fdr), lty = 4, col = "black", lwd = 0.8) +
  labs(x="log2FC",
       y="-log10(FDR)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
dat$label <- ifelse((dat$FDR < 1e-10 & dat$avg_log2FC >= 2 |
                       dat$FDR < 1e-80 & dat$avg_log2FC <= -2),
                    as.character(dat$X), "")
p + geom_text_repel(data = dat, aes(x = avg_log2FC, y = -log10(FDR), label = label),
                    size = 3, box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"),
                    segment.color = "black", 
                    show.legend = FALSE,
                    max.overlaps = 100)
