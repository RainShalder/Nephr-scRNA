if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")   #For R 4.4.0(2024-06) and above
BiocManager::version()
install.packages("remotes")
install.packages("devtools")
library(remotes)  #install github packages
library(devtools)  #install github packages

install.packages("RColorBrewer")
install.packages("tidyverse")
install.packages("xlsx")

install.packages("future")
install.packages("foreach")
install.packages("doParallel")

install.packages("Seurat")
BiocManager::install('glmGamPoi')  #Accelerate fit linear models in scRNA-Seq
BiocManager::install("DoubletFinder")
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')  #Backup route
install.packages("clustree")
install_github("immunogenomics/harmony")  #Harmony integration
install_github('immunogenomics/presto')  #Accelerate FindAllMarkers
#install_local("presto-master.zip")  #Install local package
BiocManager::install("MAST")

BiocManager::install("scater")  #For pseudobulk
install.packages("grr")
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)
BiocManager::install("DESeq2")
install.packages("magrittr")

BiocManager::install("clusterProfiler")
getOption('timeout')
options(timeout=1000)
BiocManager::install("org.Hs.eg.db")
install.packages("msigdbr")
BiocManager::install("GSVA")
BiocManager::install("GSEABase")
BiocManager::install("limma")
BiocManager::install("BiocParallel")

BiocManager::install("ComplexHeatmap")
BiocManager::install("BiocNeighbors")
devtools::install_github("sqjin/CellChat")
install.packages("NMF")

devtools::install_github("mojaveazure/seurat-disk")

