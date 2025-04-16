library(ArchR)
addArchRGenome("mm10")
addArchRThreads(threads = 8)
library(dplyr)
library(tidyr)
set.seed(1234)
setwd("/mnt/transposon1/zhangyanxiaoLab/chaiguoshi/projects/jiangmin-lab-project/work/zhang-le-ping-project/resutls/mtscATAC/archr-results")

# Function definition
archr_clustering <- function(sample_vector, sample_name) {
  
  sample_name <- as.character(sample_name)
  ArrowFiles <- sample_vector
  projMB1 <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    copyArrows = F)
  
  projMB1.IterativeLSI <- addIterativeLSI(
    ArchRProj = projMB1,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 8, 
    clusterParams = list( 
      resolution = c(0.2), 
      sampleCells = 10000, 
      n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force = T
  )
  
  projMB1.IterativeLSI.addclusters <- addClusters(
    input = projMB1.IterativeLSI,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force = T
  )
  projMB1.IterativeLSI.addclusters.addumap <- addUMAP(
    ArchRProj = projMB1.IterativeLSI.addclusters, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = T
  )
  projMB1.IterativeLSI.addclusters.addumap.addtsne <- addTSNE(
    ArchRProj = projMB1.IterativeLSI.addclusters.addumap, 
    reducedDims = "IterativeLSI", 
    name = "TSNE", 
    perplexity = 30,
    force = T
  )
  saveArchRProject(ArchRProj = projMB1.IterativeLSI.addclusters.addumap.addtsne, 
                   outputDirectory = paste0("clustering/", sample_name), 
                   load = TRUE,
                   dropCells = TRUE
  )
  return()
}