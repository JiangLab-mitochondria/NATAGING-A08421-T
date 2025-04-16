library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
set.seed(1234)

## four sample merged and clustering
TE_2852_75 <- readRDS("resutls/scRNA/seurat-results/clustering-single-sample/TE_2852_75_clean_celltype.rds")
TE_2900_0 <- readRDS("resutls/scRNA/seurat-results/clustering-single-sample/TE_2900_0_clean_celltype.rds")
TE_75_75 <- readRDS("resutls/scRNA/seurat-results/clustering-single-sample/TE_75_75_clean_celltype.rds")
TE_742_0 <- readRDS("resutls/scRNA/seurat-results/clustering-single-sample/TE_742_0_clean_celltype.rds")

all.merged <- merge(TE_2852_75, c(TE_2900_0, TE_75_75, TE_742_0),
                    add.cell.ids = c("TE_2852_75", "TE_2900_0",
                                     "TE_75_75", "TE_742_0"
                    ))
head(all.merged)
all.merged <- NormalizeData(all.merged)
library(stringr)
set.seed(1234)
g2m.genes <- str_to_title(tolower(cc.genes$g2m.genes))
s.genes <- str_to_title(tolower(cc.genes$s.genes))
head(all.merged)
all.merged <- CellCycleScoring(all.merged, s.features = s.genes, g2m.features = g2m.genes)
head(all.merged)

all.merged <- FindVariableFeatures(all.merged, selection.method = "vst", nfeatures = 2000)
all.merged <- ScaleData(all.merged)
all.merged <- RunPCA(all.merged)
all.merged <- RunUMAP(all.merged, reduction = "pca", dims = 1:30)
all.merged <- FindNeighbors(all.merged, reduction = "pca", dims = 1:30)
all.merged <- FindClusters(all.merged, resolution = 1.2)

p1 <- DimPlot(all.merged, reduction = "umap", label = TRUE) +
  NoLegend()
p1
p2 <- DimPlot(all.merged, reduction = "umap", label = F, group.by = "orig.ident")
p2

p3 <- DimPlot(all.merged, reduction = "umap", label = TRUE, group.by = "celltype") +
  NoLegend()
p3


all.merged <- ScaleData(all.merged, vars.to.regress = c("S.Score", "G2M.Score"))
all.merged <- RunPCA(all.merged)
all.merged <- RunUMAP(all.merged, reduction = "pca", dims = 1:30)
all.merged <- FindNeighbors(all.merged, reduction = "pca", dims = 1:30)
all.merged <- FindClusters(all.merged, resolution = 1.2)

p4 <- DimPlot(all.merged, reduction = "umap", label = TRUE) +
  NoLegend()
p4
p5 <- DimPlot(all.merged, reduction = "umap", label = F, group.by = "orig.ident")
p5

p6 <- DimPlot(all.merged, reduction = "umap", label = TRUE, group.by = "celltype") +
  NoLegend()
p6


saveRDS(all.merged, file = "resutls/scRNA/seurat-results/clustering-merged-samples/four_samples_merged_celltype_cellcycle_regress.rds")

## DE genes detection
cells <- readRDS("resutls/scRNA/seurat-results/clustering-merged-samples/four_samples_merged_celltype_cellcycle_regress.rds")

de_genes <- function(cell_type) {
  cell_type <- as.character(cell_type)
  celltype <- subset(cells, idents = cell_type)
  print(levels(celltype))
  Idents(object = celltype) <- "sample.type"
  print(levels(celltype))
  markers <- FindMarkers(celltype, ident.1 = "mut", ident.2 = "wt")
  markers$gene.symbol <- rownames(markers)
  return(markers)
}

levels(cells)
## Endo
marker.genes <- de_genes(cell_type = "Endo")
write.table(marker.genes,
            file="resutls/scRNA/seurat-results/FindMarkers-results/Endo-mut-vs-wt-DEgenes.txt",
            sep="\t", quote=F, row.names = F, col.names = T)
