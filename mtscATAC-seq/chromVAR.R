## kidney-3W-TE-G14102A-4475-80

library(chromVAR)
library(motifmatchr)
library(JASPAR2022)
library(Matrix)
library(SummarizedExperiment)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)
library(tidyr)
set.seed(1234)

setwd("/mnt/transposon1/zhangyanxiaoLab/chaiguoshi/projects/jiangmin-lab-project/work/zhang-le-ping-project/")


# peakfile <- "resutls/mtscATAC/cellranger-atac-results/kidney-3W-TE-G14102A-4475-80/outs/filtered_peak_bc_matrix/peaks.bed"
# peaks <- read.csv(peakfile, sep = "\t", header = F)
# table(peaks$V1)
# peaks <- subset(peaks, V1 %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
#                                  "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
#                                  "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
#                                  "chr19", "chrX", "chrY"
#                                  ))
# 
# table(peaks$V1)
# write.table(peaks,
#             file="resutls/mtscATAC/cellranger-atac-results/kidney-3W-TE-G14102A-4475-80/outs/filtered_peak_bc_matrix/kidney-3W-TE-G14102A-4475-80_peaks_filtered.bed",
#             sep="\t", quote=F, row.names = F, col.names = F)




peakfile <- "resutls/mtscATAC/cellranger-atac-results/kidney-3W-TE-G14102A-4475-80/outs/filtered_peak_bc_matrix/kidney-3W-TE-G14102A-4475-80_peaks_filtered.bed"
peaks <- getPeaks(peakfile, sort_peaks = TRUE)
peaks
peaks <- resize(peaks, width = 600, fix = "center")
peaks

bamfile <- c(
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/B-cells_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/CD_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/DCT_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/Dendritic-cells-11b_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/Endothelium_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/LOH_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/Macrophages_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/Neutrophils-1_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/Neutrophils-2_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/NK-cells_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/Pan-PT_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/PCT_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/Podocyte_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/PST_sortCoord_rename_RG.bam",
  "resutls/mtscATAC/sinto-results/kidney-3W-TE-G14102A-4475-80/cell_type_annotations/bam/Clusters5/T-cells_sortCoord_rename_RG.bam"
)

fragment_counts <- getCounts(bamfile, 
                             peaks, 
                             paired =  TRUE, 
                             by_rg = TRUE, 
                             format = "bam", 
                             colData = DataFrame(celltype = c("B-cells", "CD", "DCT",
                                                              "Dendritic-cells-11b", "Endothelium", "LOH", "Macrophages",
                                                              "Neutrophils-1", "Neutrophils-2", "NK-cells",
                                                              "Pan-PT", "PCT", "Podocyte", "PST", "T-cells"
                             )))

fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Mmusculus.UCSC.mm10)
filtering_plot <- filterSamplesPlot(fragment_counts, min_depth = 1500, 
                                    min_in_peaks = 0.15, use_plotly = FALSE)
filtering_plot
counts_filtered <- filterSamples(fragment_counts, min_depth = 1500, 
                                 min_in_peaks = 0.15, shiny = FALSE)
counts_filtered <- filterPeaks(counts_filtered, non_overlapping = TRUE)

motifs <- getJasparMotifs()
motif_ix <- matchMotifs(motifs, counts_filtered,
                        genome = BSgenome.Mmusculus.UCSC.mm10)
# computing deviations
dev <- computeDeviations(object = counts_filtered, 
                         annotations = motif_ix)
saveRDS(dev, file = "resutls/mtscATAC/chromVAR-results/kidney-3W-TE-G14102A-4475-80/kidney-3W-TE-G14102A-4475-80_chromVAR.rds")






