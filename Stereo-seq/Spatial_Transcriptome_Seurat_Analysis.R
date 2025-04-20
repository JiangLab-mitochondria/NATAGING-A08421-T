####Spatial transcriptome analysis

library(Seurat)

kid1.seuratobj<-readRDS("./sample_TE_Kid_1_bin50_seurat.rds")
kid2.seuratobj<-readRDS("./sample_TE_Kid_2_bin50_seurat.rds")


################kidney1 individual
####collected marker genes for visualization
figure<-SpatialFeaturePlot(kid1.seuratobj,features=c("Pah","Cyp4b1"))
ggsave(plot=figure,"./TE_kid1_bin50_SCT_umap_SpatialFeaturePlot_markers.pdf",width=13,height=6,units='in',dpi=350)

figure<-FeaturePlot(kid1.seuratobj,features=c("Pah","Cyp4b1"),reduction="umap")
ggsave(plot=figure,"./TE_kid1_bin50_SCT_umap_FeaturePlot_markers.pdf",width=11,height=5,units='in',dpi=350)

#resolution:0.4
kid1.seuratobj<-FindClusters(kid1.seuratobj,resolution=0.4)
kid1_Marker_Table_default<-FindAllMarkers(object=kid1.seuratobj,only.pos = TRUE)
write.csv(kid1_Marker_Table_default,"./TE_kid1_bin50_SCT_res0.4_FindAllMarker.csv",row.names=TRUE,quote=FALSE)

figure<-DotPlot(kid1.seuratobj,features=c('Slc5a8','Slc13a3','Atp11a','Slc5a12','Slc7a7','Slc12a1','Umod','Aqp2','Aqp3','Hsd11b2','Fxyd4','Tagln','Acta2','Dcn','Slc12a3','Calb1','Wnk1','Hba-a1','Hbb-bs','Igkc','Igha','Jchain','Nphs2','Podxl'),group.by="SCT_snn_res.0.4",cols=c("lightgrey","red"))+coord_flip()
ggsave(plot=figure,"./TE_kid1_bin50_SCT_umap_res0.4_markers_dotplot_colored.pdf",width=6.5,height=5,units='in',dpi=350)


################kidney2 individual
figure<-SpatialFeaturePlot(kid2.seuratobj,features=c("Cyp2e1","Cyp4b1"))
ggsave(plot=figure,"./TE_kid2_bin50_SCT_umap_SpatialFeaturePlot_markers.pdf",width=13,height=6,units='in',
       dpi=350)
figure<-FeaturePlot(kid2.seuratobj,features=c("Cyp2e1","Cyp4b1"),reduction="umap")
ggsave(plot=figure,"./TE_kid2_bin50_SCT_umap_FeaturePlot_markers.pdf",width=11,height=5,units='in',dpi=350)


kid2.seuratobj<-FindClusters(kid2.seuratobj,resolution=0.4)
kid2_Marker_default<-FindAllMarkers(object=kid2.seuratobj,only.pos = TRUE)
write.csv(kid2_Marker_default,"./TE_kid2_bin50_SCT_res0.4_FindAllMarker.csv",row.names=TRUE,quote=FALSE)
figure<-DimPlot(kid2.seuratobj,reduction='umap',group.by="SCT_snn_res.0.4",label=TRUE)
ggsave(plot=figure,"./TE_kid2_bin50_SCT_umap_res0.4.pdf",width=6,height=5,units='in',dpi=350)


Genes<-c('Cyp4b1','Cndp2','Slc22a6','Gpx3','Il31ra','Slc5a2','Igfbp5','Slc12a1','Aqp2','Clu','Slc12a3','Spp1','Napsa','Cyp7b1','Serpina1d','Mgp','Csmd1','Cd74','Egf','Umod','Slc12a1','S100g','Klk1','Calb1','Hbb-bs','Hba-a1','Hba-a2','Igkc','Igha','Jchain','Cyp11a1','Fdx1','Hsd3b1','Nphs2','Podxl','Magi2')

figure<-DoHeatmap(object=kid2.seuratobj,features=Genes,group.by="SCT_snn_res.0.4",group.colors=c('#33A02C','#FF69B4','#1F78B4','#FF4500','#B2DF8A','#A6CEE3','#A020F0','#FF7F00','#66CDAA','#DB7093','#1E90FF','#FDBF6F'))+scale_fill_gradientn(colors=c("navy","white","firebrick3"))
ggsave(plot=figure,"./TE_kid2_bin50_SCT_umap_res0.4_markers_Doheatmap_Top3genes_colored_v2.pdf",width=8,height=6.5,units='in',dpi=350)

figure<-DotPlot(kid2.seuratobj,features=c('Slc5a8','Slc13a3','Atp11a','Slc5a12','Slc7a7','Slc12a1','Umod','Aqp2','Aqp3','Hsd11b2','Fxyd4','Tagln','Acta2','Dcn','Slc12a3','Calb1','Wnk1','Hba-a1','Hbb-bs','Igkc','Igha','Jchain','Nphs2','Podxl'),group.by="SCT_snn_res.0.4")+coord_flip()
ggsave(plot=figure,"./TE_kid2_bin50_SCT_umap_res0.4_markers_dotplot.pdf",width=6.5,height=5,units='in',dpi=350)

figure<-DotPlot(kid2.seuratobj,features=c('Slc5a8','Slc13a3','Atp11a','Slc5a12','Slc7a7','Slc12a1','Umod','Aqp2','Aqp3','Hsd11b2','Fxyd4','Tagln','Acta2','Dcn','Slc12a3','Calb1','Wnk1','Hba-a1','Hbb-bs','Igkc','Igha','Jchain','Nphs2','Podxl'),group.by="SCT_snn_res.0.4",cols=c("lightgrey","red"))
ggsave(plot=figure,"./TE_kid2_bin50_SCT_umap_res0.4_markers_dotplot_colored.pdf",width=8,height=5,units='in',dpi=350)

figure<-DotPlot(kid2.seuratobj,features=c('Slc5a8','Slc13a3','Atp11a','Slc5a12','Slc7a7','Slc12a1','Umod','Aqp2','Aqp3','Hsd11b2','Fxyd4','Tagln','Acta2','Dcn','Slc12a3','Calb1','Wnk1','Hba-a1','Hbb-bs','Igkc','Igha','Jchain','Nphs2','Podxl'),group.by="SCT_snn_res.0.4",cols=c("lightgrey","red"))+coord_flip()
ggsave(plot=figure,"./TE_kid2_bin50_SCT_umap_res0.4_markers_dotplot_colored.pdf",width=6.5,height=5,units='in',dpi=350)


############################################
#merged
merged.seuratobj<-merge(kid1.seuratobj,kid2.seuratobj)

VariableFeatures(merged.seuratobj)<-unique(VariableFeatures(kid1.seuratobj),VariableFeatures(kid2.seuratobj))

merged.seuratobj<-RunPCA(merged.seuratobj)
merged.seuratobj<-RunUMAP(merged.seuratobj,dims=1:30)
merged.seuratobj<-FindNeighbors(merged.seuratobj,dims=1:30)
merged.seuratobj<-FindClusters(merged.seuratobj,resolution=0.4)
merged.seuratobj<-FindClusters(merged.seuratobj,resolution=0.6)
merged.seuratobj<-FindClusters(merged.seuratobj,resolution=0.8)
merged.seuratobj<-FindClusters(merged.seuratobj,resolution=1)

merged.seuratobj<-PrepSCTFindMarkers(merged.seuratobj,assay='SCT')

Merged_Marker_Table_default<-FindAllMarkers(object=merged.seuratobj,only.pos = TRUE)
write.csv(Merged_Marker_Table_default,"./Merged_bin50_SCT_res10.6_FindAllMarker.csv",row.names=TRUE,quote=FALSE)

#saveRDS(merged.seuratobj,"./Final_Merged_SCT_Bin50_seuratobj.rds")
merged.seuratobj<-RenameIdents(merged.seuratobj,`0`='PST',`1`='PCT',`2`='PST',`3`='LOH',`4`='Early PT',`5`='Stroma',`6`='PC',`7`='LOH',`8`='DCT',`9`='Erythrocyte',`10`='DCT',`11`='NoMarker',`12`='Podo',`13`='Plasma',`14`='Plasma',`15`='Unknown',`16`='Plasma')
merged.seuratobj@meta.data$celltypes<-Idents(merged.seuratobj)
saveRDS(merged.seuratobj,"./Final_Merged_SCT_Bin50_seuratobj.rds")

library(ggplot2)
figure<-DimPlot(merged.seuratobj,reduction='umap',group.by="celltypes",split.by="orig.ident",label=TRUE)
ggsave(plot=figure,"./merged.umap.sct.celltypes.splitBysamples.pdf",width=10,height=5,units='in',dpi=350)
figure<-SpatialDimPlot(merged.seuratobj,group.by="celltypes")
ggsave(plot=figure,"./merged.umap.sct.celltypes.SpatialDimPlot.pdf",width=13,height=6,units='in',dpi=350)

figure<-DimPlot(merged.seuratobj,reduction='umap',group.by="celltypes",label=TRUE)
ggsave(plot=figure,"./merged.umap.sct.celltypes.pdf",width=7,height=6,units='in',dpi=350)


############################
library(ggplot2)
figure<-DimPlot(merged.seuratobj,reduction='umap',group.by="SCT_snn_res.0.8",split.by="orig.ident",label=TRUE)
ggsave(plot=figure,"./merged.umap.sct.res0.8.splitBysamples.pdf",width=10,height=5,units='in',dpi=350)
figure<-DimPlot(merged.seuratobj,reduction='umap',group.by="SCT_snn_res.0.6",split.by="orig.ident",label=TRUE)
ggsave(plot=figure,"./merged.umap.sct.res0.6.splitBysamples.pdf",width=10,height=5,units='in',dpi=350)
figure<-DimPlot(merged.seuratobj,reduction='umap',group.by="SCT_snn_res.0.4",split.by="orig.ident",label=TRUE)
ggsave(plot=figure,"./merged.umap.sct.res0.4.splitBysamples.pdf",width=10,height=5,units='in',dpi=350)


#######################################################
library(Seurat)
library(ggplot2)
merged.seuratobj<-readRDS("./Final_Merged_SCT_Bin50_seuratobj.rds")
merged.seuratobj@meta.data$FinalTypes<-as.character(merged.seuratobj@meta.data$celltypes)
merged.seuratobj@meta.data$FinalTypes[merged.seuratobj@meta.data$FinalTypes=="Unknown"]<-"Urothelial cell"

figure<-DimPlot(merged.seuratobj,reduction='umap',group.by='FinalTypes',split.by='orig.ident',label=TRUE)
ggsave(plot=figure,"./merged.umap.finaltypes.splitBysample.pdf",width=10,height=5,units='in',dpi=350)

merged.seuratobj@meta.data$FinalTypes<-as.factor(merged.seuratobj@meta.data$FinalTypes)
Idents(merged.seuratobj)<-merged.seuratobj@meta.data$FinalTypes
Merged_Marker_Table_default<-FindAllMarkers(object=merged.seuratobj,only.pos = TRUE)
write.csv(Merged_Marker_Table_default,"./merged.umap.finaltypes.marker.table.csv",row.names=TRUE,quote=FALSE)

merged.seuratobj@meta.data$types<-as.character(merged.seuratobj@meta.data$orig.ident)
merged.seuratobj@meta.data$types<-gsub("sample_TE_","",merged.seuratobj@meta.data$types)
merged.seuratobj@meta.data$types<-factor(merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(merged.seuratobj)<-merged.seuratobj@meta.data$types
merged.seuratobj<-PrepSCTFindMarkers(merged.seuratobj,assay='SCT')
#Minimum UMI unchanged. Skipping re-correction.
Marker_Table_default<-FindAllMarkers(object=merged.seuratobj,only.pos = TRUE)
write.csv(Marker_Table_default,"./merged.umap.kid1.VS.kid2.DEGs.table.csv",row.names=TRUE,quote=FALSE)

##################################
merged.seuratobj<-readRDS("./Final_Merged_SCT_Bin50_seuratobj.rds")
merged.seuratobj@meta.data$FinalTypes<-as.character(merged.seuratobj@meta.data$celltypes)
merged.seuratobj@meta.data$FinalTypes[merged.seuratobj@meta.data$FinalTypes=="Unknown"]<-"Urothelial cell"

DCT.merged.seuratobj<-subset(merged.seuratobj,FinalTypes=="DCT")
DCT.merged.seuratobj@meta.data$types<-as.character(DCT.merged.seuratobj@meta.data$orig.ident)
DCT.merged.seuratobj@meta.data$types<-gsub("sample_TE_","",DCT.merged.seuratobj@meta.data$types)
DCT.merged.seuratobj@meta.data$types<-factor(DCT.merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(DCT.merged.seuratobj)<-DCT.merged.seuratobj@meta.data$types
DCT.merged.seuratobj<-PrepSCTFindMarkers(DCT.merged.seuratobj,assay='SCT')
#Found 2 SCT models. Recorrecting SCT counts using minimum median counts: 2755
Merged_Marker_Table_default<-FindAllMarkers(object=DCT.merged.seuratobj,only.pos = TRUE)

figure<-VlnPlot(DCT.merged.seuratobj,features=c("Egr1","Nr4a1","Mt1",'Spp1','Clu','mt-Co1','Ier3','Ckb','Mt2'),ncol=5,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.DCT.clusters.DEGs.VlnPlot.pdf",width=12,height=6,units='in',dpi=350)


EarlyPT.merged.seuratobj<-subset(merged.seuratobj,FinalTypes=="Early PT")
EarlyPT.merged.seuratobj@meta.data$types<-as.character(EarlyPT.merged.seuratobj@meta.data$orig.ident)
EarlyPT.merged.seuratobj@meta.data$types<-gsub("sample_TE_","",EarlyPT.merged.seuratobj@meta.data$types)
EarlyPT.merged.seuratobj@meta.data$types<-factor(EarlyPT.merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(EarlyPT.merged.seuratobj)<-EarlyPT.merged.seuratobj@meta.data$types
EarlyPT.merged.seuratobj<-PrepSCTFindMarkers(EarlyPT.merged.seuratobj,assay='SCT')
EarlyPT_Marker_Table_default<-FindAllMarkers(object=EarlyPT.merged.seuratobj,only.pos = TRUE)
write.csv(EarlyPT_Marker_Table_default,"./merged.umap.EarlyPT.clusters.Kid1.VS.Kid2.marker.table.csv",row.names=TRUE,quote=FALSE)
figure<-VlnPlot(EarlyPT.merged.seuratobj,features=c("Egr1","Clu","Spp1","Mt1"),ncol=2,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.EarlyPT.clusters.DEGs.VlnPlot.pdf",width=5,height=6,units='in',dpi=350)

sub.merged.seuratobj<-subset(merged.seuratobj,FinalTypes=="PCT")
sub.merged.seuratobj@meta.data$types<-as.character(sub.merged.seuratobj@meta.data$orig.ident)
sub.merged.seuratobj@meta.data$types<-gsub("sample_TE_","",sub.merged.seuratobj@meta.data$types)
sub.merged.seuratobj@meta.data$types<-factor(sub.merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(sub.merged.seuratobj)<-sub.merged.seuratobj@meta.data$types
sub.merged.seuratobj<-PrepSCTFindMarkers(sub.merged.seuratobj,assay='SCT')
sub_Marker_Table_default<-FindAllMarkers(object=sub.merged.seuratobj,only.pos = TRUE)
write.csv(sub_Marker_Table_default,"./merged.umap.PCT.clusters.Kid1.VS.Kid2.marker.table.csv",row.names=TRUE,quote=FALSE)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Egr1","Clu","Spp1"),ncol=2,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.PCT.clusters.DEGs.VlnPlot.pdf",width=5,height=6,units='in',dpi=350)

sub.merged.seuratobj<-subset(merged.seuratobj,FinalTypes=="PST")
sub.merged.seuratobj@meta.data$types<-as.character(sub.merged.seuratobj@meta.data$orig.ident)
sub.merged.seuratobj@meta.data$types<-gsub("sample_TE_","",sub.merged.seuratobj@meta.data$types)
sub.merged.seuratobj@meta.data$types<-factor(sub.merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(sub.merged.seuratobj)<-sub.merged.seuratobj@meta.data$types
sub.merged.seuratobj<-PrepSCTFindMarkers(sub.merged.seuratobj,assay='SCT')
sub_Marker_Table_default<-FindAllMarkers(object=sub.merged.seuratobj,only.pos = TRUE)
write.csv(sub_Marker_Table_default,"./merged.umap.PST.clusters.Kid1.VS.Kid2.marker.table.csv",row.names=TRUE,quote=FALSE)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Egr1","Clu","Spp1","Errfi1","Mgp"),ncol=3,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.PST.clusters.DEGs.VlnPlot.pdf",width=7,height=6,units='in',dpi=350)

sub.merged.seuratobj<-subset(merged.seuratobj,FinalTypes%in%c("Early PT","PCT","PST"))
sub.merged.seuratobj@meta.data$types<-as.character(sub.merged.seuratobj@meta.data$orig.ident)
sub.merged.seuratobj@meta.data$types<-gsub("sample_TE_","",sub.merged.seuratobj@meta.data$types)
sub.merged.seuratobj@meta.data$types<-factor(sub.merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(sub.merged.seuratobj)<-sub.merged.seuratobj@meta.data$types
sub.merged.seuratobj<-PrepSCTFindMarkers(sub.merged.seuratobj,assay='SCT')
sub_Marker_Table_default<-FindAllMarkers(object=sub.merged.seuratobj,only.pos = TRUE)
write.csv(sub_Marker_Table_default,"./merged.umap.PT.clusters.Kid1.VS.Kid2.marker.table.csv",row.names=TRUE,quote=FALSE)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Egr1","Clu","Spp1"),ncol=2,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.PT.clusters.DEGs.VlnPlot.pdf",width=5,height=6,units='in',dpi=350)

sub.merged.seuratobj<-subset(merged.seuratobj,FinalTypes=="LOH")
sub.merged.seuratobj@meta.data$types<-as.character(sub.merged.seuratobj@meta.data$orig.ident)
sub.merged.seuratobj@meta.data$types<-gsub("sample_TE_","",sub.merged.seuratobj@meta.data$types)
sub.merged.seuratobj@meta.data$types<-factor(sub.merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(sub.merged.seuratobj)<-sub.merged.seuratobj@meta.data$types
sub.merged.seuratobj<-PrepSCTFindMarkers(sub.merged.seuratobj,assay='SCT')
#Found 2 SCT models. Recorrecting SCT counts using minimum median counts: 2723.5
sub_Marker_Table_default<-FindAllMarkers(object=sub.merged.seuratobj,only.pos = TRUE)
write.csv(sub_Marker_Table_default,"./merged.umap.LOH.clusters.Kid1.VS.Kid2.marker.table.csv",row.names=TRUE,quote=FALSE)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Cyb5a",'Guca2b','Egf','Slc27a2','Egr1','Nr4a1'),ncol=3,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.LOH.clusters.DEGs.VlnPlot.pdf",width=7,height=6,units='in',dpi=350)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Errfi1",'Clu','Btg2','Junb','Spp1','mt-Co1','mt-Nd5','Gdf15','Mgp','Ier3','Nupr1'),ncol=6,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.LOH.clusters.DEGs.VlnPlot2.pdf",width=13,height=6,units='in',dpi=350)


sub.merged.seuratobj<-subset(merged.seuratobj,FinalTypes=="PC")
sub.merged.seuratobj@meta.data$types<-as.character(sub.merged.seuratobj@meta.data$orig.ident)
sub.merged.seuratobj@meta.data$types<-gsub("sample_TE_","",sub.merged.seuratobj@meta.data$types)
sub.merged.seuratobj@meta.data$types<-factor(sub.merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(sub.merged.seuratobj)<-sub.merged.seuratobj@meta.data$types
sub.merged.seuratobj<-PrepSCTFindMarkers(sub.merged.seuratobj,assay='SCT')
sub_Marker_Table_default<-FindAllMarkers(object=sub.merged.seuratobj,only.pos = TRUE)
write.csv(sub_Marker_Table_default,"./merged.umap.PC.clusters.Kid1.VS.Kid2.marker.table.csv",row.names=TRUE,quote=FALSE)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Acsm2",'Aqp2','Miox','Cyp4b1','Cyb5a','Aqp3','Fxyd4','Slc27a2','Cyp2j5','Slc13a3','Txnip','Hsd11b2'),ncol=6,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.PC.clusters.DEGs.VlnPlot.pdf",width=13,height=6,units='in',dpi=350)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Clu",'mt-Nd1','Egr1','mt-Co1','Spp1','Mt1','Ctsd','Ckb','Nr4a1','Gdf15','mt-Nd5','S100a14'),ncol=6,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.PC.clusters.DEGs.VlnPlot2.pdf",width=13,height=6,units='in',dpi=350)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Btg2",'Mt2','Atf5','Ly6a','Wfdc2','Junb','Mgp','S100a6','S100a1','Nupr1','Cyr61','Samd5'),ncol=6,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.PC.clusters.DEGs.VlnPlot3.pdf",width=13,height=6,units='in',dpi=350)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Sprr1a",'Ier3','Lgals3','Lcn2','Krt8','Cd24a','Chka','Kng2','Pigr','Slc12a3'),ncol=5,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.PC.clusters.DEGs.VlnPlot4.pdf",width=11,height=6,units='in',dpi=350)


sub.merged.seuratobj<-subset(merged.seuratobj,FinalTypes=="Erythrocyte")
sub.merged.seuratobj@meta.data$types<-as.character(sub.merged.seuratobj@meta.data$orig.ident)
sub.merged.seuratobj@meta.data$types<-gsub("sample_TE_","",sub.merged.seuratobj@meta.data$types)
sub.merged.seuratobj@meta.data$types<-factor(sub.merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(sub.merged.seuratobj)<-sub.merged.seuratobj@meta.data$types
sub.merged.seuratobj<-PrepSCTFindMarkers(sub.merged.seuratobj,assay='SCT')
sub_Marker_Table_default<-FindAllMarkers(object=sub.merged.seuratobj,only.pos = TRUE)
write.csv(sub_Marker_Table_default,"./merged.umap.Erythrocyte.clusters.Kid1.VS.Kid2.marker.table.csv",row.names=TRUE,quote=FALSE)

sub.merged.seuratobj<-subset(merged.seuratobj,FinalTypes=="Plasma")
sub.merged.seuratobj@meta.data$types<-as.character(sub.merged.seuratobj@meta.data$orig.ident)
sub.merged.seuratobj@meta.data$types<-gsub("sample_TE_","",sub.merged.seuratobj@meta.data$types)
sub.merged.seuratobj@meta.data$types<-factor(sub.merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(sub.merged.seuratobj)<-sub.merged.seuratobj@meta.data$types
sub.merged.seuratobj<-PrepSCTFindMarkers(sub.merged.seuratobj,assay='SCT')
sub_Marker_Table_default<-FindAllMarkers(object=sub.merged.seuratobj,only.pos = TRUE)
write.csv(sub_Marker_Table_default,"./merged.umap.Plasma.clusters.Kid1.VS.Kid2.marker.table.csv",row.names=TRUE,quote=FALSE)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Ighg2c",'Ighg1','Egr1','Nr4a1','Clu','Spp1','Errfi1'),ncol=4,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.Plasma.clusters.DEGs.VlnPlot.pdf",width=9,height=6,units='in',dpi=350)


sub.merged.seuratobj<-subset(merged.seuratobj,FinalTypes=="Podo")
sub.merged.seuratobj@meta.data$types<-as.character(sub.merged.seuratobj@meta.data$orig.ident)
sub.merged.seuratobj@meta.data$types<-gsub("sample_TE_","",sub.merged.seuratobj@meta.data$types)
sub.merged.seuratobj@meta.data$types<-factor(sub.merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(sub.merged.seuratobj)<-sub.merged.seuratobj@meta.data$types
sub.merged.seuratobj<-PrepSCTFindMarkers(sub.merged.seuratobj,assay='SCT')
sub_Marker_Table_default<-FindAllMarkers(object=sub.merged.seuratobj,only.pos = TRUE)
write.csv(sub_Marker_Table_default,"./merged.umap.Podo.clusters.Kid1.VS.Kid2.marker.table.csv",row.names=TRUE,quote=FALSE)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Ren1",'Egr1','Nr4a1','Clu','Spp1','Cyr61'),ncol=3,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.Podo.clusters.DEGs.VlnPlot.pdf",width=7,height=6,units='in',dpi=350)

sub.merged.seuratobj<-subset(merged.seuratobj,FinalTypes=="Stroma")
sub.merged.seuratobj@meta.data$types<-as.character(sub.merged.seuratobj@meta.data$orig.ident)
sub.merged.seuratobj@meta.data$types<-gsub("sample_TE_","",sub.merged.seuratobj@meta.data$types)
sub.merged.seuratobj@meta.data$types<-factor(sub.merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(sub.merged.seuratobj)<-sub.merged.seuratobj@meta.data$types
sub.merged.seuratobj<-PrepSCTFindMarkers(sub.merged.seuratobj,assay='SCT')
sub_Marker_Table_default<-FindAllMarkers(object=sub.merged.seuratobj,only.pos = TRUE)
write.csv(sub_Marker_Table_default,"./merged.umap.Stroma.clusters.Kid1.VS.Kid2.marker.table.csv",row.names=TRUE,quote=FALSE)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Acsm2",'Cyp4b1','Egr1','Nr4a1','Cyr61','Errfi1','Btg2','Thbs1','Bgn','Aebp1','Col4a1','Wfdc2','Sparc','Dcn'),ncol=7,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.Stroma.clusters.DEGs.VlnPlot.pdf",width=14,height=6,units='in',dpi=350)

sub.merged.seuratobj<-subset(merged.seuratobj,FinalTypes=="Urothelial cell")
sub.merged.seuratobj@meta.data$types<-as.character(sub.merged.seuratobj@meta.data$orig.ident)
sub.merged.seuratobj@meta.data$types<-gsub("sample_TE_","",sub.merged.seuratobj@meta.data$types)
sub.merged.seuratobj@meta.data$types<-factor(sub.merged.seuratobj@meta.data$types,levels=c("Kid_1","Kid_2"))
Idents(sub.merged.seuratobj)<-sub.merged.seuratobj@meta.data$types
sub.merged.seuratobj<-PrepSCTFindMarkers(sub.merged.seuratobj,assay='SCT')
sub_Marker_Table_default<-FindAllMarkers(object=sub.merged.seuratobj,only.pos = TRUE)
write.csv(sub_Marker_Table_default,"./merged.umap.Urothelialcells.clusters.Kid1.VS.Kid2.marker.table.csv",row.names=TRUE,quote=FALSE)
figure<-VlnPlot(sub.merged.seuratobj,features=c("Ighg2c",'8430408G22Rik','Gsdmc3','Wfdc15b','Fam25c','Ndufb3','Upk1a','Cryab','Egr1','Nr4a1','Junb','Actb','Errfi1','Clu','Hbegf','Btg2','Mt1','Dnajb1','Cyp7b1','Acsm3'),ncol=10,group.by='types',pt.size=0)
ggsave(plot=figure,"./merged.umap.Urothelialcells.clusters.DEGs.VlnPlot.pdf",width=18,height=6,units='in',dpi=350)

##########################Final colored figure
merged.seuratobj<-readRDS("./Final_Merged_SCT_Bin50_seuratobj.rds")
merged.seuratobj@meta.data$FinalTypes<-as.character(merged.seuratobj@meta.data$celltypes)
merged.seuratobj@meta.data$FinalTypes[merged.seuratobj@meta.data$FinalTypes=="Unknown"]<-"Urothelial cell"

figure2<-SpatialDimPlot(merged.seuratobj,group.by="FinalTypes",label=TRUE,label.size=2)
ggsave(plot=figure2,"./merged.umap.finaltypes.slide12.pdf",width=10,height=5,units='in',dpi=350)
figure2<-SpatialDimPlot(merged.seuratobj,group.by="FinalTypes",label=FALSE,label.size=2)
ggsave(plot=figure2,"./merged.umap.finaltypes.slide12.v2.pdf",width=10,height=5,units='in',dpi=350)
library(RColorBrewer)
brewer.pal(12,"Paired")
#[1] "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F"
# [8] "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"

figure<-DimPlot(merged.seuratobj,group.by="FinalTypes",reduction='umap',label=TRUE)
ggsave(plot=figure,"./merged.umap.finaltypes.v1.pdf",width=6,height=4.5,units='in',dpi=350)

figure<-DimPlot(merged.seuratobj,group.by="FinalTypes",reduction='umap',label=TRUE,split.by="orig.ident")
ggsave(plot=figure,"./merged.umap.finaltypes.v2.pdf",width=10,height=4.5,units='in',dpi=350)

figure<-DimPlot(merged.seuratobj,group.by="FinalTypes",reduction='umap',label=TRUE,cols=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#FF7F00",'Podo'="#CAB2D6",'PST'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"))
figure<-DimPlot(merged.seuratobj,group.by="FinalTypes",reduction='umap',label=TRUE,cols=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#CAB2D6",'Podo'="#FF7F00",'PST'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"))
ggsave(plot=figure,"./merged.umap.finaltypes.coloredPaired.pdf",width=6,height=4.5,units='in',dpi=350)

figure<-DimPlot(merged.seuratobj,group.by="FinalTypes",reduction='umap',label=TRUE,cols=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#FF7F00",'Podo'="#CAB2D6",'PST'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"),split.by="orig.ident")
figure<-DimPlot(merged.seuratobj,group.by="FinalTypes",reduction='umap',label=TRUE,cols=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#CAB2D6",'Podo'="#FF7F00",'PST'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"),split.by="orig.ident")
ggsave(plot=figure,"./merged.umap.finaltypes.coloredPaired.v2.pdf",width=10,height=4.5,units='in',dpi=350)

figure<-SpatialDimPlot(merged.seuratobj,group.by="FinalTypes",label=FALSE,label.size=2,cols=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#FF7F00",'Podo'="#CAB2D6",'PST'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"))
figure<-SpatialDimPlot(merged.seuratobj,group.by="FinalTypes",label=FALSE,label.size=2,cols=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#CAB2D6",'Podo'="#FF7F00",'PST'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"))
ggsave(plot=figure,"./merged.umap.finaltypes.spatialslide.pdf",width=10,height=5,units='in',dpi=350)

figure<-SpatialDimPlot(merged.seuratobj,group.by="FinalTypes",label=TRUE,label.size=2,cols=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#FF7F00",'Podo'="#CAB2D6",'PST'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"))
figure<-SpatialDimPlot(merged.seuratobj,group.by="FinalTypes",label=TRUE,label.size=2,cols=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#CAB2D6",'Podo'="#FF7F00",'PST'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"))
ggsave(plot=figure,"./merged.umap.finaltypes.spatialslide2.pdf",width=10,height=5,units='in',dpi=350)

Idents(merged.seuratobj)<-merged.seuratobj@meta.data$FinalTypes
Genes<-c('Nphs2','Nphs1','Podxl','Tagln', 'Acta2', 'Myh11','Sult1d1', 'Fut9', 'Slc34a1', 'Keg1','Slc5a12', 'Slc7a7', 'Slc5a2','Atp11a', 'Slc5a8', 'Slc13a3','Psca', 'Krt8', 'Upk2','Slc12a1', 'Umod','Slc12a3','Calb1','Wnk1','Pvalb','Wnk4','Hsd11b2', 'Aqp2', 'Fxyd4', 'Aqp3','Hba-a1', 'Hbb-bs', 'Hba-a2', 'Hbb-bt','Igkc','Igha','Jchain')

set.seed(2024)
subseuratobj<-subset(merged.seuratobj,downsample=300)
figure<-DoHeatmap(object=subseuratobj,features=Genes,group.by="FinalTypes")+scale_fill_gradientn(colors=c("navy","white","firebrick3"))
figure<-DoHeatmap(object=subseuratobj,features=Genes,group.by="FinalTypes")+scale_fill_gradientn(colors=c("navy","white","firebrick3"))
ggsave(plot=figure,"./merged.umap.finaltypes.markers.Doheatmap.downsample200.pdf",width=6,height=6.5,units='in',dpi=350)
figure<-DoHeatmap(object=subseuratobj,features=Genes,group.by="FinalTypes",group.colors=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#CAB2D6",'Podo'="#FF7F00",'PST'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"))+scale_fill_gradientn(colors=c("navy","white","firebrick3"))
ggsave(plot=figure,"./merged.umap.finaltypes.markers.Doheatmap.downsample200.coloredPaired.pdf",width=6,height=6.5,units='in',dpi=350)

figure<-DimPlot(merged.seuratobj,group.by="FinalTypes",reduction='umap',label=TRUE,cols=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#CAB2D6",'Podo'="#FF7F00",'PST'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"),split.by="orig.ident")
ggsave(plot=figure,"./merged.umap.finaltypes.coloredPaired.v2.pdf",width=10,height=4.5,units='in',dpi=350)


figure<-DoHeatmap(object=merged.seuratobj,features=Genes,group.by="FinalTypes")+scale_fill_gradientn(colors=c("navy","white","firebrick3"))
ggsave(plot=figure,"./merged.umap.finaltypes.markers.Doheatmap.pdf",width=7,height=6.5,units='in',dpi=350)

figure<-DoHeatmap(object=merged.seuratobj,features=Genes,group.by="FinalTypes",,group.colors=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#CAB2D6",'Podo'="#FF7F00",'PST'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"))+scale_fill_gradientn(colors=c("navy","white","firebrick3"))
ggsave(plot=figure,"./merged.umap.finaltypes.markers.Doheatmap.coloredPaired.pdf",width=7,height=6.5,units='in',dpi=350)

figure<-DotPlot(merged.seuratobj,features=Genes,group.by="FinalTypes",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
ggsave(plot=figure,"./merged.umap.finaltypes.dotplot.pdf",width=10,height=4,units='in',dpi=350)

Genes<-c('Nphs2','Podxl','Tagln', 'Acta2', 'Dcn','Sult1d1', 'Fut9', 'Slc34a1', 'Keg1','Slc5a12', 'Slc7a7', 'Slc5a2','Atp11a', 'Slc5a8', 'Slc13a3','Psca', 'Krt8', 'Upk2','Slc12a1', 'Umod','Slc12a3','Calb1','Wnk1','Pvalb','Wnk4','Hsd11b2', 'Aqp2', 'Fxyd4', 'Aqp3','Hba-a1', 'Hbb-bs', 'Hba-a2', 'Hbb-bt','Igkc','Igha','Jchain')
figure<-DotPlot(merged.seuratobj,features=Genes,group.by="FinalTypes",cols=c("lightgrey","red"))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
ggsave(plot=figure,"./merged.umap.finaltypes.dotplot.v2.pdf",width=10,height=4,units='in',dpi=350)

merged.seuratobj@meta.data$FinalTypes<-as.character(merged.seuratobj@meta.data$FinalTypes)
merged.seuratobj@meta.data$FinalTypes<-factor(merged.seuratobj@meta.data$FinalTypes,levels=c("PST","LOH","PCT","DCT","Early PT","Stroma","PC","Plasma","Erythrocyte","NoMarker","Podo","Urothelial cell"))

figure<-DoHeatmap(object=merged.seuratobj,features=Genes,group.by="FinalTypes",group.colors=c('#33A02C','#FF69B4','#1F78B4','#FF4500','#B2DF8A','#A6CEE3','#A020F0','#FF7F00','#66CDAA','#DB7093','#1E90FF','#FDBF6F'))+scale_fill_gradientn(colors=c("navy","white","firebrick3"))
ggsave(plot=figure,"./merged.umap.finaltypes.markers.Doheatmap.colored.auto.pdf",width=8,height=6.5,units='in',dpi=350)
figure<-DoHeatmap(object=merged.seuratobj,features=Genes,group.by="FinalTypes",,group.colors=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#CAB2D6",'Podo'="#FF7F00",'PST'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"))+scale_fill_gradientn(colors=c("navy","white","firebrick3"))
ggsave(plot=figure,"./merged.umap.finaltypes.markers.Doheatmap.colored.auto.v2.pdf",width=8,height=6.5,units='in',dpi=350)


merged.seuratobj@meta.data$Group<-as.character(merged.seuratobj@meta.data$orig.ident)
merged.seuratobj@meta.data$Group<-gsub("sample_TE_Kid_1","Ctrl",merged.seuratobj@meta.data$Group)
merged.seuratobj@meta.data$Group<-gsub("sample_TE_Kid_2","TrnE",merged.seuratobj@meta.data$Group)
figure<-DotPlot(merged.seuratobj,features=c('Egr1','Spp1','Clu','Nr4a1','Mt1','Ckb'),group.by="Group",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
tmp_data<-figure$data
input_data<-data.frame(row.names=unique(tmp_data$features.plot),Ctrl=tmp_data$avg.exp[tmp_data$id=="Ctrl"],TrnE=tmp_data$avg.exp[tmp_data$id=="TrnE"])
library(pheatmap)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,scale='none',filename="./merged.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=3)
annotated_col<-data.frame(row.names=c('Ctrl','TrnE'),Group=c('Ctrl','TrnE'))
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',filename="./merged.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=3)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='row',filename="./merged.umap.kid1.VS.kid2.heatmap.sclaedByrow.pdf",width=3.5,height=3)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',cellwidth=25,cellheight=25,filename="./merged.umap.kid1.VS.kid2.heatmap.v2.pdf",width=3.5,height=3)


Genes_Data<-read.csv("./merged_kid1_VS_kid2_heatmap_Genes.txt",sep='\t',header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)

tmp_seuratobj<-subset(merged.seuratobj,FinalTypes=='PST')
figure<-DotPlot(tmp_seuratobj,features=Genes_Data$Genes[Genes_Data$Celltypes=='PST'],group.by="Group",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
tmp_data<-figure$data
input_data<-data.frame(row.names=unique(tmp_data$features.plot),Ctrl=tmp_data$avg.exp[tmp_data$id=="Ctrl"],TrnE=tmp_data$avg.exp[tmp_data$id=="TrnE"])
annotated_col<-data.frame(row.names=c('Ctrl','TrnE'),Group=c('Ctrl','TrnE'))
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',filename="./merged.PST.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=3)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',cellwidth=15,cellheight=15,filename="./merged.PST.umap.kid1.VS.kid2.heatmap.v2.pdf",width=3.5,height=3)

tmp_seuratobj<-subset(merged.seuratobj,FinalTypes=='LOH')
figure<-DotPlot(tmp_seuratobj,features=Genes_Data$Genes[Genes_Data$Celltypes=='LOH'],group.by="Group",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
tmp_data<-figure$data
input_data<-data.frame(row.names=unique(tmp_data$features.plot),Ctrl=tmp_data$avg.exp[tmp_data$id=="Ctrl"],TrnE=tmp_data$avg.exp[tmp_data$id=="TrnE"])
annotated_col<-data.frame(row.names=c('Ctrl','TrnE'),Group=c('Ctrl','TrnE'))
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',filename="./merged.LOH.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=6)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',cellwidth=15,cellheight=15,filename="./merged.LOH.umap.kid1.VS.kid2.heatmap.v2.pdf",width=3.5,height=6)

tmp_seuratobj<-subset(merged.seuratobj,FinalTypes=='PCT')
figure<-DotPlot(tmp_seuratobj,features=Genes_Data$Genes[Genes_Data$Celltypes=='PCT'],group.by="Group",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
tmp_data<-figure$data
input_data<-data.frame(row.names=unique(tmp_data$features.plot),Ctrl=tmp_data$avg.exp[tmp_data$id=="Ctrl"],TrnE=tmp_data$avg.exp[tmp_data$id=="TrnE"])
annotated_col<-data.frame(row.names=c('Ctrl','TrnE'),Group=c('Ctrl','TrnE'))
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',filename="./merged.PCT.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=2.5)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',cellwidth=15,cellheight=15,filename="./merged.PCT.umap.kid1.VS.kid2.heatmap.v2.pdf",width=3.5,height=3)

tmp_seuratobj<-subset(merged.seuratobj,FinalTypes=='DCT')
figure<-DotPlot(tmp_seuratobj,features=Genes_Data$Genes[Genes_Data$Celltypes=='DCT'],group.by="Group",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
tmp_data<-figure$data
input_data<-data.frame(row.names=unique(tmp_data$features.plot),Ctrl=tmp_data$avg.exp[tmp_data$id=="Ctrl"],TrnE=tmp_data$avg.exp[tmp_data$id=="TrnE"])
annotated_col<-data.frame(row.names=c('Ctrl','TrnE'),Group=c('Ctrl','TrnE'))
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',filename="./merged.DCT.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=3.5)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',cellwidth=15,cellheight=15,filename="./merged.DCT.umap.kid1.VS.kid2.heatmap.v2.pdf",width=3.5,height=3.5)

tmp_seuratobj<-subset(merged.seuratobj,FinalTypes=='Early PT')
figure<-DotPlot(tmp_seuratobj,features=Genes_Data$Genes[Genes_Data$Celltypes=='Early PT'],group.by="Group",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
tmp_data<-figure$data
input_data<-data.frame(row.names=unique(tmp_data$features.plot),Ctrl=tmp_data$avg.exp[tmp_data$id=="Ctrl"]
                       ,TrnE=tmp_data$avg.exp[tmp_data$id=="TrnE"])
annotated_col<-data.frame(row.names=c('Ctrl','TrnE'),Group=c('Ctrl','TrnE'))
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',filename="./merged.EarlyPT.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=2.5)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',cellwidth=15,cellheight=15,filename="./merged.EarlyPT.umap.kid1.VS.kid2.heatmap.v2.pdf",width=3.5,height=3)

tmp_seuratobj<-subset(merged.seuratobj,FinalTypes=='Stroma')
figure<-DotPlot(tmp_seuratobj,features=Genes_Data$Genes[Genes_Data$Celltypes=='Stroma'],group.by="Group",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
tmp_data<-figure$data
input_data<-data.frame(row.names=unique(tmp_data$features.plot),Ctrl=tmp_data$avg.exp[tmp_data$id=="Ctrl"]
                       ,TrnE=tmp_data$avg.exp[tmp_data$id=="TrnE"])
annotated_col<-data.frame(row.names=c('Ctrl','TrnE'),Group=c('Ctrl','TrnE'))
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',filename="./merged.Stroma.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=5)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',cellwidth=15,cellheight=15,filename="./merged.Stroma.umap.kid1.VS.kid2.heatmap.v2.pdf",width=3.5,height=5)

tmp_seuratobj<-subset(merged.seuratobj,FinalTypes=='PC')
figure<-DotPlot(tmp_seuratobj,features=Genes_Data$Genes[Genes_Data$Celltypes=='PC'],group.by="Group",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
tmp_data<-figure$data
input_data<-data.frame(row.names=unique(tmp_data$features.plot),Ctrl=tmp_data$avg.exp[tmp_data$id=="Ctrl"]
                       ,TrnE=tmp_data$avg.exp[tmp_data$id=="TrnE"])
annotated_col<-data.frame(row.names=c('Ctrl','TrnE'),Group=c('Ctrl','TrnE'))
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',filename="./merged.PC.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=8)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',cellwidth=15,cellheight=15,filename="./merged.PC.umap.kid1.VS.kid2.heatmap.v2.pdf",width=3.5,height=10)

tmp_seuratobj<-subset(merged.seuratobj,FinalTypes=='Plasma')
figure<-DotPlot(tmp_seuratobj,features=Genes_Data$Genes[Genes_Data$Celltypes=='Plasma'],group.by="Group",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
tmp_data<-figure$data
input_data<-data.frame(row.names=unique(tmp_data$features.plot),Ctrl=tmp_data$avg.exp[tmp_data$id=="Ctrl"]
                       ,TrnE=tmp_data$avg.exp[tmp_data$id=="TrnE"])
annotated_col<-data.frame(row.names=c('Ctrl','TrnE'),Group=c('Ctrl','TrnE'))
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',filename="./merged.Plasma.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=3)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',cellwidth=15,cellheight=15,filename="./merged.Plasma.umap.kid1.VS.kid2.heatmap.v2.pdf",width=3.5,height=3)

tmp_seuratobj<-subset(merged.seuratobj,FinalTypes=='Erythrocyte')
figure<-DotPlot(tmp_seuratobj,features=Genes_Data$Genes[Genes_Data$Celltypes=='Erythrocyte'],group.by="Group",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
tmp_data<-figure$data
input_data<-data.frame(row.names=unique(tmp_data$features.plot),Ctrl=tmp_data$avg.exp[tmp_data$id=="Ctrl"]
                       ,TrnE=tmp_data$avg.exp[tmp_data$id=="TrnE"])
annotated_col<-data.frame(row.names=c('Ctrl','TrnE'),Group=c('Ctrl','TrnE'))
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',filename="./merged.Erythrocyte.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=5)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',cellwidth=15,cellheight=15,filename="./merged.Erythrocyte.umap.kid1.VS.kid2.heatmap.v2.pdf",width=3.5,height=8)

tmp_seuratobj<-subset(merged.seuratobj,FinalTypes=='Podo')
figure<-DotPlot(tmp_seuratobj,features=Genes_Data$Genes[Genes_Data$Celltypes=='Podo'],group.by="Group",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
tmp_data<-figure$data
input_data<-data.frame(row.names=unique(tmp_data$features.plot),Ctrl=tmp_data$avg.exp[tmp_data$id=="Ctrl"]
                       ,TrnE=tmp_data$avg.exp[tmp_data$id=="TrnE"])
annotated_col<-data.frame(row.names=c('Ctrl','TrnE'),Group=c('Ctrl','TrnE'))
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',filename="./merged.Podo.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=3)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',cellwidth=15,cellheight=15,filename="./merged.Podo.umap.kid1.VS.kid2.heatmap.v2.pdf",width=3.5,height=3)

tmp_seuratobj<-subset(merged.seuratobj,FinalTypes=='Urothelial cell')
figure<-DotPlot(tmp_seuratobj,features=Genes_Data$Genes[Genes_Data$Celltypes=='Urothelial cel'],group.by="Group",cols=c("lightgrey",'red'))+theme(axis.text.x =element_text(face='bold',size=10,angle=90,hjust=0.5,vjust=0.5))
tmp_data<-figure$data
input_data<-data.frame(row.names=unique(tmp_data$features.plot),Ctrl=tmp_data$avg.exp[tmp_data$id=="Ctrl"]
                       ,TrnE=tmp_data$avg.exp[tmp_data$id=="TrnE"])
annotated_col<-data.frame(row.names=c('Ctrl','TrnE'),Group=c('Ctrl','TrnE'))
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',filename="./merged.Urothelialcell.umap.kid1.VS.kid2.heatmap.pdf",width=3.5,height=4)
pheatmap(input_data,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=annotated_col,scale='none',cellwidth=15,cellheight=15,filename="./merged.Urothelialcell.umap.kid1.VS.kid2.heatmap.v2.pdf",width=3.5,height=4)

#############################################################################
##############20240313
merged.seuratobj@meta.data$FinalTypes2<-merged.seuratobj@meta.data$FinalTypes
idx<-rownames(merged.seuratobj@meta.data[merged.seuratobj@meta.data$SCT_snn_res.0.6==0,])
merged.seuratobj@meta.data$FinalTypes2[match(idx,rownames(merged.seuratobj@meta.data))]<-"PST1"
idx<-rownames(merged.seuratobj@meta.data[merged.seuratobj@meta.data$SCT_snn_res.0.6==2,])
merged.seuratobj@meta.data$FinalTypes2[match(idx,rownames(merged.seuratobj@meta.data))]<-"PST2"
Idents(merged.seuratobj)<-merged.seuratobj@meta.data$FinalTypes2
Merged_Marker_Table_default2<-FindAllMarkers(object=merged.seuratobj,only.pos = TRUE)
write.csv(Merged_Marker_Table_default2,"./merged.umap.finaltypes2.marker.table.csv",row.names=TRUE,quote=FALSE)

figure<-DimPlot(merged.seuratobj,group.by="FinalTypes2",reduction='umap',label=TRUE,cols=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#CAB2D6",'Podo'="#FF7F00",'PST1'='#b99cd8','PST2'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"))
ggsave(plot=figure,"./merged.umap.finaltypes2.coloredPaired.pdf",width=6,height=4.5,units='in',dpi=350)
figure<-DimPlot(merged.seuratobj,group.by="FinalTypes2",reduction='umap',label=TRUE,cols=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#CAB2D6",'Podo'="#FF7F00",'PST1'='#b99cd8','PST2'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"),split.by="orig.ident")
ggsave(plot=figure,"./merged.umap.finaltypes2.coloredPaired.splitbysamples.pdf",width=10,height=4.5,units='in',dpi=350)


library(dplyr)
TopGenes<-Merged_Marker_Table_default2%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)

merged.seuratobj@meta.data$FinalTypes2<-factor(merged.seuratobj@meta.data$FinalTypes2,levels=names(rev(sort(table(merged.seuratobj@meta.data$FinalTypes2)))))
figure<-DoHeatmap(object=merged.seuratobj,features=TopGenes$gene,group.by="FinalTypes2",group.colors=c("DCT"="#A6CEE3",'Early PT'="#1F78B4",'Erythrocyte'="#B2DF8A",'LOH'="#E31A1C",'NoMarker'="#B15928",'PC'="#33A02C",'PCT'="#FDBF6F",'Plasma'="#CAB2D6",'Podo'="#FF7F00",'PST1'='#b99cd8','PST2'="#6A3D9A",'Stroma'="#FFFF99",'Urothelial cell'="#FB9A99"))+scale_fill_gradientn(colors=c("navy","white","firebrick3"))
ggsave(plot=figure,"./merged.umap.finaltypes2.marker.Top5.Doheatmap.pdf",width=9,height=8.5,units='in',dpi=350)

##############################################################
##########20241204
library(Seurat)
merged.seuratobj<-readRDS("./Final_Merged_SCT_Bin50_seuratobj.rds")
merged.seuratobj@meta.data$FinalTypes<-as.character(merged.seuratobj@meta.data$celltypes)
merged.seuratobj@meta.data$FinalTypes[merged.seuratobj@meta.data$FinalTypes=="Unknown"]<-"Urothelial cell"

Idents(merged.seuratobj)<-merged.seuratobj@meta.data$FinalTypes
marker_table<-FindAllMarkers(object=merged.seuratobj,only.pos = TRUE)

Kid1<-readRDS("./sample_TE_Kid_1_bin50_seurat.rds")
idx<-intersect(rownames(merged.seuratobj@meta.data),rownames(Kid1@meta.data))
Kid1@meta.data$Finaltype<-merged.seuratobj@meta.data$FinalTypes[match(rownames(Kid1@meta.data),rownames(merged.seuratobj@meta.data))]
Idents(Kid1)<-Kid1@meta.data$Finaltype
kid1_marker<-FindAllMarkers(object=Kid1,only.pos = TRUE)
write.csv(kid1_marker,"./TE_Kid_1_bin50_seurat_Finaltypes_marker_table_onlypositive.csv",row.names=TRUE,quote=FALSE)

Kid2<-readRDS("./sample_TE_Kid_2_bin50_seurat.rds")
idx<-intersect(rownames(merged.seuratobj@meta.data),rownames(Kid2@meta.data))
Kid2@meta.data$Finaltype<-merged.seuratobj@meta.data$FinalTypes[match(rownames(Kid2@meta.data),rownames(merged.seuratobj@meta.data))]
Idents(Kid2)<-Kid2@meta.data$Finaltype
kid2_marker<-FindAllMarkers(object=Kid2,only.pos = TRUE)
write.csv(kid1_marker,"./TE_Kid_2_bin50_seurat_Finaltypes_marker_table_onlypositive.csv",row.names=TRUE,quote=FALSE)


