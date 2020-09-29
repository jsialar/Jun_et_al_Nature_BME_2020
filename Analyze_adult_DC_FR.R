setwd('/workdir/jsia1/chemical_reprog')
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(magrittr)
library(tidyverse)
library(future)

plan("sequential")
options(future.globals.maxSize = 10 * 1024 ^ 3)
dyn.load('/usr/local/hdf5/lib/libhdf5_hl.so.100')
options(mc.cores=32)

########
ldat_adult <- ReadVelocity(file = "/workdir/jsia1/chemical_reprog/Adult_Fibroblasts_chem_induced_processed/velocyto/Adult_Fibroblasts_chem_induced_processed.loom")
x_adult <- as.Seurat(ldat_adult)
FeatureScatter(x_adult, feature1 = "nCount_spliced", feature2 = "nFeature_spliced")
adult_induced <- subset(x_adult, subset = nFeature_spliced  > 1000 & nCount_spliced   < 50000) %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = all.genes) %>%
  RunPCA(features = VariableFeatures(object = .)) 

adult_induced=  adult_induced %>% FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution=0.5, group.singletons = F) %>% RunTSNE(dims = 1:10)
Idents(adult_induced)
FeaturePlot(adult_induced, features =c('Pax7','Top2a','Myh3','Myog'),reduction='tsne' )

DimPlot(adult_induced, label=T)

adult_induced=SetIdent(adult_induced, value='spliced_snn_res.0.5')
adult_induced=SetIdent(adult_induced, value='clusters_16')
adult_induced=StashIdent(adult_induced, 'clusters_16')

adult_induced_markers <- FindAllMarkers(adult_induced, only.pos = TRUE, 
                                          logfc.threshold = 0.25,
                                          test.use='MAST')
top10_adult_induced=adult_induced_markers %>%
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_logFC)

#DoHeatmap(neonate_non_induced, features = top10_neonate_non_induced$gene) + RotatedAxis() 
DotPlot(SetIdent(adult_induced, value='clusters_16'), features = unique(c(top10_adult_induced$gene,'Pdgfra','Pax7','Myod1')))+ 
  RotatedAxis()
FeaturePlot(adult_induced, features = 'Lyz2')
filter(neonate_induced_markers, gene=='Pax7')
filter(top10_adult_induced, cluster==7)

adult_induced=RenameIdents(adult_induced, 
                           '16'='Macrophages',
                           '15'='Epidermal cells',
                           '14'='Schwann cells',
                           '13'='Myofibroblasts',
                           '12'='Skeletal Muscle cells',
                           '11'='Epidermal cells',
                           '10'='Fibroblasts',
                           '9'='Myofibroblasts',
                           '8'='Melanocytes',
                           '7'='Endothelial cells',
                           '6'='Fibroblasts',
                           '5'='Fibroblasts',
                           '4'='Fibroblasts',
                           '3'='Sebocytes',
                           '2'='Lymph vessel cells',
                           '1'='Fibroblasts',
                           '0'='Fibroblasts')
adult_induced[['Names']]=Idents(adult_induced)
adult_induced=SetIdent(adult_induced, value='Names')
cat_counts=as.vector(table(Idents(adult_induced)))
cat_per=round(cat_counts/sum(cat_counts)*100, 1)
new_names=paste0(levels(Idents(adult_induced)),' ~ ',cat_per,'%','')
names(new_names)=levels(Idents(adult_induced))
adult_induced=RenameIdents(adult_induced, new_names)
ordered_levels=c("Epidermal cells ~ 2.7%",
                 "Fibroblasts ~ 60.9%",
                 "Macrophages ~ 0.5%",
                 "Melanocytes ~ 3.9%",
                 "Myofibroblasts ~ 4.7%",
                 "Schwann cells ~ 1.5%",
                 "Skeletal Muscle cells ~ 1.9%" ,
                 "Lymph vessel cells ~ 11.6%",
                 "Endothelial cells ~ 4.6%",
                 "Sebocytes ~ 7.8%")
Idents(adult_induced)=factor(Idents(adult_induced), levels=ordered_levels,ordered = T)
adult_induced=StashIdent(adult_induced,'Names_per')
adult_induced=SetIdent(adult_induced, value='Names_per')
DimPlot(SetIdent(adult_induced, value='Names_per'), reduction='tsne',cols='Paired', label=F)

adult_induced=SetIdent(adult_induced, value='Names')
ordered_levels_short=gsub(' ~.*','',ordered_levels)
Idents(adult_induced)=factor(Idents(adult_induced), levels=ordered_levels_short,ordered = T)
DimPlot(adult_induced, reduction='tsne',cols='Paired', label=F)




adult_induced=SetIdent(adult_induced, value='Names')

adult_induced_names_markers_coarse <- FindAllMarkers(adult_induced, only.pos = TRUE, 
                                                       logfc.threshold = 0.25,
                                                       test.use='MAST')
top10_adult_induced_coarse=adult_induced_names_markers_coarse %>%
  group_by(cluster) %>% 
  top_n(n = 3, wt = avg_logFC)
filter(top10_adult_induced_coarse, cluster=='Skeletal Muscle cells')
filter(top10_neonate_induced_coarse, gene=='Mylk')
top10_neonate_induced_coarse=mutate(top10_neonate_induced_coarse,gene=ifelse(gene=='Slpi','Mylk',gene)) %>%
  mutate(gene=ifelse(gene=='Cdkn1c','Myod1',gene))
grep('Tyrp1', top10_adult_induced_coarse$gene)


top10_adult_induced_coarse$cluster=as.character(top10_adult_induced_coarse$cluster)
top10_adult_induced_coarse=arrange(top10_adult_induced_coarse, cluster)
genelist=as_tibble(top10_adult_induced_coarse$gene) %>%
  mutate(value=ifelse(value=='H19','Myod1',value)) %>%
  mutate(value=ifelse(value=='Cdkn1c','Des',value)) %>%
  mutate(value=ifelse(value=='Tagln','Myl9',value)) %>%
  mutate(value=ifelse(value=='Hp','Scd1',value)) %>%
  add_row(value='Actg2', .before=22)

adult_induced=SetIdent(adult_induced, value='Names')
neonate_induced=SetIdent(neonate_induced, value='Names')
Idents(adult_induced)=factor(Idents(adult_induced), levels=c(sort(levels(Idents(neonate_induced)))[1:7],
                                                             'Lymph vessel cells',
                                                             'Endothelial cells',
                                                             'Sebocytes'))

top10_neonate_non_induced$cluster=as.character(top10_neonate_non_induced$cluster)
top10_neonate_non_induced=arrange(top10_neonate_non_induced, cluster)
genelist_overall=rev(c('Krt17','Krt14','Krt6a',
                   'Pdgfra','Cdh11','Dcn',
                   'Lyz2','Tyrobp','Pf4',
                   'Pmel','Dct','Tyrp1',
                   'Acta2','Myl9','Mylk',
                   'Plp1','Cryab','Postn',
                   'Myod1','Ttn','Pax7',
                   'Nts','Emcn','Fabp4',
                   'Eng','Angpt2','Gja4',
                   'Cd93','Edn1','Igfbp3',
                   'Scd1','Emb','Igfbp4'))

DotPlot(adult_induced, features = genelist_overall)+ RotatedAxis()
DotPlot(adult_induced, features = unique(genelist[[1]]))+ RotatedAxis()+
  coord_flip()
########### Muscle sub##################
adult_induced=SetIdent(adult_induced, value='Names')
muscle_id_induced=WhichCells(adult_induced, ident='Skeletal Muscle cells')


muscle_subset_induced_adult=subset(adult_induced, cells=muscle_id_induced) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>%
  ScaleData(features = all.genes) %>%
  RunPCA(features = VariableFeatures(object = .)) %>% FindNeighbors(dims = 1:30)
muscle_subset_induced_adult=FindClusters(muscle_subset_induced_adult,resolution=0.5)
Idents(muscle_subset_induced_adult)
DimPlot(muscle_subset_induced_adult)
muscle_names_no=levels(Idents(muscle_subset_induced_adult))
muscle_markers_induced_adult=lapply(muscle_names_no, function(x) FindMarkers(muscle_subset_induced_adult, 
                                                                       ident.1=x, 
                                                                       ident.2=muscle_names_no[!muscle_names_no==x],
                                                                       logfc.threshold = 0.25,
                                                                       test.use='MAST'))
muscle_markers_gene_induced=rbindlist(lapply(muscle_markers_induced, function(x) mutate(x, gene=rownames(x)) %>% top_n(n=5, avg_logFC) %>% select(gene)))
muscle_subset_induced %>% DotPlot(features=c('Pax7',unique(muscle_markers_gene_induced[[1]]) ))+
  RotatedAxis()+
  coord_flip()

muscle_markers[[1]][which(rownames(muscle_markers[[1]])=='Pax7'),]
DimPlot(muscle_subset_induced, cols='Accent')
table(Idents(muscle_subset_induced))

