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
ldat <- ReadVelocity(file = "/workdir/jsia1/chemical_reprog/Neo_Fibroblasts_chem_induced_processed/velocyto/Neo_Fibroblasts_chem_induced_processed.loom")
x <- as.Seurat(x = ldat)
neonate_induced <- subset(x, subset = nFeature_spliced  > 1000 & nCount_spliced   < 50000) %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = all.genes) %>%
  RunPCA(features = VariableFeatures(object = .)) 

neonate_induced=  neonate_induced%>% FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution=0.5, group.singletons = F) %>% RunTSNE(dims = 1:10)
Idents(neonate_induced)
FeaturePlot(neonate_induced, features =c('Pax7','Top2a','Myh3','Myog'),reduction='tsne' )
neonate_induced=SetIdent(neonate_induced, cells=lyz2pos, value='Melanocytes')
neonate_induced=subset(neonate_induced, cells=fibroblast_as_muscle, invert=T) #remove misclassified cells

DimPlot(neonate_induced, label=T)
DimPlot(neonate_induced, cells=lyz2pos, label=T)

neonate_induced=StashIdent(neonate_induced, 'clusters_15')
neonate_induced=SetIdent(neonate_induced, value='clusters_15')
DimPlot(neonate_induced, label=T, reduction='tsne')
neonate_induced_markers <- FindAllMarkers(neonate_induced, only.pos = TRUE, 
                                              logfc.threshold = 0.25,
                                              test.use='MAST')
top10_neonate_induced=neonate_induced_markers %>%
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_logFC)

#DoHeatmap(neonate_non_induced, features = top10_neonate_non_induced$gene) + RotatedAxis() 
DotPlot(SetIdent(neonate_induced, value='clusters_15'), features = unique(c(top10_neonate_induced$gene,'Pdgfra','Pax7','Myod1')))+ 
  RotatedAxis()
FeaturePlot(neonate_induced, features = 'Lyz2')
filter(neonate_induced_markers, gene=='Pax7')
filter(top10_neonate_induced, cluster==8)
neonate_induced=SetIdent(neonate_induced, value='clusters_15')
table(Idents(neonate_induced))

neonate_induced=RenameIdents(neonate_induced, 
                             '15'='Melanocytes',
                             '14'='Schwann cells',
                             '13'='Epidermal cells',
                             '12'='Fibroblasts',
                             '11'='Vascular cells',
                             '10'='Skeletal Muscle cells',
                             '9'='Myofibroblasts',
                             '8'='Fibroblasts',
                             '7'='Fibroblasts',
                             '6'='Skeletal Muscle cells',
                             '5'='Skeletal Muscle cells',
                             '4'='Fibroblasts',
                             '3'='Fibroblasts',
                             '2'='Skeletal Muscle cells',
                             '1'='Fibroblasts',
                             '0'='Fibroblasts')
  
Idents(neonate_induced, cells=lyz2pos)='Macrophages'
table(Idents(neonate_induced))
neonate_induced[['Names']]=Idents(neonate_induced)
neonate_induced=SetIdent(neonate_induced, value='Names')
cat_counts=as.vector(table(Idents(neonate_induced)))
cat_per=round(cat_counts/sum(cat_counts)*100, 1)
new_names=paste0(levels(Idents(neonate_induced)),' ~ ',cat_per,'%','')
names(new_names)=levels(Idents(neonate_induced))
neonate_induced=RenameIdents(neonate_induced, new_names)
Idents(neonate_induced)=factor(Idents(neonate_induced), levels=sort(levels(Idents(neonate_induced))),ordered = T)
neonate_induced=StashIdent(neonate_induced,'Names_per')
neonate_induced=SetIdent(neonate_induced, value='Names_per')
DimPlot(neonate_induced, reduction='tsne',cols='Paired', label=F)



neonate_induced=SetIdent(neonate_induced, value='Names')
neonate_induced=RenameIdents(neonate_induced, 'Melanocytes cells'='Melanocytes') %>% StashIdent('Names')
neonate_induced_names_markers_coarse <- FindAllMarkers(neonate_induced, only.pos = TRUE, 
                                                logfc.threshold = 0.25,
                                                test.use='MAST')
top10_neonate_induced_coarse=neonate_induced_names_markers_coarse %>%
  group_by(cluster) %>% 
  top_n(n = 3, wt = avg_logFC)
filter(top10_neonate_induced_coarse, cluster=='Myofibroblasts')
filter(top10_neonate_induced_coarse, gene=='Mylk')
top10_neonate_induced_coarse=mutate(top10_neonate_induced_coarse,gene=ifelse(gene=='Slpi','Mylk',gene)) %>%
  mutate(gene=ifelse(gene=='Cdkn1c','Myod1',gene))
  
neonate_induced=SetIdent(neonate_induced, value='Names')
Idents(neonate_induced)=factor(Idents(neonate_induced), levels=sort(levels(Idents(neonate_induced))),ordered = T)
DotPlot(neonate_induced, features = genelist_overall)+ RotatedAxis()

# top10_neonate_induced_coarse$cluster=as.character(top10_neonate_induced_coarse$cluster)
# top10_neonate_induced_coarse=arrange(top10_neonate_induced_coarse, cluster)
# DotPlot(neonate_induced, features = unique(top10_neonate_induced_coarse$gene))+ RotatedAxis()

############ Sub fibroblasts
fibroblast_id_induced=WhichCells(neonate_induced, ident='Fibroblasts')

fibrobast_subset_induced=subset(neonate_induced, cells=fibroblast_id_induced) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>%
  ScaleData(features = all.genes) %>%
  RunPCA(features = VariableFeatures(object = .)) %>% FindNeighbors(dims = 1:30)
fibrobast_subset_induced=FindClusters(fibrobast_subset_induced,resolution=0.1)
Idents(fibrobast_subset_induced)
fibroblasts_names_no_induced=levels(Idents(fibrobast_subset_induced))
fibroblast_markers_induced=lapply(fibroblasts_names_no_induced, function(x) FindMarkers(fibrobast_subset_induced, 
                                                  ident.1=x, 
                                                  ident.2=fibroblasts_names_no_induced[!fibroblasts_names_no_induced==x],
                                                  logfc.threshold = 0.25,
                                                  test.use='MAST'))
library(data.table)
fibroblast_markers_induced_gene=rbindlist(lapply(fibroblast_markers_induced, function(x) mutate(x, gene=rownames(x)) %>% top_n(n=3, avg_logFC) %>% select(gene)))
fibrobast_subset_induced %>% DotPlot(features=unique(fibroblast_markers_induced_gene[[1]]))+
  RotatedAxis()+
  coord_flip()

DimPlot(fibrobast_subset_induced, cols='Set1')
fibroblast_as_muscle=WhichCells(fibrobast_subset_induced, idents = 4)

fibrobast_subset_induced=subset(fibrobast_subset_induced, cells = fibroblast_as_muscle, invert=T)
DimPlot(fibrobast_subset_induced, cols='Set1')

fibrobast_subset_induced=FindVariableFeatures(fibrobast_subset_induced, selection.method = "vst", nfeatures = 1000) %>%
  ScaleData(features = all.genes) %>%
  RunPCA(features = VariableFeatures(object = .)) %>% FindNeighbors(dims = 1:30)
fibrobast_subset_induced=FindClusters(fibrobast_subset_induced,resolution=0.1)
fibroblasts_names_no_induced=levels(Idents(fibrobast_subset_induced))
fibroblast_markers_induced=lapply(fibroblasts_names_no_induced, function(x) FindMarkers(fibrobast_subset_induced, 
                                                                                        ident.1=x, 
                                                                                        ident.2=fibroblasts_names_no_induced[!fibroblasts_names_no_induced==x],
                                                                                        logfc.threshold = 0.25,
                                                                                        test.use='MAST'))

######### Muscle sub###############
neonate_induced=SetIdent(neonate_induced, value='Names')
muscle_id_induced=WhichCells(neonate_induced, ident='Skeletal Muscle cells')


muscle_subset_induced=subset(neonate_induced, cells=muscle_id_induced) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>%
  ScaleData(features = all.genes) %>%
  RunPCA(features = VariableFeatures(object = .)) %>% FindNeighbors(dims = 1:30)
muscle_subset_induced=FindClusters(muscle_subset_induced,resolution=0.2)
Idents(muscle_subset_induced)
DimPlot(muscle_subset_induced)
muscle_names_no=levels(Idents(muscle_subset_induced))
muscle_markers_induced=lapply(muscle_names_no, function(x) FindMarkers(muscle_subset_induced, 
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

####### RNA velocity #############

fibrobast_subset_vel_induced <- RunVelocity(object=fibrobast_subset_induced, ncores = 32)
ident.colors <- RColorBrewer::brewer.pal(length(levels(fibrobast_subset_vel_induced)), 'Set1')
names(ident.colors) <- levels(fibrobast_subset_vel_induced)
cell.colors <- ident.colors[Idents(fibrobast_subset_vel_induced)]
names(cell.colors) <- colnames(fibrobast_subset_vel_induced)
embedding_results_induced_fib=show.velocity.on.embedding.cor(emb = Embeddings(object = fibrobast_subset_vel_induced, reduction = "umap"), 
                                                 vel = Tool(object = fibrobast_subset_vel_induced, slot = "RunVelocity"), 
                               cell.colors = ac(x = cell.colors, alpha = 0.5), scale='sqrt',n.cores =32,
                               cex = 0.8, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 100,  
                               do.par = FALSE, cell.border.alpha = 0.1)

show.velocity.on.embedding.cor(emb = Embeddings(object = fibrobast_subset_vel_induced, reduction = "umap"), 
                               vel = Tool(object = fibrobast_subset_vel_induced, slot = "RunVelocity"), 
                               cell.colors = ac(x = cell.colors, alpha = 0.2), scale='sqrt', n.cores =32,
                               cex = 0.5, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 50,  
                               do.par = FALSE, cell.border.alpha = 0.1, arrow.scale=2,cc=embedding_results_induced_fib$cc)

muscle_subset_vel_induced <- RunVelocity(object=muscle_subset_induced, ncores = 32)
ident.colors <- RColorBrewer::brewer.pal(length(levels(muscle_subset_induced)), 'Set1')
names(ident.colors) <- levels(muscle_subset_vel_induced)
cell.colors <- ident.colors[Idents(muscle_subset_vel_induced)]
names(cell.colors) <- colnames(muscle_subset_vel_induced)
embedding_results_induced_muscle=show.velocity.on.embedding.cor(emb = Embeddings(object = muscle_subset_vel_induced, reduction = "umap"), 
                                                 vel = Tool(object = muscle_subset_vel_induced, slot = "RunVelocity"), 
                                                 cell.colors = ac(x = cell.colors, alpha = 0.2), n.cores =32,
                                                 cex = 0.5, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 50,  
                                                 do.par = FALSE, cell.border.alpha = 0.1, arrow.scale=2)

show.velocity.on.embedding.cor(emb = Embeddings(object = muscle_subset_vel_induced, reduction = "umap"), 
                               vel = Tool(object = muscle_subset_vel_induced, slot = "RunVelocity"), 
                               cell.colors = ac(x = cell.colors, alpha = 0.2), scale='sqrt', n.cores =32,
                               cex = 0.5, show.grid.flow = T, min.grid.cell.mass = 0.5, grid.n = 50,  
                               do.par = FALSE, cell.border.alpha = 0.1, arrow.scale=2,cc=embedding_results_induced_muscle$cc)
