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
x_stromal <- as.Seurat(ldat_stromal)
dim(x_stromal)

FeatureScatter(x_stromal, feature1 = "nCount_spliced", feature2 = "nFeature_spliced")
stromal <- subset(x_stromal, subset = nFeature_spliced  > 1000 & nCount_spliced   < 75000) %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = all.genes) %>%
  RunPCA(features = VariableFeatures(object = .)) 

stromal=  stromal %>% FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution=0.5, group.singletons = F) %>%
  RunTSNE(dims = 1:10)

stromal=RunTSNE(stromal, dims = 1:20)

SetIdent(adult_induced, value='Names_test')
Idents(stromal)
FeaturePlot(stromal, features =c('Pax7','Top2a','Myh3','Myog'),reduction='tsne' )

DimPlot(stromal, label=T)

stromal=StashIdent(stromal, 'clusters_16')

stromal=SetIdent(stromal, value='clusters_16')

stromal_markers <- FindAllMarkers(stromal, only.pos = TRUE, 
                                        logfc.threshold = 0.25,
                                        test.use='MAST')
top10_stromal=stromal_markers %>%
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_logFC)
cluster_6_genes=stromal_markers %>% filter(cluster=='6') %>% top_n(20, wt=avg_logFC)
#DoHeatmap(neonate_non_induced, features = top10_neonate_non_induced$gene) + RotatedAxis() 
DotPlot(SetIdent(stromal, value='clusters_16'), features = c('Pdgfra','Pecam1','Cd19','Cd28','Tnmd','Myl9','Myod1','Cd68','Ckm'))+ 
  RotatedAxis()
DotPlot(SetIdent(stromal, value='clusters_16'), features = unique(c('Pax7', 'Pdgfra',top10_stromal$gene)))+ 
  RotatedAxis()
DotPlot(SetIdent(stromal, value='clusters_16'), features = unique(c('Pdgfra',cluster_6_genes$gene)))+ 
  RotatedAxis()

FeaturePlot(stromal, features = 'Lyz2')
filter(neonate_induced_markers, gene=='Pax7')
filter(top10_adult_induced, cluster==7)
stromal=SetIdent(stromal, value='clusters_16')
stromal=RenameIdents(stromal,'16'='Smooth muscle cells',
                     '15'='Chondrocytes',
                     '14'='Schwann cells',
                     '13'='Fibroblasts',
                     '12'='Macrophages',
                     '11'='Myofibroblasts',
                     '10'='Tenocytes',
                     '9'='Skeletal muscle cells' ,
                     '8'='Macrophages',
                     '7'='Chondrocytes',
                     '6'='Schwann cells',
                     '5'='Skeletal muscle cells',
                     '4'='Endothelial',
                     '3'='Tenocytes',
                     '2'='Fibroblasts',
                     '1'='Endothelial',
                     '0'='Fibroblasts')
stromal[['Names']]=Idents(stromal)
stromal=SetIdent(stromal, value='Names')
cat_counts=as.vector(table(Idents(stromal)))
cat_per=round(cat_counts/sum(cat_counts)*100, 1)
new_names=paste0(levels(Idents(stromal)),' ~ ',cat_per,'%','')
names(new_names)=levels(Idents(stromal))
stromal=RenameIdents(stromal, new_names)
stromal[['Names_per']]=Idents(stromal)

DimPlot(SetIdent(stromal, value='Names_test'), reduction='tsne',label=T)
DimPlot(SetIdent(stromal, value='Names_per'), reduction='tsne',cols='Paired', label=F)



stromal=SetIdent(stromal, value='Names')
stromal_names_markers_coarse <- FindAllMarkers(stromal, only.pos = TRUE, 
                                                     logfc.threshold = 0.25,
                                                     test.use='MAST')
top10_stromal_coarse=stromal_names_markers_coarse %>%
  group_by(cluster) %>% 
  top_n(n = 3, wt = avg_logFC)
filter(top10_adult_induced_coarse, cluster=='Skeletal Muscle cells')
filter(top10_neonate_induced_coarse, gene=='Mylk')

genelist=as_tibble(top10_stromal_coarse$gene) %>%
  mutate(value=ifelse(value=='Gdf15','Angpt2',value)) %>%
  mutate(value=ifelse(value=='Ednrb','Eng',value)) %>%
  mutate(value=ifelse(value=='Col1a2','Tnmd',value)) %>%
  mutate(value=ifelse(value=='Gsn','Pdgfra',value)) %>%
  mutate(value=ifelse(value=='Lgals1','Pax7',value)) %>%
  filter(!value %in% c('Chrna1','Ankrd1','Gal'))
  


DotPlot(stromal, features = unique(genelist[[1]]))+ RotatedAxis()

DotPlot(stromal, features = top10_stromal_coarse$gene)+ RotatedAxis()+
  coord_flip()

fibroblast_id_induced=WhichCells(stromal, ident='Fibroblasts')

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
stromal=SetIdent(stromal, value='Names')
muscle_id_stromal=WhichCells(stromal, ident='Skeletal muscle cells')


muscle_subset_stromal=subset(stromal, cells=muscle_id_stromal) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>%
  ScaleData(features = all.genes) %>%
  RunPCA(features = VariableFeatures(object = .)) %>% FindNeighbors(dims = 1:30)
muscle_subset_stromal=FindClusters(muscle_subset_stromal,resolution=0.2)
Idents(muscle_subset_stromal)
DimPlot(muscle_subset_stromal)
muscle_names_no=levels(Idents(muscle_subset_stromal))
muscle_markers_stromal=lapply(muscle_names_no, function(x) FindMarkers(muscle_subset_stromal, 
                                                                       ident.1=x, 
                                                                       ident.2=muscle_names_no[!muscle_names_no==x],
                                                                       logfc.threshold = 0.25,
                                                                       test.use='MAST'))
muscle_markers_gene_stromal=rbindlist(lapply(muscle_markers_stromal, function(x) mutate(x, gene=rownames(x)) %>% top_n(n=10, avg_logFC) %>% select(gene)))
muscle_subset_stromal %>% DotPlot(features=c(unique(muscle_markers_gene_stromal[[1]]),'Pax7'))+
  RotatedAxis()
  

muscle_markers[[1]][which(rownames(muscle_markers[[1]])=='Pax7'),]
DimPlot(muscle_subset_stromal, cols='Accent')
table(Idents(muscle_subset_stromal))

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
