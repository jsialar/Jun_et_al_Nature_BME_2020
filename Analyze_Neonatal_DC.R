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
ldat_non_induced <- ReadVelocity(file = "/workdir/jsia1/chemical_reprog/Neo_Fibroblasts_no_chem_processed/velocyto/Neo_Fibroblasts_no_chem_processed.loom")
x_non_induced <- as.Seurat(x = ldat_non_induced)
neonate_non_induced_2 <- subset(x_non_induced, subset = nFeature_spliced  > 1000 & nCount_spliced   < 50000) %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = all.genes) %>%
  RunPCA(features = VariableFeatures(object = .)) 

neonate_non_induced_2=  neonate_non_induced_2%>% FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution=0.5) %>% RunTSNE(dims = 1:10)

FeaturePlot(neonate_non_induced_2, features =c('Pax7','Top2a','Myh3','Myog'))
DimPlot(neonate_non_induced_2, label=F)

neonate_non_induced_2_markers <- FindAllMarkers(neonate_non_induced_2, only.pos = TRUE, 
                                              logfc.threshold = 0.25,
                                              test.use='MAST')
top10_neonate_non_induced=neonate_non_induced_2_markers %>%
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_logFC)

#DoHeatmap(neonate_non_induced, features = top10_neonate_non_induced$gene) + RotatedAxis() 
DotPlot(neonate_non_induced_2, features = c(unique(top10_neonate_non_induced$gene),'Pax7','Pdgfra'))+ RotatedAxis()

filter(neonate_non_induced_2_markers, gene=='Pax7')
filter(top10_neonate_non_induced, cluster==8)


#neonate_non_induced_2=StashIdent(neonate_non_induced_2, 'clusters_15')
neonate_non_induced_2=SetIdent(neonate_non_induced_2, value = 'clusters_15')
neonate_non_induced_2=RenameIdents(neonate_non_induced_2, 
                                   '15'='Macrophages',
                                   '14'='Epidermal cells',
                                   '13'='Melanocytes',
                                   '12'='Fibroblasts',
                                   '11'='Vascular cells',
                                   '10'='Schwann cells',
                                   '9'='Epidermal cells',
                                   '8'='Skeletal Muscle cells',
                                   '7'='Fibroblasts',
                                   '6'='Fibroblasts',
                                   '5'='Myofibroblasts',
                                   '4'='Fibroblasts',
                                   '3'='Fibroblasts',
                                   '2'='Epidermal cells',
                                   '1'='Fibroblasts',
                                   '0'='Fibroblasts')
neonate_non_induced_2[['Names']]=Idents(neonate_non_induced_2)
cat_counts=as.vector(table(Idents(neonate_non_induced_2)))
cat_per=round(cat_counts/sum(cat_counts)*100, 1)
new_names=paste0(levels(Idents(neonate_non_induced_2)),' ~ ',cat_per,'%','')
names(new_names)=levels(Idents(neonate_non_induced_2))
neonate_non_induced_2=RenameIdents(neonate_non_induced_2, new_names)
Idents(neonate_non_induced_2)=factor(Idents(neonate_non_induced_2), levels=sort(levels(Idents(neonate_non_induced_2))),ordered = T)
neonate_induced=StashIdent(neonate_non_induced_2,'Names_per')
neonate_induced=SetIdent(neonate_non_induced_2, value='Names_per')
DimPlot(neonate_non_induced_2,reduction='tsne', cols='Paired', label=F)




neonate_non_induced_2=SetIdent(neonate_non_induced_2, value='Names')

neonate_non_induced_2_names_markers <- FindAllMarkers(neonate_non_induced_2, only.pos = TRUE, 
                                                logfc.threshold = 0.25,
                                                test.use='MAST')
top10_neonate_non_induced=neonate_non_induced_2_names_markers %>%
  group_by(cluster) %>% 
  top_n(n = 3, wt = avg_logFC)
filter(top10_neonate_non_induced, cluster=='Myofibroblasts')
filter(neonate_non_induced_2_names_markers, gene=='Mylk')
top10_neonate_non_induced=mutate(top10_neonate_non_induced,gene=ifelse(gene=='Igfbp7','Mylk',gene))  %>%
  mutate(gene=ifelse(gene=='H19','Eng',gene))

neonate_non_induced_2=SetIdent(neonate_non_induced_2, value='Names')
# neonate_non_induced_2=RenameIdents(neonate_non_induced_2,'Epidermal cells'='Epi',
#              'Fibroblasts'='Fib',
#              'Macrophages'='Mac',
#              'Melanocytes'='Mel',
#              'Myofibroblasts'='Myo',
#              'Schwann cells'='Schw',
#              'Skeletal Muscle cells'='SkM',
#              'Vascular cells'='Vas')

Idents(neonate_non_induced_2)=factor(Idents(neonate_non_induced_2), levels=sort(levels(Idents(neonate_non_induced_2))),ordered = T)
DotPlot(neonate_non_induced_2, features = genelist_overall)+ RotatedAxis()

#top10_neonate_non_induced$cluster=as.character(top10_neonate_non_induced$cluster)
#top10_neonate_non_induced=arrange(top10_neonate_non_induced, cluster)
#DotPlot(neonate_non_induced_2, features = unique(top10_neonate_non_induced$gene))+ RotatedAxis()


############## Fibroblast sub
fibroblast_id=WhichCells(neonate_non_induced_2, ident='Fibroblasts')

fibrobast_subset=subset(neonate_non_induced_2, cells=fibroblast_id) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>%
  ScaleData(features = all.genes) %>%
  RunPCA(features = VariableFeatures(object = .)) %>% FindNeighbors(dims = 1:30)
fibrobast_subset=FindClusters(fibrobast_subset,resolution=0.1)
table(Idents(fibrobast_subset))
DimPlot(fibrobast_subset, label = T)
fibrobast_subset=subset(fibrobast_subset, cells=WhichCells(fibrobast_subset, idents = 4, invert = T))
fibroblasts_names_no=levels(Idents(fibrobast_subset))
fibroblast_markers=lapply(fibroblasts_names_no, function(x) FindMarkers(fibrobast_subset, 
                                                  ident.1=x, 
                                                  ident.2=fibroblasts_names_no[!fibroblasts_names_no==x],
                                                  logfc.threshold = 0.25,
                                                  test.use='MAST'))
library(data.table)
fibroblast_markers_gene=rbindlist(lapply(fibroblast_markers, function(x) mutate(x, gene=rownames(x)) %>% top_n(n=3, avg_logFC) %>% select(gene)))
fibrobast_subset %>% DotPlot(features=unique(fibroblast_markers_gene[[1]]))+
  RotatedAxis()+
  coord_flip()
DimPlot(fibrobast_subset, cols='Set1')


######### Muscle sub###############
neonate_non_induced_2=SetIdent(neonate_non_induced_2, value='Names')
muscle_id=WhichCells(neonate_non_induced_2, ident='Skeletal Muscle cells')


muscle_subset=subset(neonate_non_induced_2, cells=muscle_id) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>%
  ScaleData(features = all.genes) %>%
  RunPCA(features = VariableFeatures(object = .)) %>% FindNeighbors(dims = 1:30)
muscle_subset=FindClusters(muscle_subset,resolution=0.2)
Idents(muscle_subset)
muscle_names_no=levels(Idents(muscle_subset))
muscle_markers=lapply(muscle_names_no, function(x) FindMarkers(muscle_subset, 
                                                                        ident.1=x, 
                                                                        ident.2=muscle_names_no[!muscle_names_no==x],
                                                                        logfc.threshold = 0.25,
                                                                        test.use='MAST'))
muscle_markers_gene=rbindlist(lapply(muscle_markers, function(x) mutate(x, gene=rownames(x)) %>% top_n(n=5, avg_logFC) %>% select(gene)))
muscle_subset %>% DotPlot(features=c('Pax7',unique(muscle_markers_gene[[1]]) ))+
  RotatedAxis()+
  coord_flip()

muscle_markers[[1]][which(rownames(muscle_markers[[1]])=='Pax7'),]
DimPlot(muscle_subset, cols='Set1')
table(Idents(muscle_subset))
##### Find unique markers for progenitor cells#######
differentiating_muscle_id=WhichCells(muscle_subset,idents =1)
neonate_non_induced_nomuscle=subset(neonate_non_induced_2, cells=differentiating_muscle_id, invert=T)
idents_nomuscle=as.character(unique(neonate_non_induced_nomuscle[['Names']][[1]]))
idents_nomuscle_=idents_nomuscle[idents_nomuscle!='Skeletal Muscle cells']

muscle_progenitor_markers=lapply(idents_nomuscle_, function(x) FindMarkers(neonate_non_induced_nomuscle, 
                                                               ident.1='Skeletal Muscle cells', 
                                                               ident.2=x,
                                                               logfc.threshold = 0.25,
                                                               test.use='MAST'))
muscle_progenitor_markers=lapply(muscle_progenitor_markers, function(x) {x$gene=rownames(x) 
x})
names(muscle_progenitor_markers)=idents_nomuscle_
muscle_progenitor_markers_=rbindlist(muscle_progenitor_markers, idcol='cluster')

candidates=muscle_progenitor_markers_ %>%
  filter(p_val_adj<0.01) %>%
  group_by(gene) %>%
  mutate(n=n(), max_p=max(p_val_adj ), min_logfc=min(avg_logFC)) %>%
  ungroup() %>%
  filter(n==7, min_logfc>1)
DoHeatmap(neonate_non_induced_nomuscle, features = candidates$gene) 


neonate_non_induced_2=SetIdent(neonate_non_induced_2, value = 'clusters_15')
neonate_non_induced_2=RenameIdents(neonate_non_induced_2, 
                                   '15'='Macrophages',
                                   '14'='Epidermal cells 1',
                                   '13'='Melanocytes',
                                   '12'='Fibroblasts 1',
                                   '11'='Vascular cells',
                                   '10'='Schwann cells',
                                   '9'='Epidermal cells 2',
                                   '8'='Skeletal Muscle cells',
                                   '7'='Fibroblasts 2',
                                   '6'='Fibroblasts 3',
                                   '5'='Myofibroblasts',
                                   '4'='Fibroblasts 4',
                                   '3'='Fibroblasts 5',
                                   '2'='Epidermal cells 3',
                                   '1'='Fibroblasts 6',
                                   '0'='Fibroblasts 7')

neonate_non_induced_2=StashIdent(neonate_non_induced_2, 'Names_sub')
neonate_non_induced_2=SetIdent(neonate_non_induced_2, value='Names_sub')
neonate_non_induced_2=RenameIdents(neonate_non_induced_2, 'Skeletal Muscle cells' = 'SkMu (Progenitors)')
Idents(neonate_non_induced_2, cells=differentiating_muscle_id)='SkMu (Differentiating)'
neonate_non_induced_2=StashIdent(neonate_non_induced_2, 'Names_sub_muscles')
table(Idents(neonate_non_induced_2))
idents_nomuscle_2=as.character(unique(neonate_non_induced_2[['Names_sub_muscles']][[1]]))
idents_nomuscle_2_=idents_nomuscle_2[idents_nomuscle_2!='SkMu (Progenitors)']
neonate_non_induced_2=SetIdent(neonate_non_induced_2, value='Names_sub_muscles')
muscle_progenitor_markers_2=lapply(idents_nomuscle_2_, function(x) FindMarkers(neonate_non_induced_2, 
                                                                           ident.1='SkMu (Progenitors)', 
                                                                           ident.2=x,
                                                                           logfc.threshold = 0.25,
                                                                           test.use='MAST'))
muscle_progenitor_markers_2=lapply(muscle_progenitor_markers_2, function(x) {x$gene=rownames(x) 
x})
neonate_non_induced_2=SetIdent(neonate_non_induced_2, value='Names_sub_muscles')
abbreved=gsub('Fibroblasts','Fb',Idents(neonate_non_induced_2)) %>%
  gsub('Epidermal cells','Epi',.) %>%
  gsub('Schwann cells','Schw',.) %>%
  gsub('Vascular cells','Vas',.) %>%
  gsub('Myofibroblasts','Myo',.) %>%
  gsub('Macrophages','Mac',.) %>%
  gsub('Melanocytes','Mel',.) %>%
  gsub('SkMu \\(Progenitors\\)','SkM(Prog)',.) %>%
  gsub('SkMu \\(Differentiating\\)','SkM(Diff)',.)
  
neonate_non_induced_2[['Names_sub_muscles_abrev']]=abbreved
names(muscle_progenitor_markers_2)=idents_nomuscle_2_
muscle_progenitor_markers_2_=rbindlist(muscle_progenitor_markers_2, idcol='cluster')

candidates_2=muscle_progenitor_markers_2_ %>%
  filter(p_val_adj<0.01) %>%
  group_by(gene) %>%
  mutate(n=n(), max_p=max(p_val_adj ), min_logfc=min(avg_logFC)) %>%
  ungroup() %>%
  filter(n==16, min_logfc>0.4) %>% 
  distinct(gene, .keep_all = T) %>%
  arrange(desc(min_logfc),max_p)

muscle_progenitor_markers_2_ %>% 
  filter(gene %in% candidates_2$gene, gene=='Pax7')
options(max.print=100)

muscle_progenitor_markers_2_%>%
  filter(gene=='Pax7')
DoHeatmap(neonate_non_induced_2, features = candidates_2$gene, disp.min = -Inf, disp.max = Inf) 
DotPlot(neonate_non_induced_2, features = candidates_2$gene)+ RotatedAxis()

#Manual violin plot
violin_plot=as.data.frame(GetAssayData(neonate_non_induced_2)) %>% 
  rownames_to_column('gene') %>%
  filter(gene %in% candidates_2$gene)
violin_plot_2=violin_plot%>%
  gather(key='cellid', value='expression',-gene) %>%
  inner_join(y=neonate_non_induced_2@meta.data['Names_sub_muscles'] %>% rownames_to_column('cellid'), 
             by='cellid')
ggplot(data=violin_plot_2, aes(x=Names_sub_muscles, y=expression))+geom_violin(scale = 'area')+facet_grid(gene~.)

FeaturePlot(neonate_non_induced_2, reduction='tsne',features = candidates_2$gene)
FeaturePlot(neonate_non_induced_2, reduction='tsne',
            cols=c('green','red'),
            features = candidates_2$gene)

neonate_non_induced_2=SetIdent(neonate_non_induced_2,value='Names_sub_muscles_abrev')
vinplot_list=VlnPlot(neonate_non_induced_2, features = candidates_2$gene, pt.size=0, ncol =1, sort=F, combine=F)
vinplot_list2=lapply(vinplot_list, function(x)x+theme(legend.position = "none",
                                                      axis.title= element_blank())+
                       scale_x_discrete(limits=sort(levels(Idents(neonate_non_induced_2))))+
                       coord_flip())

library(cowplot)
plot_grid(plotlist=vinplot_list2, nrow=1)

###### RNA velocity #############

fibrobast_subset_vel <- RunVelocity(object=fibrobast_subset, ncores = 32)
ident.colors <- RColorBrewer::brewer.pal(length(levels(fibrobast_subset_vel)), 'Set1')
names(ident.colors) <- levels(fibrobast_subset_vel)
cell.colors <- ident.colors[Idents(fibrobast_subset_vel)]
names(cell.colors) <- colnames(fibrobast_subset_vel)
embedding_results=show.velocity.on.embedding.cor(emb = Embeddings(object = fibrobast_subset_vel, reduction = "umap"), 
                                                 vel = Tool(object = fibrobast_subset_vel, slot = "RunVelocity"), 
                               cell.colors = ac(x = cell.colors, alpha = 0.5), n.cores =32,
                               cex = 0.8, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 100,  
                               do.par = FALSE, cell.border.alpha = 0.1)

show.velocity.on.embedding.cor(emb = Embeddings(object = fibrobast_subset_vel, reduction = "umap"), 
                               vel = Tool(object = fibrobast_subset_vel, slot = "RunVelocity"), 
                               cell.colors = ac(x = cell.colors, alpha = 0.2), scale='sqrt', n.cores =32,
                               cex = 0.5, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 50,  
                               do.par = FALSE, cell.border.alpha = 0.1, arrow.scale=2,cc=embedding_results$cc)

muscle_subset_vel <- RunVelocity(object=muscle_subset, ncores = 32)
ident.colors <- RColorBrewer::brewer.pal(length(levels(muscle_subset_vel)), 'Set1')
names(ident.colors) <- levels(muscle_subset_vel)
cell.colors <- ident.colors[Idents(muscle_subset_vel)]
names(cell.colors) <- colnames(muscle_subset_vel)
embedding_results=show.velocity.on.embedding.cor(emb = Embeddings(object = muscle_subset_vel, reduction = "umap"), 
                                                 vel = Tool(object = muscle_subset_vel, slot = "RunVelocity"), 
                                                 cell.colors = ac(x = cell.colors, alpha = 0.2), n.cores =32,
                                                 cex = 0.5, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 50,  
                                                 do.par = FALSE, cell.border.alpha = 0.1, arrow.scale=2)

show.velocity.on.embedding.cor(emb = Embeddings(object = muscle_subset_vel, reduction = "umap"), 
                               vel = Tool(object = muscle_subset_vel, slot = "RunVelocity"), 
                               cell.colors = ac(x = cell.colors, alpha = 0.2), scale='sqrt', n.cores =32,
                               cex = 0.5, show.grid.flow = T, min.grid.cell.mass = 0.5, grid.n = 50,  
                               do.par = FALSE, cell.border.alpha = 0.1, arrow.scale=2,cc=embedding_results$cc)
