##### Integrate cells from all skin samples (exclude adult hindlimb stromal cells) ##############
neonate_non_induced_2[['orig.ident']]='Neonate Ctrl'
neonate_induced[['orig.ident']]='Neonate FR'
adult_induced[['orig.ident']]='Adult FR'

skin.anchors <- FindIntegrationAnchors(object.list = list(neonate_non_induced_2, neonate_induced, adult_induced),
                                           k.filter = 50,dims = 1:30)
skin.integrated <- IntegrateData(anchorset = skin.anchors, dims = 1:30)

DefaultAssay(skin.integrated) <- "integrated"

skin.integrated[['Names']]=gsub('Myofibroblasts ','Myofibroblasts',skin.integrated[['Names']][[1]])
skin.integrated[['Names']]=gsub('Vascular smooth muscle cells','Myofibroblasts',skin.integrated[['Names']][[1]])
levels_integrated=c()
# Run the standard workflow for visualization and clustering
skin.integrated <- ScaleData(skin.integrated) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  RunPCA(features = VariableFeatures(object = .)) 
skin.integrated <- RunTSNE(skin.integrated, reduction = "pca", dims = 1:10)

skin.integrated=SetIdent(skin.integrated, value='orig.ident')
skin.integrated=RenameIdents(skin.integrated,'Neonate Ctrl'='Neonate -FR',
                               'Neonate FR'='Neonate +FR')
skin.integrated=StashIdent(skin.integrated, 'orig.ident')
integrated_levels=c("Epidermal cells",
                    "Fibroblasts",
                    "Macrophages",
                    "Melanocytes",
                    "Myofibroblasts",
                    "Schwann cells",
                    "Skeletal Muscle cells",
                    "Vascular cells",
                    "Lymph vessel cells",
                    "Endothelial cells",
                    "Sebocytes")
skin.integrated=SetIdent(skin.integrated, value='Names')
Idents(skin.integrated)=factor(Idents(skin.integrated), levels=integrated_levels,ordered = T)

DimPlot(skin.integrated, reduction = "umap", split.by='orig.ident',label = F, 
              repel = TRUE, cols = 'Paired', ncol=1, pt.size=0.4)

############## Integrate all muscle subsets from all 4 samples################
muscle_subset[['orig.ident']]='Neonate Ctrl'
muscle_subset_induced[['orig.ident']]='Neonate FR'
muscle_subset_induced_adult[['orig.ident']]='Adult FR'
muscle_subset_stromal[['orig.ident']]='Stromal'

muscle_list=lapply(list(muscle_subset, 
            muscle_subset_induced, 
            muscle_subset_induced_adult,
            muscle_subset_stromal), function(x) StashIdent(x,'muscle_clusters'))
muscle_all.anchors <- FindIntegrationAnchors(object.list = muscle_list,k.filter = 50,dims = 1:30)
muscle_all.integrated <- IntegrateData(anchorset = muscle_all.anchors, dims = 1:30)


DefaultAssay(muscle_all.integrated) <- "integrated"

for (x in unique(muscle_all.integrated[['orig.ident']])[[1]]){
muscle_all.integrated$muscle_clusters[muscle_all.integrated$orig.ident==x]=
  paste(x,muscle_all.integrated$muscle_clusters[muscle_all.integrated$orig.ident==x])}

#Get DE genes for all muscle tissues
de_genes_ctrl=lapply(muscle_markers, function(x) mutate(x, gene=rownames(x)) %>% filter(p_val_adj<0.01))
sapply(de_genes_ctrl, nrow)
DE_genes_list_neo_ctrl=rbindlist(lapply(muscle_markers, function(x) mutate(x, gene=rownames(x)) %>% top_n(n=500, avg_logFC) %>% select(gene)))
DE_genes_neonate_ctrl=unique(DE_genes_list_neo_ctrl[[1]])

de_genes_induced_adult=lapply(muscle_markers_induced_adult, function(x) mutate(x, gene=rownames(x)) %>% filter(p_val_adj<0.01))
sapply(de_genes_induced_adult, nrow)
DE_genes_list_induced_adult=rbindlist(lapply(muscle_markers_induced_adult, function(x) mutate(x, gene=rownames(x)) %>% top_n(n=50, avg_logFC) %>% select(gene)))
DE_genes_induced_adult=unique(DE_genes_list_induced_adult[[1]])

DE_genes_muscle_skin=unique(c(DE_genes_stromal, DE_genes_neonate, DE_genes_neonate_ctrl, DE_genes_induced_adult))

# Run the standard workflow for visualization and clustering
muscle_all.integrated <- ScaleData(muscle_all.integrated)
muscle_all.integrated <- RunPCA(muscle_all.integrated, features =  DE_genes_muscle_skin, verbose = FALSE)
muscle_all.integrated <- RunUMAP(muscle_all.integrated, reduction = "pca", dims = 1:10)
DimPlot(muscle_all.integrated, reduction = "umap", group.by = "muscle_clusters")
DimPlot(muscle_all.integrated, reduction = "umap" )
FeaturePlot(muscle_all.integrated, reduction = "umap", feature='Pax7')
VlnPlot(muscle_all.integrated, features =c('Pax7', 'Myod1'))

# Find clusters
muscle_all.integrated=  muscle_all.integrated %>% FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution=0.1, group.singletons = F) 
Idents(muscle_all.integrated)
muscle_all.integrated=StashIdent(muscle_all.integrated, 'integrated_clusters_7')
pData(HSMM_muscle)=muscle_all.integrated@meta.data
muscle_all.integrated_markers <- FindAllMarkers(muscle_all.integrated, only.pos = TRUE, 
                                        logfc.threshold = 0.25,
                                        test.use='MAST')

top10_muscle_all.integrated=muscle_all.integrated_markers %>%
  group_by(cluster) %>% 
  top_n(n = 8, wt = avg_logFC )

DotPlot(muscle_all.integrated, features = unique(c(top10_muscle_all.integrated$gene,'Pax7','Myod1')))+ 
  RotatedAxis()
#Make seurat object
seurat_object=muscle_all.integrated
data <- GetAssayData(seurat_object, slot='scale.data')
pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data)
head(seurat_object@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
HSMM_muscle <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              expressionFamily = uninormal())

#Compute trajectory
HSMM_muscle <- setOrderingFilter(HSMM_muscle, DE_genes_muscle_skin)
#HSMM_muscle <- setOrderingFilter(HSMM_muscle, VariableFeatures(muscle_all.integrated))

HSMM_muscle <- reduceDimension(HSMM_muscle, max_components = 2,
                               method = 'DDRTree', norm_method ='none', pseudo_expr =0, scaling=F,
                               relative_expr =F)
HSMM_muscle <- orderCells(HSMM_muscle)
HSMM_muscle <- orderCells(HSMM_muscle, root_state = 4)

#Get DE genes among pseudotime clusters
pseudotime_cluster_tb=pData(HSMM_muscle) %>% 
  rownames_to_column() %>%
  filter((State == 1 & Pseudotime>15)|(State == 5 & Pseudotime>18)|(State == 4 & Pseudotime<2)) %>%
  group_by(State) 

muscle_all.integrated_sub=subset(muscle_all.integrated, cells=pseudotime_cluster_tb$rowname)
Idents(muscle_all.integrated_sub)=pseudotime_cluster_tb$State
muscle_all.integrated_sub_markers=FindAllMarkers(muscle_all.integrated_sub, only.pos = TRUE, 
               logfc.threshold = 0.25,
               test.use='MAST')
top20_pseudtimeclusters=muscle_all.integrated_sub_markers %>%
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_logFC )


#Plot trajectory
plot_cell_trajectory(HSMM_muscle,
                     show_tree=F,
                     show_backbone = F,
                     show_branch_points=F,
                     markers=c('Mki67','Pax7','Cdkn1a','Myog'),
                     use_color_gradient=T,
                     markers_linear=T,
                     cell_size=0.5,
                     option='B')+theme(text=element_text(size=20))

HSMM_muscle_BEAM=BEAM(HSMM_muscle, cores = 1)
trace(plot_genes_branched_heatmap, edit='T')

plot_genes_branched_heatmap(HSMM_muscle_2[top20_pseudtimeclusters$gene,],
                            num_clusters = 3,
                            cores = 1,
                            branch_labels=c('Differentiating','Proliferating'),
                            show_rownames = T,
                            norm_method='log')

##### Compare -FR and + FR pseudotime clusters
subset_earlier_pseudotime=pData(HSMM_muscle) %>% 
  rownames_to_column() %>%
  filter(State==1&(orig.ident=="Neonate Ctrl"|orig.ident=="Neonate FR")) %>% 
  filter(Pseudotime<1) %>%
  select(rowname)

#Integrated
pseudotime_0=subset(muscle_all.integrated, cells=subset_earlier_pseudotime[[1]])
pseudotime_0=SetIdent(pseudotime_0, value='orig.ident')
pseudotime_DE=FindAllMarkers(pseudotime_0, only.pos = T, 
               test.use='MAST')
pseudotime_DE %>% group_by(cluster) %>% tally
DoHeatmap(pseudotime_0, features=pseudotime_DE$gene)

#Merged
neonatal_muscle.merged=merge(muscle_subset, muscle_subset_induced)
neonatal_muscle.merged=ScaleData(neonatal_muscle.merged)
pseudotime_0_merged=subset(neonatal_muscle.merged, cells=subset_earlier_pseudotime[[1]])
pseudotime_0_merged=SetIdent(pseudotime_0_merged, value='orig.ident')
pseudotime_DE_merged=FindAllMarkers(pseudotime_0_merged, only.pos = T, 
                             test.use='MAST')
write.csv(filter(pseudotime_DE_merged, avg_logFC>0.5), 'pseudotime_DE_merged.txt')
pseudotime_DE_merged %>% filter(avg_logFC>0.5) %>%
  group_by(cluster) %>% tally
top_genes_pseudotime_0=pseudotime_DE_merged %>% group_by(cluster) %>% top_n(20, wt=avg_logFC)

DotPlot(pseudotime_0_merged, features =top_genes_pseudotime_0$gene)+ 
  RotatedAxis()+theme(axis.text=element_text(size=12))
table(pseudotime_0_merged@meta.data$orig.ident)

##### Compare -FR and + FR pseudotime clusters between endogenous and induced proliferating

unique(pData(HSMM_muscle)['orig.ident'])
subset_earlier_pseudotime=pData(HSMM_muscle) %>% 
  rownames_to_column() %>%
  filter(State==5&(orig.ident=="Stromal"|orig.ident=="Neonate FR")) %>% 
  filter(Pseudotime>15) 

#Merged
neonatal_muscle.merged=merge(muscle_subset_stromal, muscle_subset_induced)
neonatal_muscle.merged=ScaleData(neonatal_muscle.merged)
pseudotime_0_merged=subset(neonatal_muscle.merged, cells=subset_earlier_pseudotime[[1]])
pseudotime_0_merged=SetIdent(pseudotime_0_merged, value='orig.ident')
pseudotime_DE_merged=FindAllMarkers(pseudotime_0_merged, only.pos = T, 
                                    test.use='MAST')
write.csv(filter(pseudotime_DE_merged, avg_logFC>0.5), 'pseudotime_DE_merged.txt')
pseudotime_DE_merged %>% filter(avg_logFC>0.5) %>%
  group_by(cluster) %>% tally
top_genes_pseudotime_0=pseudotime_DE_merged %>% group_by(cluster) %>% top_n(20, wt=avg_logFC)

DotPlot(pseudotime_0_merged, features =top_genes_pseudotime_0$gene)+ 
  RotatedAxis()+theme(axis.text=element_text(size=11))
table(pseudotime_0_merged@meta.data$orig.ident)

#######goseq
options(max.print = 50)
library(goseq)
pseudotime_DE
#de_genes=as.integer(rownames(pseudotime_0_merged) %in% 
#                      filter(pseudotime_DE, cluster=='Neonate FR')$gene)
de_genes=as.integer(rownames(pseudotime_0_merged) %in% 
                      filter(pseudotime_DE_merged, cluster=='Neonate FR', avg_logFC>0.5)$gene)
names(de_genes)=rownames(pseudotime_0_merged)
sum(de_genes)
pwf=nullp(de_genes,"mm9","geneSymbol")
GO.wall=goseq(pwf,"mm9","geneSymbol", test.cats='GO:BP', method='Hypergeometric')
GO.wall$padjusted=p.adjust(GO.wall$over_represented_pvalue)

write.csv(GO.wall[GO.wall$padjusted<0.01,], 'GO_results_proliferating.csv')
GO.wall[GO.wall$padjusted<0.01,]

#######Compare -FR and +FR all cells after integration
FR_names=c('Neonate FR', 'Neonate Ctrl')

FR_allcells_markers=FindMarkers(SetIdent(muscle_all.integrated, value='orig.ident'), 
                                 ident.1=FR_names[1], 
                                 ident.2=FR_names[2],
                                 logfc.threshold = 0,
                                 test.use='MAST')

rownames_to_column(FR_allcells_markers) %>% filter(p_val_adj<0.01&avg_logFC>0.25)
