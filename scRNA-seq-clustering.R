library(Seurat)
library(tidyverse)
library(Matrix)


cells_per_sample = function(object) {
  # Create a dataframe from a table of number of cells per SampleID to create
  # a bar plot
  df = as.data.frame(rev(table(seurat.obj$SampleID)))
  colnames(df) = c('SampleID', 'Number_of_Cells')
  ggplot(df, aes(x=reorder(SampleID, -Number_of_Cells), y=Number_of_Cells, 
                 fill=SampleID)) +
    geom_bar(stat='identity') +
    scale_y_continuous(expand=c(0,0)) +
    NoLegend() + 
    RotatedAxis() + 
    ylab(expression(N[cells])) + 
    xlab('Sample ID') + 
    ggtitle(paste('Number of Cells:', sum(df$Number_of_Cells)))
}

# collect the matrix files from the alignment output
align.output.dir = 'align_out_files/'
files = list.files(align.output.dir)
# collect the file base names for all the samples
base.names = str_replace_all(files[grepl('_matrix.mtx', files)], 
                             '_matrix.mtx', '')
# make Seurat objects of all the samples, then combine them
object.list = lapply(base.names, function(file) {
  # read the sparse matrix, barcodes, and genes files
  mat.file = readMM(paste0(align.output.dir, file, '_matrix.mtx'))
  barcodes = read.csv(paste0(align.output.dir, file, '_barcodes.tsv'), header=F,
                      sep='\t')
  genes = read.csv(paste0(align.output.dir, file, '_genes.tsv'), header=F, 
                   sep='\t')
  
  # assign col and row names to the matrix file from the genes and barcodes
  colnames(mat.file) = barcodes$V1
  rownames(mat.file) = genes$V2
  
  current.object = CreateSeuratObject(counts = mat.file, project = 'scRNA-seq')
  current.object@meta.data$SampleID = file
  current.object
})

# combine all the objects into one
seurat.obj = merge(x=object.list[[1]], y=object.list[2:length(object.list)])

# visualize the count of cells per sample
cells_per_sample(seurat.obj)

# Visualize based upon mitochondrial gene content, 
# genes per cell (nFeatures_RNA), and UMIs per cell (nCount_RNA) to do QC
seurat.obj[["mitochondrial.content"]] = PercentageFeatureSet(seurat.obj, 
                                                             pattern = "^MT-")
VlnPlot(seurat.obj, features='mitochondrial.content', group.by='SampleID', 
        ncol=1, pt.size=0)
VlnPlot(seurat.obj, features='nFeature_RNA', group.by='SampleID', ncol=1, 
        pt.size=0)
VlnPlot(seurat.obj, features='nCount_RNA', group.by='SampleID', ncol=1, 
        pt.size=0)

# subset based upon violin plots and re-plot cells per sample to see the change
seurat.obj = subset(seurat.obj, nCount_RNA >= 300 & nCount_RNA <= 2750)
seurat.obj = subset(seurat.obj, mitochondrial.content >= 15)

cells_per_sample(seurat.obj)

# perform normalization, scaling, and feature selection for further analysis
seurat.obj = NormalizeData(seurat.obj)
seurat.obj = ScaleData(seurat.obj, features=rownames(seurat.obj))
seurat.obj = FindVariableFeatures(seurat.obj, nfeatures=2000)
# clean unused memory
gc()

# plot variable features
LabelPoints(VariableFeaturePlot(seurat.obj), 
            points=head(VariableFeatures(seurat.obj), 15), repel=T, xnudge=0, 
            ynudge=0)

# run PCA and visualize
seurat.obj = RunPCA(seurat.obj, features=VariableFeatures(seurat.obj), npcs=50)
DimPlot(seurat.obj, reduction='pca', group.by='SampleID')
# elbow plot to choose number of PCs for further analysis
ElbowPlot(seurat.obj, ndims=50)

# cluseter and UMAP non-linear dimensionality reduction
seurat.obj = FindNeighbors(seurat.obj, dims=1:22)
seurat.obj = FindClusters(seurat.obj)
seurat.obj = RunUMAP(seurat.obj, dims=1:22)
DimPlot(seurat.obj, reduction='umap', group.by='seurat_clusters', label=T) + 
  NoLegend() + 
  ggtitle('UMAP Dimensionality Reduction by Clusters')

# Identify cluster biomarkers, save them and plot a heatmap
markers = FindAllMarkers(seurat.obj, logfc.threshold=2, test.use='MAST')
write.table(markers, 'biomarkers.csv', row.names = F)
top = markers %>% group_by(cluster) %>% top_n(5, wt=avg_log2FC) %>% .$gene

DoHeatmap(seurat.obj, features=top, group.by='seurat_clusters', label=F)
