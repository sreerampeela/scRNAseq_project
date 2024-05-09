## script reference @ 
# https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_standard_workflow.R
## using R version 4.3.3


## loading libraries
library(Seurat)
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(gridExtra)

#set the working directory
setwd("F:/scRNAseq_learning")

# load the 20k data file
sample_data <- Read10X_h5("20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
str(sample_data) # overview of the object
cts <-  sample_data$`Gene Expression` # extract gene counts

# creating Seurat object
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC_project", 
                                       min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj


## QC
View(nsclc.seurat.obj@meta.data)
# % MT reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, 
                                                         pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
        , ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering -----------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & 
                             nFeature_RNA < 2500 & percent.mt < 5)
# 3. Normalize data
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)

# 4. Identify highly variable features --------------
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, 
                                         selection.method = "vst", 
                                         nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
##error when setting repel = T; 
#used suggestions from git issue https://github.com/satijalab/seurat/issues/1590
LabelPoints(plot = plot1, points = top10, repel = F, xnudge = 0, ynudge = 0)


# 5. Scaling -------------
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

str(nsclc.seurat.obj)

# 6. Perform Linear dimensionality reduction --------------
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, 
                           features = VariableFeatures(object = nsclc.seurat.obj))

# visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)


# determine dimensionality of the data
ElbowPlot(nsclc.seurat.obj)


# 7. Clustering ------------
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)


# understanding resolution
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)

plt1 <- DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

## trying with 10 dims
nsclc.seurat.obj10 <- FindNeighbors(nsclc.seurat.obj, dims = 1:10)
nsclc.seurat.obj10 <- FindClusters(nsclc.seurat.obj10, resolution = c(0.1,0.3, 0.5, 0.7, 1))
plt2 <- DimPlot(nsclc.seurat.obj10, group.by = "RNA_snn_res.0.5", label = TRUE)

plt4 <- DimPlot(nsclc.seurat.obj10, group.by = "RNA_snn_res.0.3", label = TRUE)
plt3 <- DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.3", label = TRUE)

plt5 <- DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.7", label = TRUE)
plt6 <- DimPlot(nsclc.seurat.obj10, group.by = "RNA_snn_res.0.7", label = TRUE)

plt7 <- DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)
plt8 <- DimPlot(nsclc.seurat.obj10, group.by = "RNA_snn_res.0.1", label = TRUE)

grid.arrange(plt7, plt8, plt1, plt2, plt3, plt4, plt5, plt6, ncol = 2, nrow = 4)
# setting identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)

# non-linear dimensionality reduction --------------
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
# other dimensions
nsclc.seurat.obj10 <- RunUMAP(nsclc.seurat.obj10, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
plt_umap1 <- DimPlot(nsclc.seurat.obj, reduction = "umap", label = T)
plt_umap2 <- DimPlot(nsclc.seurat.obj10, reduction = "umap", label = T)
grid.arrange(plt_umap1,plt_umap2 , ncol = 2, nrow = 1)
