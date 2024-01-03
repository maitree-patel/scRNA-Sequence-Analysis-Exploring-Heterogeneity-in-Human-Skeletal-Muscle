library(Seurat)
library(SeuratDisk)
library(tidyverse)

muscle <- Read10X(data.dir = "/Users/maitreepatel/Desktop/R/scRNA analysis/GSE214544_RAW")
head(muscle)
str(muscle)

muscle.seurat <- CreateSeuratObject(counts = muscle,
                                    project = "skeletal_muscle",
                                    min.cells = 3,
                                    min.features = 200)
str(muscle.seurat)

#1. Quality Control
muscle.seurat[["percent.mt"]] <- PercentageFeatureSet(muscle.seurat,
                                                      pattern = "^MT-")
muscle.seurat@meta.data
#visualizing QC metrics
VlnPlot(muscle.seurat,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

FeatureScatter(muscle.seurat,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")
#some cells with really high counts
#we have a correlation of 0.9

FeatureScatter(muscle.seurat,
               feature1 = "nCount_RNA",
               feature2 = "percent.mt") +
  geom_smooth(method = "lm")

muscle.seurat <- subset(muscle.seurat,
                        subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

#viewing metrics plot again
VlnPlot(muscle.seurat,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

FeatureScatter(muscle.seurat,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

#2. Normalizing data
muscle.seurat <- NormalizeData(muscle.seurat,
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)

#3. Identifying highly variable features
muscle.seurat <- FindVariableFeatures(muscle.seurat)
top10 <- head(VariableFeatures(muscle.seurat), 10)  

#visualizing the highly variable features, with top 10 labelled
plot.var <- VariableFeaturePlot(muscle.seurat)
LabelPoints(plot = plot.var,
            points = top10,
            repel = TRUE)

#4. scaling data
#getting all genes
all.features <- rownames(muscle.seurat)

muscle.seurat <- ScaleData(muscle.seurat,
                           features = all.features)

#5. Performing linear dimensionality reduction
muscle.seurat <- RunPCA(muscle.seurat,
                        features = VariableFeatures(muscle.seurat))

DimPlot(muscle.seurat,
        reduction = "pca")

DimHeatmap(muscle.seurat,
           dims = 1:15,
           cells = 500,
           balanced = TRUE)

ElbowPlot(muscle.seurat)

#clustering cells
muscle.seurat <- FindNeighbors(muscle.seurat,
                               dims = 1:20)
muscle.seurat <- FindClusters(muscle.seurat,
                              resolution = 0.05)
muscle.seurat@meta.data 
#checking which resolution to keep through plotting

DimPlot(muscle.seurat,
        group.by = "RNA_snn_res.0.05",
        reduction = "pca",
        label = TRUE)

muscle.seurat <- RunUMAP(muscle.seurat,
                         dims = 1:20)

DimPlot(muscle.seurat,
        reduction = "umap")
#extremely heterogenous muscle
#total of 23 cell clusters

#Indentifying cluster markers
cluster0.markers <- FindMarkers(muscle.seurat, 
                                ident.1 = 0,
                                logfc.threshold = 0.25, 
                                test.use = "roc", 
                                only.pos = TRUE)

head(cluster0.markers, 
     n = 5)

VlnPlot(muscle.seurat,
        features = c(rownames(cluster0.markers)[1],
                     rownames(cluster0.markers)[2]))

cluster1.markers <- FindMarkers(muscle.seurat, 
                                ident.1 = 1,
                                logfc.threshold = 0.25, 
                                test.use = "roc", 
                                only.pos = TRUE)

VlnPlot(muscle.seurat,
        features = c(rownames(cluster1.markers)[1],
                     rownames(cluster1.markers)[2]))

#finding all markers
muscle.allmarkers <- FindAllMarkers(muscle.seurat,
                                    only.pos = T,
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25)

x <- muscle.allmarkers %>%
  group_by(cluster) %>% #grouping by cluster
  top_n(n = 1, wt = avg_log2FC) #visualized by the top gene based on logfc (i.e. higherst expressed gene)

FeaturePlot(muscle.seurat,
            features = x$gene)

#cluster 4
cluster4.markers <- FindMarkers(muscle.seurat, 
                                ident.1 = 4,
                                logfc.threshold = 0.25, 
                                test.use = "roc", 
                                only.pos = TRUE)

VlnPlot(muscle.seurat,
        features = c(rownames(cluster4.markers)[1],
                     rownames(cluster4.markers)[2]))

FeaturePlot(muscle.seurat,
            features = "TNNT1")
head(cluster4.markers,
     n = 5)
which(rownames(cluster4.markers)=="MYL2")


#cluster 9 
cluster9.markers <- FindMarkers(muscle.seurat, 
                                ident.1 = 9,
                                logfc.threshold = 0.25, 
                                test.use = "roc", 
                                only.pos = TRUE)
VlnPlot(muscle.seurat,
        features = c(rownames(cluster9.markers)[1],
                     rownames(cluster9.markers)[2]))

FeaturePlot(muscle.seurat,
            features = "FTH1")
FeaturePlot(muscle.seurat,
            features = "PROK2")
head(cluster9.markers,
     n = 20)


VlnPlot(muscle.seurat,
        features = "PROK2")
which(rownames(cluster9.markers)=="PROK2")

#cluster 8
cluster8.markers <- FindMarkers(muscle.seurat, 
                                ident.1 = 8,
                                logfc.threshold = 0.25, 
                                test.use = "roc", 
                                only.pos = TRUE)
VlnPlot(muscle.seurat,
        features = c(rownames(cluster8.markers)[1],
                     rownames(cluster8.markers)[2])) 
head(cluster8.markers,
     n = 20)
