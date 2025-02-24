# scRNA Sequence Analysis Exploring Heterogeneity in Human Skeletal Muscle
scRNA-Sequence Analysis was performed as a way to analyze the transcriptional profile of human skeletal muscles using Seurat package in R.

## Objective 
scRNA-Sequencing differs from Bulk RNA Sequencing in providing a finer resolution of gene expression within cells enabling us to explore  the hetergeneity between cells. From analysing baseline transcriptome of the Human skeletal muscles (https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles), scRNA-seq analysis was further used to expllore the single-cell transcriptome. The main aim was to understand the heterogeneity in skeletal muscle cells by loading, pre-processing, clustering and marker identification from an individual sample from the study "Single-cell sequencing deconvolutes cellular responses to exercise in human skeletal muscle". 

## Data and Methodology
scRNA-Sequence data was accessed through Gene Expression Omnibus (GEO) through accession number GSE214544. Data from an individual before exercising was used in this analysis to be consistent with analyzing baseline skeletal muscle. The data was analyzed in R using Seurat, SeuratDisk, Tidyverse and SingleR packages to meet the above mentioned objectives. The dataset consists of 21509 features and 15564 cells. The code has been made available as a part of this repository.

## Exploring Results
### Data Pre-processing
#### 1. Quality Control
The plot below shows the QC metrics for the filtered data. QC is performed to filter out low quality gene and cell counts which include:

- Cells with low features
- Cells with high percentage of mitochondrial transcript are indicative of dying or poor quality cells.

Therefore, these conditions were used to filter out poor quality gene and cell instances which were assessed though various metric visualizations.

<img width="481" alt="image" src="https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/23698df4-be42-425e-a6d0-7d72b88d9d3f">

The figure visualizes QC metrics, from the plot, a threshold for each metric is given below, the following instances were included in downstream analysis: 
- Cells with a percentage of mitochondrial reads lower than 10%
- Cells with features more than 200 and less than 4000

#### 2. Normalization 
Normalization of the data was performed to enable comparison between cells and make more interpretable visualizations.

#### 3. Identifying Highly Variable Features
Not all features expressed in the cell are significant in terms of biological signal. To include the features that significantly vary from the average expression levels across the cells, these were identified to include downstream and visualized.

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/5361df52-18bd-4e7c-8466-c277b67083df)

The plot shows the highly variable features in red with the top 10 features labelled. To name a few characteristic features, the second most variable gene TNNI2 encodes the fast-twitch skeletal muscle protein responsible for regulation of striated muscle contraction. TNNC2 gene is encodes another troponin protein complex with 3 subunits also regulating muscle contraction. The top 10 consist of the many features conastituting the various structural and functional muscle cell groups for example ANKRD2 is involved in muscle stress response. They also include proteins interacting with macrophages.

#### 4. Scaling Data
Before going into dimensionality reduction (the next step) we scale the data to account for the variation by normalizing the mean expression to be 0 (with data points clustered around it which makes for a better visualization as well).

#### 5. Dimensionality Reduction
Dimensionality reduction is an importnat part of scRNA-Seq analysis. Since we compare thousands of cells across thousands of genes, we have a highly dimensional dataset. We therefore use dimensionality reduction methods, here Principle Component Analysis (PCA), to reduce it to a lower dimensional space. Our main aim is to account for as much variation (biological signal) as possible while reducing dimensionality.

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/4fc83352-2b3f-4c9c-b9ad-a8b611041c52)

From our PCA plot (using PC1 and PC2), we see 6 major groups just by eye. In order to select the number of PCs for our downstream analysis, we visualize the standard deviation or percent variance explained by each PC. The choice of number of PCs to include is made to maximize the variance explained and minimize dimensions. For this analysis, 20 PCs were included accounting.

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/68d20c0d-16a3-44cf-9af8-a1e4a15d97c0)

#### 5. Clustering
Clustering of cells based on expression levels was performed and visualized using PCA and UMAP respectively.

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/34ca86df-35a7-44fd-96c5-42e82bc7e9e9)

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/e72422c4-d7f0-40d9-8fa6-60e60c223780)

#### 6. Marker Identification
After clustering, to provide information about each cluster, we identify markers (expressed genes) that are charcteristic of a cluster. Markers for the cluster were identified and visualized with the top most expressed gene in each cluster as visualized below. 

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/4e89c43a-1680-41ce-bede-90f0b5237086)

Individual cluster marker identification was also carried out since this was an exploratory project to unravel the diversity in the skeletal muscle cells. The plot below shows Decorin gene (DCN) as the top expressed in cluster 0. DCN encodes one of the proteins from the leucine-rich family of glycoproteins. 

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/baeca970-21eb-4298-bcaa-70d22217e167)

Cluster 1 has the ACTA2 gene differentially expressed, as shown below. The ACTA2 gene or the Actin Alpha 2, Smooth Muscle gene encodes a actin protein that is highly conserved in eukaryotes and is responsible for cell biological functioning like cell structure and mobility. Mutations in this gene is linked with many diseases including coronary artery disease and stroke.

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/949c24c4-6638-4fd9-bf8c-f7ce19add1e3)

Another characteristic cluster was the cluster 4 visualized in the two plot below. It is the Apolipoprotein E gene encodes a major apoprotein involved in lipoprotein catabolism. Mutations in this gene are linked to Alzheimer's disease and cardiovascular disease.

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/803ed806-3d4d-4053-b527-837dde4c91fd)

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/ad0983d1-11f4-47c9-9eb2-742ff05e0b6d)

### Manual Cell Type Annotation
Due to the unavailability of a single-cell reference databse for skeletal muscles, manual cell type annotation was carried out through gene overlay plots and marker identification based on gene expression in cluster. This was then compared with and manually labelled using The Human Protein Atlas as the reference database. 

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/3d3f002c-1ab0-4045-ad3d-fbdb691ac4e4)


### Conclusion
Gene expression information for specific cell-types are important in understanding it in terms of the diversity of cell types that exist in samples, here skeletal muscles. This further helps us in the understanding and diagnosis of diseases and answer important molecular biology related questions.























