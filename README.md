# scRNA Sequence Analysis Exploring Heterogeneity in Human Skeletal Muscle
scRNA-Sequence Analysis was performed as a way to analyze the transcriptional profile of human skeletal muscles using Seurat package in R.

## Objective 
scRNA-Sequencing differs from Bulk RNA Sequencing in provided a finer resolution of gene expression within cells enabling us to explore  the hetergeneity between cells. From analysing baseline transcriptome of the Human skeletal muscles (https://github.com/maitree-patel/RNASeq-Exploratory-Data-Analysis-of-Human-Skeletal-Muscles), scRNA_Sequence analysis was further used to expllore the single-cell transcriptome. The main aim was to understand the heterogeneity in Skeletal muscle cells by loading, pre-processing, clustering and marker identification from one sample from the study "Single-cell sequencing deconvolutes cellular responses to exercise in human skeletal muscle". 

## Data and Methodology
scRNA-Sequence data was accessed through Gene Expression Omnibus (GEO) through accession number GSE214544. Data from an individual before exercising was used in this analysis to be consistent with analyzing baseline skeletal muscle. The data was analyzed in R using Seurat, SeuratDisk, Tidyverse and SingleR packages to meet the above mentioned objectives. The dataset consists of 21509 features and 15564 cells. The code has been made available as a part of this repository.

## Exploring Results
### Data Pre-processing
#### 1. Quality Control
The plot below shows the QC metrics for the filtered data. QC is performed to filter out low quality which include:

- Cells with low features
- Cells with high percentage of mitochondrial transcript are indicative of dying or poor quality cells.
Therefore, these conditions were used to filter out poor quality cells which were assessed though various metric visualizations.

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/e8cf4784-c8ac-4efa-8a0d-09d4b9e2e8b4)

The figure shows all cells with a percentage of mitochondrial reads lower than 10% and cells with features more than 200 and less than 4000 as included for downstream analysis. Fixed thresholds were applied by visualizing the QC metrics to apply to the data.

#### 2. Normalization 
Normalization of the data was performed to enable comparison between cells and make more interpretable visualizations.

#### 3. Identifying Highly Variable Features
Not all features expressed in the cell are significant in terms of biological signal. To include the features that significantly vary from the average expression levels across the cells, these were identified to include downstream and visualized.

![image](https://github.com/maitree-patel/scRNA-Sequence-Analysis-Exploring-Heterogeneity-in-Human-Skeletal-Muscle/assets/134908239/185dc2de-187d-452a-94b9-6003fb2e064f)

The plot shows the highly variable features in red with the top 10 features labelled. To name a few characteristic features, the second most variable gene TNNI2 encodes the fast-twitch skeletal muscle protein responsible for regulation of striated muscle contraction. TNNC2 gene is encodes another troponin protein complex with 3 subunits also regulating muscle contraction. The top 10 consist of the many features conastituting the various structural and functional muscle cell groups for example ANKRD2 is involved in muscle stress response. They also include proteins interacting with macrophages.











