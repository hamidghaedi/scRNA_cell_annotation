# scRNA cell annotation




In this repo, we will be doing:

(i) Assigning cell labels from reference data: or Automatic cell annotation for bladder cancer dataset used in [scRNA-seq analysis](https://github.com/hamidghaedi/scRNA_seq-analysis) using the reference data from [Tabula Sapiens](https://www.science.org/doi/full/10.1126/science.abl4896?af=R). 

(ii) Assigning cell labels to sub-cluster of epithelial cells from gene sets


# (i) Assigning cell labels using a ref dataset

```r 
# loading libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)
```

## Download Tabula Sapiens for bladder 

```r
# Set the destination directory for saving the downloaded file
destination_dir <- getwd()

# Set the URL of the file to be downloaded
file_url <- "https://ndownloader.figshare.com/files/27388874"

# Define the file name for saving the downloaded file
file_name <- "TS_Bladder.h5ad"  

# Download the file and save it to the specified directory
#download.file(file_url, file.path(destination_dir, file_name))
```

## Reading and processing ref data
in python:
```python 
# import libs
import scanpy as sc
from scipy import io
import os

# Set working directory
os.chdir("C://Users/qaedi/OneDrive - Queen's University/Documents/scRNA_cell_annotation")

# creating a dir
!mkdir matrix_files

# reading h5ad file
adata = sc.read_h5ad('TS_Bladder.h5ad')

adata
#AnnData object with n_obs × n_vars = 21568 × 58833
#    obs: 'Annotation', 'Predictability', 'Manually Annotated', 'Donor', 'Method', 'Organ', #'Compartment', 'Anatomical Information'
#    var: 'gene_symbol', 'ensembl_id', 'gene_length'
#    uns: 'Annotation_colors', 'Compartment_colors', 'Donor_colors', 'Manually Annotated_colors', #'Method_colors', 'Organ_colors', 'Propagated.Annotationcollapsed_colors', '_scvi', 'donor_colors', #'leiden', 'method_colors', 'neighbors', 'tissue_colors', 'umap'
#    obsm: 'X_umap'
#    layers: 'counts', 'raw_counts'
#    obsp: 'connectivities', 'distances'

# Ensuring that the raw count is the main matrix
adata.X = adata.layers['raw_counts']

# generating files needed to create Seurat object

with open('matrix_files/barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
        f.write(item + '\n')
with open('matrix_files/features.tsv', 'w') as f:
    for item in ['\t'.join([x,x,'Gene Expression']) for x in adata.var_names]:
        f.write(item + '\n')
io.mmwrite('matrix_files/matrix', adata.X.T)

# gzipping the files
import gzip
import glob

file_list = glob.glob("matrix_files/*")
for file_path in file_list:
    with open(file_path, "rb") as file_in:
        with gzip.open(file_path + ".gz", "wb") as file_out:
            file_out.write(file_in.read())
            
```

Now we have a directory with all needed files to create a Seurat object:

```r 
library(Seurat)
library(SingleR)
library(patchwork)
library(cowplot)


# create seurat object



## locating files
raw_data <- Read10X(data.dir = "matrix_files/")

## reading metadata
metadata <- read.csv("ref_metadata.csv")
rownames(metadata) <- metadata$X

## creating Seurat object
su <- CreateSeuratObject(counts = raw_data, meta.data = metadata)

## retain cells that are manually annotated
su <- subset(su, subset = Manually.Annotated == TRUE)

## creating a group for batch id
su$batch_id <- paste0(su$Donor, "_", su$Method)

#____________________Exploring source of variation______________#
## Normalizing the counts
su <- NormalizeData(object = su, normalization.method = "LogNormalize")

# Find variable feature
su <- FindVariableFeatures(su, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = TRUE)
                     
## Scaling
su <- ScaleData(su)

## Perform PCA
su <- RunPCA(su)

## Plot the PCA colored by cell cycle phase
no_split <- DimPlot(su,
        reduction = "pca",
        group.by= "batch_id")
        
with_split <- DimPlot(su,
        reduction = "pca",
        group.by= "batch_id",
        split.by= "batch_id")
        
no_split + with_split
```
![PCA_batch_id.png](https://github.com/hamidghaedi/scRNA_cell_annotation/blob/main/image/PCA_batch_id.png)

So there is no difference between batches , so we can proceed with file as is. 

```r
## Converting Seurat to sce object
refSce <- as.SingleCellExperiment(su)
```
## Reading and processing query data

```r
hsu <- readRDS("~/scRNA/github/harmonized_seurat.RDS")

# Adding cluster data
hsu$clusters <- Idents(hsu)

query_su <- CreateSeuratObject(counts = GetAssayData(object = hsu, assay = "RNA"), meta.data = hsu@meta.data)

# convert to sce
querySce <- as.SingleCellExperiment(query_su)
```

## Running `singleR` to predict cell types

```r
pred <- SingleR(test=querySce, ref=refSce, 
    labels=refSce$Annotation, de.method="wilcox")

# Checking rownames
all(rownames(query_su@meta.data) == rownames(data.frame(pred)))

# Adding predicted classes to metadata
query_su$singleR_pred <- data.frame(pred)$pruned.labels
```
### Prediction result visualization
```r
library(pheatmap)
tab <- table(query_su$clusters, query_su$singleR_pred)

pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))
```
As it showed in the below plot, the overlap between Tabula Sapiens labels and our manual labels is significant:

![heatmap_manual_labels_TS_labels.png](https://github.com/hamidghaedi/scRNA_cell_annotation/blob/main/image/heatmap_manual_labels_TS_labels.png)

