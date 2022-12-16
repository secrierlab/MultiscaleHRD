#####
## Preprocessing Qian et al 2020 BC dataset
#####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(Seurat)
library(dplyr)
library(cowplot)

# Load Qian et al. 2020 data
bc.data <- Read10X(data.dir = 
                     '~/Documents/PhD/Data/scRNASeq/export/BC_counts/')

# Initialize the Seurat object with the raw (non-normalized) data
#   init: 33694 genes, 44024 cells

## 8.3 Filtering low-quality cells
counts_per_cell <- Matrix::colSums(bc.data)
counts_per_gene <- Matrix::rowSums(bc.data)
genes_per_cell <- Matrix::colSums(bc.data>0)
cells_per_gene <- Matrix::rowSums(bc.data>0)

# 8.3.1 Summary counts for genes and cells
hist(log10(counts_per_cell+1), main='counts per cell', col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
hist(log10(counts_per_gene+1), main='counts per gene', col='wheat')
hist(log10(cells_per_gene+1), main='cells per gene', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat'); title('counts vs genes per cell')

# 8.3.2 Plot cells ranked by their number of detected genes
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')


## 8.4 Beginning with Seurat

# 8.4.1 Creating a seurat object

# Keep genes expressed in >= 3 cells (~.1% of the data).
# Keep all cells with at least 200 detected genes
#   now: 26040 genes, 38735 cells
seurat <- CreateSeuratObject(counts = bc.data, 
                             min.cells = 3, min.features = 200,
                             project = '10x_bc', assay = 'RNA')

## Load metadata
anno <- read.csv('~/Documents/PhD/Data/scRNASeq/2103-Breastcancer_metadata.csv', header = TRUE)
meta.data <- seurat@meta.data
all(rownames(meta.data) == anno$Cell)
seurat$CellFromTumor <- anno$CellFromTumor
seurat$PatientNumber <- anno$PatientNumber
seurat$CellType <- anno$CellType

# # Keep only cells from tumours
# seurat <- subset(x = seurat, subset = CellFromTumor == 'TRUE')

## 8.5 Preprocessing Step 1: Filter out low-quality cells (in addition to 1. 200 minimum features)

# Common metric for judging damaged cells: relative expression of mitochondrially derived genes
#   Tell-tale sign of cell stress: widespread RNA degradation after apoptosis

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = '^MT-', x = rownames(x = seurat@assays$RNA@data), value=TRUE)
percent.mito <- Matrix::colSums(seurat@assays$RNA@data[mito.genes, ])/Matrix::colSums(seurat@assays$RNA@data)
seurat <- AddMetaData(object = seurat, metadata = percent.mito, col.name = 'percent.mito')
VlnPlot(object = seurat, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mito'))

# Plot correlations RNA counts and other features
par(mfrow=c(1,2))
FeatureScatter(object = seurat, feature1 = 'nCount_RNA', feature2 = 'percent.mito')
FeatureScatter(object = seurat, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')

# Filter out cells with unique gene counts > 6,000 or mitochondrial content > 15%
seurat <- subset(x = seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                   percent.mito > -Inf & percent.mito < .15)

## 8.6.1 Preprocessing Step 2: Expression normalization
seurat <- NormalizeData(object = seurat, normalization.method = 'LogNormalize',
                        scale.factor = 10000)

# Save scRNA-seq data
expr.data_qian2020 <- as.matrix(GetAssayData(seurat, slot = 'data'))
save(expr.data_qian2020, file = '~/Documents/PhD/Data/scRNASeq/exprData_Qian2020.Rdata')

