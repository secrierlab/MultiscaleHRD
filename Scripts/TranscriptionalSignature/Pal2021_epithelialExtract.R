#####
## Extraction of EPCAM+ clusters from the Pal et al. 2021 BC cohort
#####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(Matrix)
library(Seurat)
library(dplyr)
library(ggplot2)

setwd('~/Downloads/GSE161529_RAW/')

# Save matrix and barcode file names
files.matrix <- list.files(path = '~/Downloads/GSE161529_RAW', pattern = 'matrix', full.names = TRUE)
files.barcodes <- list.files(path = '~/Downloads/GSE161529_RAW', pattern = 'barcodes', full.name = TRUE)

# Remove normal samples and save tumour sample names
files.matrix <- files.matrix[!grepl(pattern = '_N', files.matrix)]
files.barcodes <- files.barcodes[!grepl(pattern = '_N', files.barcodes)]

samples.name <- sapply(files.matrix, function(x) 
  strsplit(strsplit(x, split='-')[[1]][1], split = '/')[[1]][6])

# features = gene names
features <- read.delim('~/Downloads/GSE161529_features.tsv.gz',
                       header = FALSE, sep ='\t')

# Initialise list to save EPCAM+ clusters:
#   This script is designed to run directly on a laptop
#   Consequently, it can take a while and so the data has been split into five
#   For i in x:y:
#     We extract the matrix and barcode files for sample i
#     This is followed by standard scRNA-seq preprocessing steps
#     Each sample is clustered and EPCAM expression across each cluster is checked
#     Clusters where  median(EPCAM) != 0 are labelled EPCAM+ and saved
#       These data are then added to cells.EPCAMpos
#     If there are no EPCAM+ clusters, we add NULL to the list
#     (samples 12 and 18 were removed because they crashed my computer)
cells.EPCAMpos <- list()

# for (i in 1:10) {
# for (i in c(11,13:17, 19:20)) {
# for (i in 21:30) {
# for (i in 31:40) {
for (i in 41:length(files.matrix)) {
  
  print(paste0('Analysing sample ', i, ' of ', length(files.matrix), ': ', Sys.time()))
  
  # Sort matrix files
  mtx <- readMM(files.matrix[i])
  mtx <- as.matrix(mtx)
  
  barcodes <- read.delim(files.barcodes[i], header = FALSE, sep='\t')
  
  rownames(mtx) <- features$V2
  colnames(mtx) <- barcodes$V1
  
  # Convert to Seurat object
  seurat.mtx <- CreateSeuratObject(counts = mtx[!duplicated(rownames(mtx)), ],
                                   meta.data = colnames(mtx[!duplicated(rownames(mtx)), ]) %>% data.frame(),
                                   min.cells = 3, min.features = 200,
                                   project = samples.name[i])
  
  # QC and selecting cells for further analysis
  seurat.mtx[['percent.mt']] <- PercentageFeatureSet(seurat.mtx, pattern='^MT-')
  
  # Normalizing the data
  seurat.mtx <- NormalizeData(seurat.mtx, normalization.method = 'LogNormalize', scale.factor = 10000)
  
  # Identification of highly variable features
  seurat.mtx <- FindVariableFeatures(seurat.mtx, selection.method = 'vst',
                                     nfeatures = 2000)
  
  # Scaling the data
  all.genes <- rownames(seurat.mtx)
  seurat.mtx <- ScaleData(seurat.mtx, features = all.genes)
  
  # Perform linear dimensional reduction
  seurat.mtx <- RunPCA(seurat.mtx, features = VariableFeatures(object = seurat.mtx))
  # DimPlot(seurat.mtx, reduction = 'pca')
  # DimHeatmap(seurat.mtx, dims = 1, cells = 500, balanced = TRUE)
  
  # Cluster the cells
  seurat.mtx <- FindNeighbors(seurat.mtx, dims = 1:10)
  seurat.mtx <- FindClusters(seurat.mtx, resolution = 0.5)
  
  # Run non-linear dimensional reduction
  seurat.mtx <- RunUMAP(seurat.mtx, dims = 1:10)
  
  # Identify relevant clusters
  df.umap <- data.frame(
    cluster = Idents(seurat.mtx),
    EPCAM = as.numeric(seurat.mtx@assays$RNA@counts['EPCAM', ])
  )
  
  df.umap.median <- df.umap %>%
    group_by(cluster) %>% summarise(median_EPCAM = median(EPCAM))
  
  if (sum(df.umap.median$median_EPCAM) == 0) {
    cells.EPCAMpos[[samples.name[i]]] <- NULL
  } else {
    
    clusters.EPCAMpos <- df.umap.median$cluster[df.umap.median$median_EPCAM > 0]
    
    # Extract cells in EPCAM+ clusters and save in list
    seurat.mtx.EPCAMpos <- seurat.mtx[, Idents(seurat.mtx) %in% clusters.EPCAMpos]
    
    cells.EPCAMpos[[samples.name[i]]] <- seurat.mtx.EPCAMpos
    
  }
  
}
Sys.time()

# Removes samples with no EPCAM+ clusters
cells.notNULL <- Filter(Negate(is.null), cells.EPCAMpos)
samples.remaining <- names(cells.notNULL)

# Merge Seurat objects
cells.merged <- merge(cells.notNULL[[1]],
                      y = cells.notNULL[-1],
                      add.cell.ids = samples.remaining)
table(cells.merged$orig.ident)

# Save merged Seurat objects
#     NOTE: Ensure filename is changed to reflect which samples were analysed
save(cells.merged, file = '~/Documents/PhD/Data/scRNASeq/Pal2021_EPCAMmarking/Pal2021_EPCAM_41to45.Rdata')

