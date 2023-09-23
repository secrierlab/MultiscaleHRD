#####
## Plot Heatmap of TCGA-BRCA training cohort for chosen signature
#####

# Load libraries
library(wesanderson)
library(pheatmap)

# Load and process data
load('~/Projects/HRD_TranscriptionalSignature/Data/TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.Rdata')

samples <- sapply(rownames(Z.tumor_training), function(x) substr(x, 1, 12))
Z.tumor_training <- log2(Z.tumor_training + 1)
Z.tumor_training <- apply(Z.tumor_training, 2, scale)
rownames(Z.tumor_training) <- samples

load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata') # Group Annotation
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('BRCA_status', 'HRD')]
ann_tcga_train <- ann_tcga[rownames(Z.tumor_training), ]

# Load signature and subset expression data
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
sig <- signature.centroid.list$ElasticNet_alpha0.25

Z.tumor_training <- Z.tumor_training[,rownames(sig)]

ann_cols <- list(
  BRCA_status = c('BRCA1' = 'blue', 'BRCA2' = 'red', 'none' = 'white'),
  HRD = c('HRD' = wes_palette('GrandBudapest1')[2], 'HR-proficient' = wes_palette('GrandBudapest1')[1])
)

# Set colour range
paletteLength <- 100
myColor <- colorRampPalette(c('navy', 'darkblue', 'white', 'red', 'darkred'))(paletteLength)
myBreaks <- c(seq(min(Z.tumor_training), -3, length.out = ceiling(paletteLength/4)+1),
              seq(-3, 0, length.out = floor(paletteLength/4))[-1],
              seq(0, 3, length.out = floor(paletteLength/4))[-1],
              seq(3, max(Z.tumor_training), length.out = floor(paletteLength/4))[-1])

pheatmap(t(Z.tumor_training), show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = ann_tcga_train, annotation_colors = ann_cols,
         color = myColor, breaks = myBreaks, clustering_method = 'average'
         , filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_Heatmap.pdf'
         )
