#####
## Apply transcriptional signature to Chung et al. 2017 bulk and scRNAseq
#####

# Load libraries
library(dplyr)
library(ggpubr)

# Load signature
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid.hrd <- signature.centroid.list$ElasticNet_alpha0.25

# Load expression data, match with signature, and log2-normalise
expr.chung <- read.table('~/Data/scRNASeq/Chung2017/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt',h=T)

genes.intersect <- intersect(rownames(centroid.hrd), expr.chung$gene_name)

centroid.hrd <- centroid.hrd[genes.intersect, ]
expr.chung <- expr.chung[match(genes.intersect, expr.chung$gene_name), ]
rownames(expr.chung) <- expr.chung$gene_name
expr.chung <- expr.chung[,-c(1:3)]

expr.chung <- log2(expr.chung + 1)

# Separate into bulk and scRNAseq
expr.chung_bulk <- expr.chung[,grepl(pattern = 'Pooled', colnames(expr.chung))]
expr.chung_sc <- expr.chung[,15:ncol(expr.chung)]

# Extract only tumour cells from single cell data
chung_info <- read.table('~/Data/scRNASeq/Chung2017/GSE75688_final_sample_information.txt',h=T)
cells.tumor <- chung_info$sample[chung_info$type == 'SC' & chung_info$index == 'Tumor']

expr.chung_sc <- expr.chung_sc[,cells.tumor]

# Calculate HRD scores
hrd_scores.bk <- data.frame(
  sample = sapply(colnames(expr.chung_bulk), function(x) strsplit(x,split='_')[[1]][1]),
  HRD = apply(expr.chung_bulk, 2, function(x) cor(x,centroid.hrd$HRD)),
  HR_proficient = apply(expr.chung_bulk, 2, function(x) cor(x,centroid.hrd$HR_proficient))
)
hrd_scores.bk$HRD_score_bulk <- hrd_scores.bk$HRD - hrd_scores.bk$HR_proficient

hrd_scores.sc <- data.frame(
  cell = colnames(expr.chung_sc),
  HRD = apply(expr.chung_sc, 2, function(x) cor(x,centroid.hrd$HRD)),
  HR_proficient = apply(expr.chung_sc, 2, function(x) cor(x,centroid.hrd$HR_proficient))
)
hrd_scores.sc$HRD_score_sc <- hrd_scores.sc$HRD - hrd_scores.sc$HR_proficient
hrd_scores.sc$sample <- sapply(hrd_scores.sc$cell, function(x) strsplit(x,split='_')[[1]][1])

# Match bulk and single-cell HRD scores
hrd_scores.sc_summary <- hrd_scores.sc %>%
  group_by(sample) %>% summarise(mean_HRD_score_sc = mean(HRD_score_sc))

hrd_scores.df <- merge(x = hrd_scores.bk[,c('sample','HRD_score_bulk')],
                       y = hrd_scores.sc_summary)

# Plot results
g_chung <- ggplot(hrd_scores.df, aes(x = HRD_score_bulk, y = mean_HRD_score_sc)) +
  geom_point() + geom_smooth(method = 'lm', col = 'darkred') + stat_cor() +
  theme_minimal() +
  xlab('Bulk HRD') + ylab('Mean single-cell HRD')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure6/Chung_BulkSingleCell.pdf',
       plot = g_chung, width = 4, height = 3)
