#####
## Application of transcriptional signature to single-cell and bulk sequencing from Chung et al. 2017
####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# Load in data, remove lymph node sequencing, and seperate into bulk and single cells
dat <- read.table('~/Documents/PhD/Data/scRNASeq/Chung2017/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt', header = TRUE)
dat <- dat[!duplicated(dat$gene_name), ]

bulk <- dat[, grepl(pattern = 'Pooled', names(dat))]
names(bulk) <- sapply(names(bulk), function(x) strsplit(x,split='_')[[1]][1])
rownames(bulk) <- dat$gene_name

singleCell <- dat[,18:ncol(dat)]
singleCell <- singleCell[,!grepl(pattern = 'Re', names(singleCell))]
rownames(singleCell) <- dat$gene_name

# Use sample metadata to extract tumour cells
info <- read.table('~/Documents/PhD/Data/scRNASeq/Chung2017/GSE75688_final_sample_information.txt', header = TRUE)
tumor_cells <- info$sample[info$type == 'SC' & info$index == 'Tumor']
singleCell <- singleCell[,tumor_cells]

## Load in HRD template and calculate sample- and cell-wise HRD scores
load('Results/TranscriptionalSignature/Templates/HRD/template_HRD_ElasticNet.Rdata')

# Extract intersecting genes (n=114)
genes.intersect <- intersect(rownames(bulk), rownames(template.hrd))

template.hrd <- template.hrd[genes.intersect, ]
bulk <- bulk[genes.intersect, ]
singleCell <- singleCell[genes.intersect, ]

# Calculate HRD scores
scores_bulk <- data.frame(
  sample = colnames(bulk), 
  HRD = apply(bulk, 2, function(x) cor(x, template.hrd$HRD)),
  HR_prof = apply(bulk, 2, function(x) cor(x, template.hrd$HR_proficient))
)
scores_bulk$HRD_score_bulk <- scores_bulk$HRD - scores_bulk$HR_prof

scores_singleCell <- data.frame(
  sample = sapply(colnames(singleCell), function(x) strsplit(x,split='_')[[1]][1]),
  HRD = apply(singleCell, 2, function(x) cor(x, template.hrd$HRD)),
  HR_prof = apply(singleCell, 2, function(x) cor(x, template.hrd$HR_proficient))
)
scores_singleCell$HRD_score <- scores_singleCell$HRD - scores_singleCell$HR_prof

# Summarise single cells scores and match to bulk scores
scores_singleCell_mean <- scores_singleCell %>%
  group_by(sample) %>% summarise(mean_HRD_score = mean(HRD_score))

scores <- merge(x = scores_bulk[,c('sample','HRD_score_bulk')],
                y = scores_singleCell_mean)

g_chung2017 <- ggplot(scores, aes(x = HRD_score_bulk, y = mean_HRD_score)) +
  geom_point() + geom_smooth(method = 'lm') + stat_cor() +
  theme_minimal() +
  xlab('Bulk HRD') + ylab('Mean Single-cell HRD')
ggsave(filename = 'Figures/Figure7/Chung2017_singleCell_vs_Bulk.pdf',
       plot = g_chung2017, width = 4, height = 3)
