#####
## Apply signatures to Bassez et al 2021 BC dataset
#####

setwd('~/Data/scRNASeq/Bassez2021/')

# Load libraries
library(dplyr)
library(ggplot2)
library(anndata)

# Load expression and metadata and subset for cancer cells
load('exprData_Bassez2021.Rdata')

meta.bassez <- read.csv('1872-BIOKEY_metaData_cohort1_web.csv', header = TRUE)
meta.bassez <- meta.bassez[match(colnames(expr.data_bassez2021), meta.bassez$Cell), ]
expr.cancer <- expr.data_bassez2021[,meta.bassez$cellType == 'Cancer_cell']

# Load signature and subset expression data
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid.hrd <- signature.centroid.list$ElasticNet_alpha0.25

genes.intersect <- intersect(rownames(centroid.hrd), rownames(expr.cancer))
centroid.hrd <- centroid.hrd[genes.intersect, ]
expr.cancer <- expr.cancer[genes.intersect, ]

nonZero_genes <- apply(expr.cancer, 2, function(x) sum(x>0))
hist(nonZero_genes, breaks = 50)

# Calculate HRD scores across cells
hrdScores_bassez2021 <- data.frame(
  Cell = colnames(expr.cancer),
  Sample = sapply(colnames(expr.cancer), function(x) strsplit(x,split='_')[[1]][1]),
  HRD = apply(expr.cancer, 2, function(x) cor(x,centroid.hrd$HRD)),
  HR_proficient = apply(expr.cancer, 2, function(x) cor(x,centroid.hrd$HR_proficient))
)
hrdScores_bassez2021 <- hrdScores_bassez2021[!is.na(hrdScores_bassez2021$HRD), ]
hrdScores_bassez2021$HRD_score <- hrdScores_bassez2021$HRD - hrdScores_bassez2021$HR_proficient

# Prep for CellphoneDB analysis
cpdb.meta <- data.frame(
  Cell = meta.bassez$Cell,
  CellType = meta.bassez$cellType
)
hrdScores_bassez2021$HRD_group <- ifelse(hrdScores_bassez2021$HRD_score > 0, 'HRD', 'HR-proficient')
cpdb.meta <- merge(x = cpdb.meta, y = hrdScores_bassez2021[,c('Cell','HRD_group')], all.x = TRUE)
cpdb.meta$HRD_group[is.na(cpdb.meta$HRD_group)] <- ''
cpdb.meta$cell_type <- apply(cpdb.meta, 1, function(x) paste0(x[2],x[3],collapse = '_'))
cpdb.meta <- cpdb.meta[cpdb.meta$cell_type != 'Cancer_cell', ]

write.table(cpdb.meta[,c('Cell','cell_type')], file = '~/Projects/HRD_TranscriptionalSignature/CellphoneDB/Bassez2021/Bassez2021_meta.txt',
            row.names = FALSE, quote = FALSE, sep = '\t')

cpdb.expr <- as.data.frame(expr.data_bassez2021)
cpdb.expr <- cpdb.expr[,cpdb.meta$Cell]
cpdb.expr$Gene <- rownames(cpdb.expr)
cpdb.expr <- cpdb.expr[,c(ncol(cpdb.expr),1:(ncol(cpdb.expr)-1))]

write.table(cpdb.expr, file = '~/Projects/HRD_TranscriptionalSignature/CellphoneDB/Bassez2021/Bassez2021_counts.txt',
            row.names = FALSE, quote = FALSE, sep = '\t')

