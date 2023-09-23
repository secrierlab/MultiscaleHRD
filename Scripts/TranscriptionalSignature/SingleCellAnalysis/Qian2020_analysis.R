#####
## Apply signatures to Qian et al 2020 BC dataset
#####

setwd('~/Data/scRNASeq/Qian2020/')

# Load libraries
library(dplyr)
library(ggplot2)
library(anndata)
library(wesanderson)

# Load expression and metadata and subset for cancer cells
load('exprData_Qian2020.Rdata')

meta.qian <- read.csv('2103-Breastcancer_metadata.csv', header = TRUE)
meta.qian <- meta.qian[match(colnames(expr.data_qian2020), meta.qian$Cell), ]
expr.cancer <- expr.data_qian2020[,meta.qian$CellType == 'Cancer']

# Load signature and subset expression data
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid.hrd <- signature.centroid.list$ElasticNet_alpha0.25

genes.intersect <- intersect(rownames(centroid.hrd), rownames(expr.cancer))
centroid.hrd <- centroid.hrd[genes.intersect, ]
expr.cancer <- expr.cancer[genes.intersect, ]

# nonZero_genes <- apply(expr.cancer, 2, function(x) sum(x>0))
# hist(nonZero_genes, breaks = 50)

# Calculate HRD scores across cells
hrdScores_qian2020 <- data.frame(
  Cell = colnames(expr.cancer),
  Sample = sapply(colnames(expr.cancer), function(x) strsplit(x,split='_')[[1]][1]),
  HRD = apply(expr.cancer, 2, function(x) cor(x,centroid.hrd$HRD)),
  HR_proficient = apply(expr.cancer, 2, function(x) cor(x,centroid.hrd$HR_proficient))
)
hrdScores_qian2020 <- hrdScores_qian2020[!is.na(hrdScores_qian2020$HRD), ]
hrdScores_qian2020$HRD_score <- hrdScores_qian2020$HRD - hrdScores_qian2020$HR_proficient

# Factor results in decreasing average HRD score
qian2020_hrdSummary <- hrdScores_qian2020 %>%
  group_by(Sample) %>% summarise(mean_HRD = mean(HRD_score)) %>%
  arrange(desc(mean_HRD))
hrdScores_qian2020$Sample <- factor(hrdScores_qian2020$Sample,
                                    levels = qian2020_hrdSummary$Sample)

# Plot density plots
ggplot(hrdScores_qian2020, aes(x = HRD_score)) +
  geom_density(fill = 'lightblue') + facet_wrap(~Sample)

# Compare BRCA1-/- vs Luminal A samples
hrdScores_plot <- hrdScores_qian2020
hrdScores_plot$label <- NA
hrdScores_plot$label[hrdScores_plot$Sample == 'sc5rJUQ033'] <- 'BRCA1-/- TNBC'
hrdScores_plot$label[hrdScores_plot$Sample == 'sc5rJUQ064'] <- 'Lum A-like'
hrdScores_plot <- hrdScores_plot[!is.na(hrdScores_plot$label), ]

g_plotTwo <- ggplot(hrdScores_plot, aes(x = HRD_score, fill = label)) +
  geom_density(alpha = 0.4) + theme_minimal() +
  theme(legend.position = 'top',
        legend.title = element_blank()) +
  ylab('density(cells)') +
  scale_fill_manual(values = c(wes_palette('Moonrise3')[c(1,5)]))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure6/Qian2020_BRCA1_LumA.pdf',
       plot = g_plotTwo, width = 4, height = 3)

# Prep for CellphoneDB analysis
cpdb.meta <- data.frame(
  Cell = meta.qian$Cell,
  CellType = meta.qian$CellType
)
hrdScores_qian2020$HRD_group <- ifelse(hrdScores_qian2020$HRD_score > 0, 'HRD', 'HR-proficient')
cpdb.meta <- merge(x = cpdb.meta, y = hrdScores_qian2020[,c('Cell','HRD_group')], all.x = TRUE)
cpdb.meta$HRD_group[is.na(cpdb.meta$HRD_group)] <- ''
cpdb.meta$cell_type <- apply(cpdb.meta, 1, function(x) paste0(x[2],x[3],collapse = '_'))
cpdb.meta <- cpdb.meta[cpdb.meta$cell_type != 'Cancer', ]

write.table(cpdb.meta[,c('Cell','cell_type')], file = '~/Projects/HRD_TranscriptionalSignature/CellphoneDB/Qian2020/Qian2020_meta.txt',
            row.names = FALSE, quote = FALSE, sep = '\t')

cpdb.expr <- as.data.frame(expr.data_qian2020)
cpdb.expr <- cpdb.expr[,cpdb.meta$Cell]
cpdb.expr$Gene <- rownames(cpdb.expr)
cpdb.expr <- cpdb.expr[,c(ncol(cpdb.expr),1:(ncol(cpdb.expr)-1))]

write.table(cpdb.expr, file = '~/Projects/HRD_TranscriptionalSignature/CellphoneDB/Qian2020/Qian2020_counts.txt',
            row.names = FALSE, quote = FALSE, sep = '\t')

