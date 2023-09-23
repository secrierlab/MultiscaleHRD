#####
## Apply signatures to Qian et al 2020 BC dataset
#####

setwd('~/Data/scRNASeq/Qian2020/')

# Load libraries
library(dplyr)
library(ggplot2)
library(anndata)
library(wesanderson)
library(Seurat)

# Load expression and metadata and subset for cancer cells
load('exprData_Qian2020.Rdata')

meta.qian <- read.csv('2103-Breastcancer_metadata.csv', header = TRUE)
meta.qian <- meta.qian[match(colnames(expr.data_qian2020), meta.qian$Cell), ]

# Load signature and subset expression data
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid.hrd <- signature.centroid.list$ElasticNet_alpha0.25
expr.hrd <- expr.data_qian2020[rownames(centroid.hrd), meta.qian$Cell]

# Gene inclusion plotting
expr.nonZero <- data.frame(
  Cell = colnames(expr.hrd), 
  Sample = sapply(colnames(expr.hrd), function(x) strsplit(x,split='_')[[1]][1]),
  prop_GenesExpressed = apply(expr.hrd, 2, function(x) 100*mean(x>0))
)
g_nonZero <- ggplot(expr.nonZero, aes(x = prop_GenesExpressed, fill = Sample)) +
  geom_density(alpha = 0.4) +
  theme_minimal() + theme(legend.position = 'none') +
  geom_vline(xintercept = mean(expr.nonZero$prop_GenesExpressed), col = 'red', linetype = 'dashed') +
  xlab('% Genes Expressed / Cell')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_QianNonZeroGenes.pdf',
       plot = g_nonZero)

expr.nonZero_summary <- expr.nonZero %>%
  group_by(Sample) %>% summarise(meanExpression = mean(prop_GenesExpressed))
g_nonZeroSummary <- ggplot(expr.nonZero_summary, aes(x = meanExpression)) +
  geom_histogram(fill = 'lightblue', color = 'darkblue') +
  theme_minimal() + xlab('Mean % Genes Expressed Across Cells / Sample')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_QianNonZeroGeneSummary.pdf',
       plot = g_nonZeroSummary)

# Calculate HRD scores across cells
hrdScores_qian2020 <- data.frame(
  Cell = colnames(expr.hrd),
  CellType = meta.qian$CellType,
  Sample = sapply(colnames(expr.hrd), function(x) strsplit(x,split='_')[[1]][1]),
  HRD = apply(expr.hrd, 2, function(x) cor(x,centroid.hrd$HRD)),
  HR_proficient = apply(expr.hrd, 2, function(x) cor(x,centroid.hrd$HR_proficient))
)
hrdScores_qian2020 <- hrdScores_qian2020[!is.na(hrdScores_qian2020$HRD), ]
hrdScores_qian2020$HRD_score <- hrdScores_qian2020$HRD - hrdScores_qian2020$HR_proficient

# Factor results in decreasing average HRD score
qian2020_hrdSummary <- hrdScores_qian2020[hrdScores_qian2020$CellType == 'Cancer',] %>%
  group_by(Sample) %>% summarise(mean_HRD = mean(HRD_score)) %>%
  arrange(desc(mean_HRD))
hrdScores_qian2020$Sample <- factor(hrdScores_qian2020$Sample,
                                    levels = qian2020_hrdSummary$Sample)
hrdScores_qian2020$Cancer <- hrdScores_qian2020$CellType == 'Cancer'

mu <- hrdScores_qian2020 %>%
  group_by(Sample, Cancer) %>% summarise(medianHRD = median(HRD_score))

g_tme <-ggplot(mu, aes(Cancer, medianHRD, fill=Cancer)) +
  geom_boxplot() +
  geom_point() + geom_line(aes(group = Sample)) +
  theme_minimal() +
  theme(legend.position = 'none') + ylab('median HRD score') +
  scale_fill_manual(values = wes_palette('Darjeeling2')) +
  ggtitle('Qian et al.')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_HRDinTME_Qian.pdf',
       plot = g_tme, width = 4, height = 4)

# Map clinical data and proportion of cells with HRD > 0
qian2020_hrd_props <- hrdScores_qian2020[hrdScores_qian2020$Cancer,] %>%
  group_by(as.character(Sample)) %>% summarise(prop_HRD = 100*mean(HRD_score > 0))
qian2020_hrd_props$BC_subtype = factor(c(
  'HER2', 'TN', 'B1_TN', 'TN', 'TN', 'HER2',
  'TN', 'HER2', 'Lum_HER2', 'TN', 'TN', 'TN', 'LumB', 'LumA'
), levels = c('B1_TN','Lum_HER2','HER2',
              'TN','LumA','LumB'))
qian2020_hrd_props <- qian2020_hrd_props[order(qian2020_hrd_props$prop_HRD), ]

names(qian2020_hrd_props)[1] <- 'Sample'
qian2020_hrd_props$Sample <- factor(qian2020_hrd_props$Sample,
                                    levels = qian2020_hrd_props$Sample)

g_qian2020_HRDprops <- ggplot(data = qian2020_hrd_props, aes(x = Sample, y = prop_HRD, fill = BC_subtype)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab('% HRD cells') +
  scale_fill_manual(values = c(wes_palette('Moonrise3'),
                               wes_palette('Moonrise1')[1]))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure6/Qian2020_propHRDscores.pdf',
       plot = g_qian2020_HRDprops, width = 4, height = 4)

# Plot density plots
hrdScores_qian2020 <- merge(x = hrdScores_qian2020, y = qian2020_hrd_props[,c('Sample','BC_subtype')])
g_densities <- ggplot(hrdScores_qian2020[hrdScores_qian2020$Cancer, ], aes(x = HRD_score, fill = BC_subtype)) +
  geom_density(alpha = 0.5) + 
  scale_fill_manual(values = c(wes_palette('Moonrise3'),
                               wes_palette('Moonrise1')[1])) +
  theme(legend.position = 'top') +
  facet_wrap(~Sample, ncol = 4)
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_HRDdensities_Qian.pdf',
       plot = g_densities)

## Qian2020 UMAP plotting

# Create new Seurat object with normalised data
qian_seurat <- CreateSeuratObject(counts = expr.data_qian2020[, meta.qian$Cell[meta.qian$CellType == 'Cancer']],
                                  project = 'qian2020', min.cells = 3, min.features = 200)

# Identify highly variable features and scale data
qian_seurat <- FindVariableFeatures(qian_seurat, selection.method = 'vst',
                                    nfeatures = 2000)
qian_seurat <- ScaleData(qian_seurat)

# Perform linear dimensionality reduction, clustering, and UMAP
qian_seurat <- RunPCA(qian_seurat, features = VariableFeatures(object = qian_seurat))
qian_seurat <- FindNeighbors(qian_seurat, dims = 1:10)
qian_seurat <- FindClusters(qian_seurat, resolution = 0.5)
qian_seurat <- RunUMAP(qian_seurat, dims = 1:10)

# UMAP plotting
umap.data <- as.data.frame(qian_seurat@reductions$umap@cell.embeddings)
umap.data <- merge(x = umap.data, y = hrdScores_qian2020[,c('Sample', 'HRD_score')], by=0)
umap.data <- merge(x = umap.data, y = qian2020_hrd_props[,c('Sample','BC_subtype')], by = 'Sample')

# Plot UMAP coordinates coloured by breast cancer subtype
g_qian2020_BC_subtype <- ggplot(umap.data, aes(x = UMAP_1, y = UMAP_2, color = BC_subtype)) +
  geom_point(size = .2) + theme_minimal() +
  scale_color_manual(values = c(wes_palette('Moonrise3'),
                                wes_palette('Moonrise1')[1])) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'top') +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure7/Qian2020_UMAP_BC_subtype.pdf',
       plot = g_qian2020_BC_subtype, width = 4, height = 4)

# Plot UMAP coordinates coloured by HRD score
g_qian2020_hrdScores <- ggplot(umap.data, aes(x = UMAP_1, y = UMAP_2, color = HRD_score)) +
  geom_point(size = .2) + theme_minimal() +
  scale_color_gradientn(colors = c(rep(wes_palette('GrandBudapest1')[1],2),
                                   'gray90',
                                   rep(wes_palette('GrandBudapest1')[2],2)),
                        values = c(0,.3,
                                   abs(min(umap.data$HRD_score))/(abs(min(umap.data$HRD_score)) + max(umap.data$HRD_score)),
                                   .7,1)) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'top')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure7/Qian2020_UMAP_HRDscores.pdf',
       plot = g_qian2020_hrdScores, width = 4, height = 4)
