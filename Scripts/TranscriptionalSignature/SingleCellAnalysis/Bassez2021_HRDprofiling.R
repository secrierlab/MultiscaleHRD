#####
## Apply signatures to bassez et al 2021 BC dataset
#####

setwd('~/Data/scRNASeq/Bassez2021/')

# Load libraries
library(dplyr)
library(ggplot2)
library(anndata)
library(wesanderson)
library(Seurat)

# Load expression and metadata and subset for cancer cells
load('exprData_Bassez2021.Rdata')

meta.bassez <- read.csv('1872-BIOKEY_metaData_cohort1_web.csv', header = TRUE)
meta.bassez <- meta.bassez[match(colnames(expr.data_bassez2021), meta.bassez$Cell), ]

# Load signature and subset expression data
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid.hrd <- signature.centroid.list$ElasticNet_alpha0.25
expr.hrd <- expr.data_bassez2021[rownames(centroid.hrd), meta.bassez$Cell]

# Gene inclusion plotting
expr.nonZero <- data.frame(
  Cell = colnames(expr.hrd), 
  Sample = sapply(colnames(expr.hrd), function(x) paste0(strsplit(x,split='_')[[1]][1:2], collapse='_')),
  prop_GenesExpressed = apply(expr.hrd, 2, function(x) 100*mean(x>0))
)
g_nonZero <- ggplot(expr.nonZero, aes(x = prop_GenesExpressed, fill = Sample)) +
  geom_density(alpha = 0.4) +
  theme_minimal() + theme(legend.position = 'none') +
  geom_vline(xintercept = mean(expr.nonZero$prop_GenesExpressed), col = 'red', linetype = 'dashed') +
  xlab('% Genes Expressed / Cell')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_BassezNonZeroGenes.pdf',
       plot = g_nonZero, width = 5.35, height = 4.79)

expr.nonZero_summary <- expr.nonZero %>%
  group_by(Sample) %>% summarise(meanExpression = mean(prop_GenesExpressed))
g_nonZeroSummary <- ggplot(expr.nonZero_summary, aes(x = meanExpression)) +
  geom_histogram(fill = 'lightblue', color = 'darkblue') +
  theme_minimal() + xlab('Mean % Genes Expressed Across Cells / Sample')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_BassezNonZeroGeneSummary.pdf',
       plot = g_nonZeroSummary, width = 5.35, height = 4.79)

# Calculate HRD scores across cells
hrdScores_bassez2021 <- data.frame(
  Cell = colnames(expr.hrd),
  CellType = meta.bassez$cellType,
  Sample = sapply(colnames(expr.hrd), function(x) paste0(strsplit(x,split='_')[[1]][1:2],collapse = '_')),
  HRD = apply(expr.hrd, 2, function(x) cor(x,centroid.hrd$HRD)),
  HR_proficient = apply(expr.hrd, 2, function(x) cor(x,centroid.hrd$HR_proficient))
)
hrdScores_bassez2021 <- hrdScores_bassez2021[!is.na(hrdScores_bassez2021$HRD), ]
hrdScores_bassez2021$HRD_score <- hrdScores_bassez2021$HRD - hrdScores_bassez2021$HR_proficient

# Factor results in decreasing average HRD score
bassez2021_hrdSummary <- hrdScores_bassez2021[hrdScores_bassez2021$CellType == 'Cancer_cell',] %>%
  group_by(Sample) %>% summarise(mean_HRD = mean(HRD_score)) %>%
  arrange(desc(mean_HRD))
hrdScores_bassez2021$Sample <- factor(hrdScores_bassez2021$Sample,
                                    levels = bassez2021_hrdSummary$Sample)
hrdScores_bassez2021$Cancer <- hrdScores_bassez2021$CellType == 'Cancer_cell'

mu <- hrdScores_bassez2021 %>%
  group_by(Sample, Cancer) %>% summarise(medianHRD = median(HRD_score))

g_tme <-ggplot(mu, aes(Cancer, medianHRD, fill=Cancer)) +
  geom_boxplot() +
  geom_point() + geom_line(aes(group = Sample)) +
  theme_minimal() +
  theme(legend.position = 'none') + ylab('median HRD score') +
  scale_fill_manual(values = wes_palette('Darjeeling2')) +
  ggtitle('Bassez et al.')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_HRDinTME_Bassez.pdf',
       plot = g_tme, width = 4, height = 4)

# Map clinical data and proportion of cells with HRD > 0
bassez2021_hrd_props <- hrdScores_bassez2021[hrdScores_bassez2021$Cancer,] %>%
  group_by(as.character(Sample)) %>% summarise(prop_HRD = 100*mean(HRD_score > 0))
names(bassez2021_hrd_props)[1] <- 'Sample'
bassez2021_hrd_props$id <- sapply(bassez2021_hrd_props$Sample,
                                  function(x) as.numeric(strsplit(x,split='_')[[1]][2]))
bassez2021_hrd_props <- bassez2021_hrd_props[order(bassez2021_hrd_props$id), ]

bassez2021_hrd_props$BC_subtype = factor(c(
  'ER-HER2+','TNBC','TNBC','TNBC','TNBC',
  'ER-HER2+','TNBC','ER-HER2+','ER+HER2+','TNBC',
  'TNBC','ER+HER2-','ER+HER2-','ER+HER2-','TNBC',
  'TNBC','ER+HER2-','ER+HER2-','TNBC','ER+HER2-',
  'ER+HER2-','ER+HER2-','TNBC','ER+HER2-','TNBC',
  'ER+HER2+','TNBC','ER+HER2-','ER+HER2-','ER+HER2-','ER+HER2-'
), levels = c('ER+HER2+','ER+HER2-','ER-HER2+','TNBC'))
bassez2021_hrd_props <- bassez2021_hrd_props[order(bassez2021_hrd_props$prop_HRD), ]

bassez2021_hrd_props$Sample <- factor(bassez2021_hrd_props$Sample,
                                    levels = bassez2021_hrd_props$Sample)

g_bassez2021_HRDprops <- ggplot(data = bassez2021_hrd_props, aes(x = Sample, y = prop_HRD, fill = BC_subtype)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab('% HRD cells') +
  scale_fill_manual(values = c(wes_palette('Moonrise3')[c(5,2:4)]))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure6/Bassez2021_propHRDscores.pdf',
       plot = g_bassez2021_HRDprops, width = 4, height = 4)

# Plot density plots
hrdScores_bassez2021 <- merge(x = hrdScores_bassez2021, y = bassez2021_hrd_props[,c('Sample','BC_subtype')])
g_densities <- ggplot(hrdScores_bassez2021[hrdScores_bassez2021$Cancer, ], aes(x = HRD_score, fill = BC_subtype)) +
  geom_density(alpha = 0.5) + 
  scale_fill_manual(values = c(wes_palette('Moonrise3')[c(5,2:4)])) +
  theme(legend.position = 'top') +
  facet_wrap(~Sample, ncol = 4)
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_HRDdensities_Bassez.pdf',
       plot = g_densities, height = 8, width = 5.35)

## Bassez2021 UMAP plotting

# Create new Seurat object with normalised data
bassez_seurat <- CreateSeuratObject(counts = expr.data_bassez2021[, meta.bassez$Cell[meta.bassez$cellType == 'Cancer_cell']],
                                    project = 'bassez2021', min.cells = 3, min.features = 200)

# Identify highly variable features and scale data
bassez_seurat <- FindVariableFeatures(bassez_seurat, selection.method = 'vst',
                                    nfeatures = 2000)
bassez_seurat <- ScaleData(bassez_seurat)

# Perform linear dimensionality reduction, clustering, and UMAP
bassez_seurat <- RunPCA(bassez_seurat, features = VariableFeatures(object = bassez_seurat))
bassez_seurat <- FindNeighbors(bassez_seurat, dims = 1:10)
bassez_seurat <- FindClusters(bassez_seurat, resolution = 0.5)
bassez_seurat <- RunUMAP(bassez_seurat, dims = 1:10)

# UMAP plotting
umap.data <- as.data.frame(bassez_seurat@reductions$umap@cell.embeddings)
rownames(hrdScores_bassez2021) <- hrdScores_bassez2021$Cell
umap.data <- merge(x = umap.data, y = hrdScores_bassez2021[,c('Sample', 'HRD_score')], by=0)
umap.data <- merge(x = umap.data, y = bassez2021_hrd_props[,c('Sample','BC_subtype')], by = 'Sample')

# Plot UMAP coordinates coloured by breast cancer subtype
g_bassez2021_BC_subtype <- ggplot(umap.data, aes(x = UMAP_1, y = UMAP_2, color = BC_subtype)) +
  geom_point(size = .2) + theme_minimal() +
  scale_color_manual(values = c(wes_palette('Moonrise3')[c(5,2:4)])) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'top') +
  guides(color = guide_legend(override.aes = list(size = 5),
                              nrow = 2))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure6/bassez2021_UMAP_BC_subtype.pdf',
       plot = g_bassez2021_BC_subtype, width = 4, height = 4)

# Plot UMAP coordinates coloured by HRD score
g_bassez2021_hrdScores <- ggplot(umap.data, aes(x = UMAP_1, y = UMAP_2, color = HRD_score)) +
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
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure6/bassez2021_UMAP_HRDscores.pdf',
       plot = g_bassez2021_hrdScores, width = 4, height = 4)

