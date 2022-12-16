#####
## Application of transcriptional signature to single-cell sequencing from Qian et al. 2020
#####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(ggplot2)
library(dplyr)
library(wesanderson)
library(Seurat)

# Load data and subset for tumour cells
load('~/Documents/PhD/Data/scRNASeq/exprData_Qian2020.Rdata')
cell.info <- read.csv('~/Documents/PhD/Data/scRNASeq/2103-Breastcancer_metadata.csv')
tumor_cells <- cell.info$Cell[cell.info$CellType == 'Cancer']

expr.tumor <- expr.data_qian2020[, colnames(expr.data_qian2020) %in% tumor_cells]

# Load template and extract intersecting genes (n=113)
load('Results/TranscriptionalSignature/Templates/HRD/template_HRD_ElasticNet.Rdata')

genes.intersect <- intersect(rownames(template.hrd), rownames(expr.tumor))
template.hrd <- template.hrd[genes.intersect, ]
expr.tumor <- expr.tumor[genes.intersect, ]

# Check gene dropouts across cells and samples
prop_genes_per_cell <- apply(expr.tumor, 2, function(x) sum(x!=0)/length(x))
summary(prop_genes_per_cell)

expr.df <- data.frame(t(expr.tumor))
expr.df$sample <- sapply(rownames(expr.df), 
                         function(x) strsplit(x,split='_')[[1]][1])
expr.df <- expr.df %>% group_by(sample) %>%
  summarise_all(sum)
expr.df$GeneExpressed <- apply(expr.df[,-1], 1, 
                               function(x) sum(x!=0)/length(x))
g_include <- ggplot(expr.df, aes(x = sample, y = GeneExpressed)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank()) +
  ylab('% Genes Expressed') + ggtitle('Qian2020 Dropout Rates')
ggsave(filename = 'Figures/Supplementary/Qian2020_GeneDropout.pdf',
       plot = g_include, width = 4, height = 4)

# Calculate cell-wide HRD scores
qian2020_hrd <- data.frame(
  Sample = sapply(colnames(expr.tumor), function(x) strsplit(x,split='_')[[1]][1]),
  HRD = apply(expr.tumor, 2, function(x) cor(x, template.hrd$HRD)),
  HR_proficient = apply(expr.tumor, 2, function(x) cor(x, template.hrd$HR_proficient))
)
qian2020_hrd$HRD_score <- qian2020_hrd$HRD - qian2020_hrd$HR_proficient

# Summarise sample-wide HRD scores and plot with features
qian2020_hrd <- qian2020_hrd[!is.na(qian2020_hrd$HRD_score), ]
qian2020_hrd_props <- qian2020_hrd %>%
  group_by(Sample) %>% 
  summarise(prop_HRD = sum(HRD_score > 0)/n())
qian2020_hrd_props$BC_subtype = factor(c(
  'HER2', 'TN', 'B1_TN', 'TN', 'TN', 'HER2',
  'TN', 'HER2', 'Lum_HER2', 'TN', 'TN', 'TN', 'LumB', 'LumA'
), levels = c('B1_TN','Lum_HER2','HER2',
              'TN','LumA','LumB'))
qian2020_hrd_props <- qian2020_hrd_props[order(qian2020_hrd_props$prop_HRD), ]

qian2020_hrd_props$Sample <- factor(qian2020_hrd_props$Sample,
                                    levels = qian2020_hrd_props$Sample)

g_qian2020_HRDprops <- ggplot(data = qian2020_hrd_props, aes(x = Sample, y = prop_HRD, fill = BC_subtype)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab('% HRD cells') +
  scale_fill_manual(values = c(wes_palette('Moonrise3'),
                               wes_palette('Moonrise1')[1]))
ggsave(filename = 'Figures/Figure7/Qian2020_propHRDscores.pdf',
       plot = g_qian2020_HRDprops, width = 4, height = 4)

# Plot full distributions for the BRCA1-/- and StageII LumA samples
qian2020_hrd.plot <- qian2020_hrd[qian2020_hrd$Sample %in% c('sc5rJUQ033', 'sc5rJUQ064'),]
g_plotTwo <- ggplot(qian2020_hrd.plot, aes(x = HRD_score, fill = Sample)) +
  geom_density(alpha = 0.5) + theme_minimal() +
  theme(legend.position = 'top') +
  ylab('density(cells)') +
  scale_fill_manual(values = c(wes_palette('Moonrise3')[1],
                               wes_palette('Moonrise3')[5]))
ggsave(filename = 'Figures/Figure7/Qian2020_BRCA1_LumA_plots.pdf',
       plot = g_plotTwo, width = 4, height = 3)

# Plot full distributions
qian2020_hrd.plot2 <- merge(x = qian2020_hrd, y = qian2020_hrd_props[,-2])

qian2020_hrd.plot2$Sample <- factor(qian2020_hrd.plot2$Sample,
                                    levels = qian2020_hrd_props$Sample[14:1])

g_full <- ggplot(qian2020_hrd.plot2, aes(x = HRD_score, fill = BC_subtype)) +
  geom_density() + theme_minimal() + 
  scale_fill_manual(values = c(wes_palette('GrandBudapest1'),
                               wes_palette('GrandBudapest2')[1:2])) +
  facet_wrap( ~ Sample)
ggsave(filename = 'Figures/Supplementary/Qian2020_plotFull.pdf',
       plot = g_full, width = 7, height = 5)

## Qian2020 UMAP plotting

# Create new Seurat object with normalised data
qian_seurat <- CreateSeuratObject(counts = expr.data_qian2020[, colnames(expr.data_qian2020) %in% tumor_cells],
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
umap.data <- merge(x = umap.data, y = qian2020_hrd[,c('Sample', 'HRD_score')], by=0)
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
ggsave(filename = 'Figures/Figure7/Qian2020_UMAP_BC_subtype.pdf',
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
ggsave(filename = 'Figures/Figure7/Qian2020_UMAP_HRDscores.pdf',
       plot = g_qian2020_hrdScores, width = 4, height = 4)
