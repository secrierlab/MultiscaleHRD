#####
## Cancer cell extraction, and application of transcriptional signature to single-cell sequencing from Pal et al. 2021
#####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

# setwd('~/Documents/PhD/Data/scRNASeq/Pal2021_EPCAMmarking/')

# Load EPCAM+ clusters and merge together Seurat objects
merged_cells.list <- list()
for (i in 1:length(list.files())) {
  print(paste0(i, ': ', Sys.time()))
  load(list.files(path = '~/Documents/PhD/Data/scRNASeq/Pal2021_EPCAMmarking/',
                  full.names = TRUE)[i])
  merged_cells.list[[i]] <- cells.merged
}

cells.merged <- merge(merged_cells.list[[1]],
                      y = merged_cells.list[-1])
rm(merged_cells.list)
table(cells.merged$orig.ident)
length(table(cells.merged$orig.ident))

# Identify highly variable features
cells.merged <- FindVariableFeatures(cells.merged, selection.method = 'vst',
                                     nfeatures = 2000)

# Scaling the data (only on variable genes due to computational limits)
# all.genes <- rownames(cells.merged)
cells.merged <- ScaleData(cells.merged)
# rm(all.genes)

# Perform linear dimensional reduction
cells.merged <- RunPCA(cells.merged, features = VariableFeatures(object = cells.merged))

# Cluster merged Seurat object
cells.merged <- FindNeighbors(cells.merged, dims = 1:10)
cells.merged <- FindClusters(cells.merged, resolution = 0.5)

# Run non-linear dimensional reduction
cells.merged <- RunUMAP(cells.merged, dims = 1:10)

# Plot UMAP
umap.data <- as.data.frame(cells.merged@reductions$umap@cell.embeddings)
umap.data$sample <- cells.merged$orig.ident
umap.data$cluster <- Idents(cells.merged)
umap.data$normal <- umap.data$cluster %in% c(2,14)

ggplot(umap.data, aes(x = UMAP_1, y = UMAP_2, color = sample)) +
  geom_point(size = .2) + theme_minimal()
ggplot(umap.data, aes(x = UMAP_1, y = UMAP_2, color = normal)) +
  geom_point(size = .2) + theme_minimal()

# Extract genes in HRD template
load('~/Documents/GitHub/HRD_classification/Results/TranscriptionalSignature/Templates/HRD/template_HRD_ElasticNet.Rdata')
genes.intersect <- intersect(rownames(cells.merged), rownames(template.hrd))

template.hrd <- template.hrd[genes.intersect, ]
cells.merged_small <- cells.merged[genes.intersect, !umap.data$normal]
cells.merged_small <- as.data.frame(cells.merged_small@assays$RNA@data)

# Check gene dropouts across cells and samples
prop_genes_per_cell <- apply(cells.merged_small, 1, function(x) sum(x!=0)/length(x))
summary(prop_genes_per_cell)

expr.df <- data.frame(t(cells.merged_small))
expr.df$sample <- sapply(rownames(expr.df), function(x)
  paste(strsplit(x, split='_')[[1]][1:2], collapse = '_'))
expr.df <- expr.df %>% group_by(sample) %>%
  summarise_all(sum)
expr.df$GeneExpressed <- apply(expr.df[,-1], 1, function(x)
  sum(x!=0)/length(x))
g_include <- ggplot(expr.df, aes(x = sample, y = GeneExpressed)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank()) +
  ylab('% Genes Expressed') + ggtitle('Pal2021 Dropout Rates')
ggsave(filename = 'Figures/Supplementary/Pal2021_GeneDropout.pdf',
       plot = g_include, width = 4, height = 4)

# Collate HRD scores
pal2021_res <- data.frame(
  Cell = colnames(cells.merged_small),
  Sample = sapply(colnames(cells.merged_small), function(x)
    paste(strsplit(x, split='_')[[1]][1:2], collapse = '_')),
  BC_subtype = sapply(colnames(cells.merged_small), function(x)
    strsplit(x, split='_')[[1]][2]),
  HRD = apply(cells.merged_small, 2, function(x) cor(x, template.hrd$HRD)),
  HR_proficient = apply(cells.merged_small, 2, function(x) cor(x, template.hrd$HR_proficient))
)
pal2021_res <- pal2021_res[!is.na(pal2021_res$HRD), ]
pal2021_res$HRD_score <- pal2021_res$HRD - pal2021_res$HR_proficient

# Summarise data
pal2021_summary <- pal2021_res %>%
  group_by(Sample, BC_subtype) %>%
  summarise(prop_HRDcells = sum(HRD_score > 0)/n(),
            cell_count = n()) %>%
  arrange(prop_HRDcells)
pal2021_summary <- pal2021_summary[pal2021_summary$cell_count > 20, ]
pal2021_summary$Sample <- factor(pal2021_summary$Sample, levels = pal2021_summary$Sample)
g_pal2021_waterfall <- ggplot(pal2021_summary, aes(x = Sample, y = prop_HRDcells, fill = BC_subtype)) + 
  geom_bar(stat = 'identity') + theme_minimal() + theme(axis.text.x = element_blank(),
                                                        axis.title.x = element_blank()) +
  scale_fill_manual(values = wes_palette('Moonrise3')) +
  ylab('% HRD cells')
ggsave(filename = 'Figures/Figure7/Pal2021_propHRDscores.pdf',
       plot = g_pal2021_waterfall)

pal2021_res$Sample <- factor(pal2021_res$Sample, 
                             levels = pal2021_summary$Sample[length(pal2021_summary$Sample):1])
pal2021_res <- pal2021_res[!is.na(pal2021_res$Sample), ]
g_pal2021_full <- ggplot(pal2021_res, aes(x = HRD_score, fill = BC_subtype)) +
  geom_density() + theme_minimal() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  facet_wrap(~Sample, ncol = 7)
ggsave(filename = 'Figures/Supplementary/Pal2021_plotFull.pdf',
       plot = g_pal2021_full, width = 9, height = 5)

# Lazy: plot with original UMAP coordinates
umap.scores <- merge(x = umap.data, y = pal2021_res[,c(3,6)], by=0)

g_pal2021_hrdScores <- ggplot(umap.scores, aes(x = UMAP_1, y = UMAP_2, color = HRD_score)) +
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
ggsave(filename = '~/Documents/GitHub/HRD_classification/Figures/Figure7/Pal2021_UMAP_HRDscores.pdf',
       plot = g_pal2021_hrdScores, width = 4, height = 4)


g_pal2021_BC_subtype <- ggplot(umap.scores, aes(x = UMAP_1, y = UMAP_2, color = BC_subtype)) +
  geom_point(size = .2) + theme_minimal() +
  scale_color_manual(values = wes_palette('Moonrise3')) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'top') +
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = '~/Documents/GitHub/HRD_classification/Figures/Figure7/Pal2021_UMAP_SampleType.pdf',
       plot = g_pal2021_BC_subtype, width = 4, height = 4)
