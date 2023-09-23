setwd('~/Data/TCGA')

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(tibble)
library(decoupleR)
library(ggplot2)

query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)
# GDCdownload(query)
data <- GDCprepare(query = query)
data <- data[, data$sample_type == 'Primary Tumor']
data <- data[, !duplicated(data$patient)]

# Match with HRD classification
load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')
# ann_tcga$Patient <- rownames(ann_tcga)
ann_tcga$HRD <- sapply(ann_tcga$HRD_prob, function(x) ifelse(x>=0.79,'HRD','HR-proficient'))
ann_tcga <- ann_tcga[!is.na(ann_tcga$ER_status), ]
# ann_tcga <- ann_tcga[ann_tcga$ER_status == 'Negative', ]
ann_tcga <- ann_tcga[,c('Patient','HRD')]

patients.intersect <- intersect(ann_tcga$Patient, data$patient)
ann_tcga <- ann_tcga[match(patients.intersect, ann_tcga$Patient), ]
data <- data[, match(patients.intersect, data$patient)]

data$HRD_status <- factor(ann_tcga$HRD, levels = c('HR-proficient','HRD'))

# Construct a DESeqDataSet data object
dds <- DESeqDataSet(data, design = ~ HRD_status)

# Gene count filtering
dds <- dds[rowSums(counts(dds)) > 10, ]

# Normalisation
dds <- estimateSizeFactors(dds)

# DIFFERENTIAL GENE EXPRESSION ANALYSIS
dds_DGE <- DESeq(dds)
dds_DGE_results <- results(dds_DGE)

dds_statVals <- data.frame(
  EnsemblID = rownames(dds_DGE_results),
  stat = dds_DGE_results$stat
)
dds_statVals$ID <- rowData(dds)$gene_name[match(rownames(rowData(dds)), dds_statVals$EnsemblID)]
dds_statVals <- dds_statVals[!duplicated(dds_statVals$ID), ]
rownames(dds_statVals) <- NULL

deg <- dds_statVals %>%
  select(ID, stat) %>%
  column_to_rownames(var = 'ID') %>%
  as.matrix()

counts <- assay(dds)
counts <- counts[dds_statVals$EnsemblID, ]
rownames(counts) <- dds_statVals$ID

counts_logNorm <- log2(counts + 1)
colnames(counts_logNorm) <- sapply(colnames(counts_logNorm), function(x) substr(x,1,12))

design <- ann_tcga
names(design) <- c('sample','condition')

# Run fucking progeny
net <- get_progeny(organism = 'human', top = 100)

# Run mlm
contrast_acts <- run_mlm(mat = deg, net = net, .source = 'source',
                         .target = 'target', .mor = 'weight', minsize = 5)
contrast_acts

# Plot
g_progeny <- ggplot(contrast_acts, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab('Pathways') + ylab('Enrichment in HRD Samples')
ggsave(filename = '~/Projects/HRD_MutationalSignature/Figures/Figure2/HRD_decoupleR.pdf',
       plot = g_progeny, height = 4)
