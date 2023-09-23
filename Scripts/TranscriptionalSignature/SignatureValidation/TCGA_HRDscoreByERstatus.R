#####
## Compare HRD transcriptional scores across ER-status
#####

# Load libraries
library(pROC)
library(ggplot2)
library(ggpubr)
library(wesanderson)

load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')
ann_tcga <- ann_tcga[ann_tcga$ER_status %in% c('Negative','Positive'),]

# Load reference testing data
load('~/Data/TCGA/TCGA_BRCA.ExprDeconvolution_050823_p0.79_testing.Rdata')
rownames(ann_tcga) <- ann_tcga$Patient

samples.intersect <- intersect(substr(rownames(Z.tumor_testing),1,12), rownames(ann_tcga))

ann_tcga_test <- ann_tcga[samples.intersect, ]

# Apply to non-deconvoluted samples
library(TCGAbiolinks)
library(SummarizedExperiment)

setwd('~/Data/TCGA')
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts',
  barcode = rownames(Z.tumor_testing)
)
# GDCdownload(query)
expr.test <- GDCprepare(query = query)
expr.test <- expr.test[,expr.test$sample_type == 'Primary Tumor']
expr.test <- expr.test[!duplicated(rowData(expr.test)$gene_name) &
                         !is.na(rowData(expr.test)$gene_name), ]

library(SummarizedExperiment)
expr.tumor_testing <- assay(expr.test, 'fpkm_uq_unstrand')
rownames(expr.tumor_testing) <- rowData(expr.test)$gene_name
colnames(expr.tumor_testing) <- sapply(colnames(expr.tumor_testing),
                                       function(x) substr(x,1,12))
expr.tumor_testing <- log2(expr.tumor_testing+1)

expr.tumor_testing <- t(expr.tumor_testing)
expr.tumor_testing <- expr.tumor_testing[rownames(ann_tcga_test), ]

# Load signature centroids
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
signature.of.interest <- signature.centroid.list$ElasticNet_alpha0.25

# Calculate HRD scores for testing data and match with ann_tcga_testing
hrdScore_func <- function(expr) {
  expr.hrd = expr[,rownames(signature.of.interest)]
  cor_hrd = apply(expr.hrd, 1, function(x) cor(x,signature.of.interest$HRD))
  cor_hrproficient = apply(expr.hrd, 1, function(x) cor(x,signature.of.interest$HR_proficient))
  hrdscores = cor_hrd - cor_hrproficient
  return(hrdscores)
}

res.df <- data.frame(
  Patient = rownames(expr.tumor_testing),
  HRD_score = hrdScore_func(expr.tumor_testing)
)

res.df <- merge(x = res.df, y = ann_tcga_test[,c('Patient','HRD','ER_status')])

# Calculate AUC values for positive and negative ER status
auc.erNeg <- roc(HRD ~ HRD_score, data = res.df[res.df$ER_status == 'Negative',])
auc.erPos <- roc(HRD ~ HRD_score, data = res.df[res.df$ER_status == 'Positive',])

# Plot results by ER status
g_erPos <- ggboxplot(data = res.df[res.df$ER_status == 'Positive',],
                     x = 'HRD', y = 'HRD_score', add = 'jitter', color = 'HRD')+
  scale_color_manual(values = wes_palette('GrandBudapest1')) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank()) +
  stat_compare_means() +
  ggtitle(paste0('ER-positive: AUC = ',round(auc.erPos$auc,2)))

g_erNeg <- ggboxplot(data = res.df[res.df$ER_status == 'Negative',],
                     x = 'HRD', y = 'HRD_score', add = 'jitter', color = 'HRD')+
  scale_color_manual(values = wes_palette('GrandBudapest1')) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank()) +
  stat_compare_means() +
  ggtitle(paste0('ER-negative: AUC = ',round(auc.erNeg$auc,2)))

g_join <- ggarrange(plotlist = list(g_erPos, g_erNeg))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/SupplementaryFigures/Supp_TCGA_HRDscoreByERstatus.pdf',
       plot = g_join, width = 7, height = 4)
