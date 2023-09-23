#####
## Apply final HRD transcriptional signature to TCGA-BRCA test cohort
#####

# Load libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggpubr)
library(wesanderson)
library(tidyr)
library(pROC)

# Load data
load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata') # Group Annotation
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('Patient', 'HRD', 'BRCA_status')]
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('PALB2', 'RAD51C')), ]

ann_tcga$group <- ann_tcga$BRCA_status
ann_tcga$group[ann_tcga$HRD == 'HRD' & ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga$group[ann_tcga$group == 'none'] <- 'HR-proficient'

ann_tcga$HRD <- factor(ann_tcga$HRD, levels = c('HR-proficient','HRD'))
ann_tcga$group <- factor(ann_tcga$group, levels = c('HR-proficient', 'HRD_BRCA+',
                                                    'BRCA1', 'BRCA2'))

load('~/Data/TCGA/TCGA_BRCA.ExprDeconvolution_050823_p0.79_testing.Rdata') # reference dataset to define test cohort

# Define test annotations
samples.intersect <- intersect(substr(rownames(Z.tumor_testing),1,12), rownames(ann_tcga))
ann_tcga_test <- ann_tcga[samples.intersect, ]

# Obain non-deconvoluted test cohort via TCGAbiolinks

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

# Extract FPKM-normalised cohort, log2 normalise and match to annotation
expr.tumor_testing <- assay(expr.test, 'fpkm_uq_unstrand')
rownames(expr.tumor_testing) <- rowData(expr.test)$gene_name
colnames(expr.tumor_testing) <- sapply(colnames(expr.tumor_testing),
                                       function(x) substr(x,1,12))
expr.tumor_testing <- log2(expr.tumor_testing+1)

expr.tumor_testing <- t(expr.tumor_testing)
expr.tumor_testing <- expr.tumor_testing[rownames(ann_tcga_test), ]

# Load and extract the ElasticNet_alpha0.25 signature centroid
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
sig <- signature.centroid.list$ElasticNet_alpha0.25

# Calculate and plot HRD scores (including AUC values)
results_hrd <- data.frame(
  Patient = rownames(expr.tumor_testing),
  HRD_score = apply(expr.tumor_testing[,rownames(sig)], 1, function(x) cor(x,sig$HRD)) -
    apply(expr.tumor_testing[,rownames(sig)], 1, function(x) cor(x,sig$HR_proficient))
)

results_hrd <- merge(x = results_hrd, y = ann_tcga_test)

roc(HRD ~ HRD_score, data = results_hrd) # prints AUC estimate

g_hrdvshrd <- ggboxplot(data = results_hrd, x = 'HRD', y = 'HRD_score', fill = 'HRD') +
  stat_compare_means() + ylab('HRD score') +
  theme(legend.position = 'top', legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_HRDvsHRD.pdf',
       plot = g_hrdvshrd, width = 4, height = 4)

group_comparisons <- list(c('HR-proficient','BRCA2'),
                          c('HR-proficient','BRCA1'),
                          c('HR-proficient','HRD_BRCA+'))
g_hrdvsbrca <- ggboxplot(data = results_hrd, x = 'group', y = 'HRD_score', fill = 'group') +
  stat_compare_means(comparisons = group_comparisons) + ylab('HRD score') +
  theme(legend.position = 'top', legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_HRDvsBRCA.pdf',
       plot = g_hrdvsbrca, width = 5, height = 4)

# Calculate and plot BRCA defect-specific scores

results_brca <- data.frame(
  Patient = rownames(expr.tumor_testing),
  BRCA1 = apply(expr.tumor_testing[,rownames(sig)], 1, function(x) cor(x,sig$BRCA1)),
  BRCA2 = apply(expr.tumor_testing[,rownames(sig)], 1, function(x) cor(x,sig$BRCA2)),
  HRD_BRCApos = apply(expr.tumor_testing[,rownames(sig)], 1, function(x) cor(x,sig$HRD_BRCApos)),
  HR_proficient = apply(expr.tumor_testing[,rownames(sig)], 1, function(x) cor(x,sig$HR_BRCA_proficient))
)

results_brca <- merge(x = results_brca, y = ann_tcga_test)

results_brca.plot <- results_brca[,c(1:5,8)] %>%
  pivot_longer(cols = -c(Patient, group), names_to = 'Signature', values_to = 'Score')
results_brca.plot$Signature <- factor(results_brca.plot$Signature,
                                      levels = c('BRCA1','BRCA2',
                                                 'HRD_BRCApos','HR_proficient'))

g_brcavsbrca <- ggboxplot(data = results_brca.plot, x = 'group', y = 'Score', fill = 'group') +
  theme(legend.position = 'top', legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  facet_wrap(~Signature, scales = 'free')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_BRCAvsBRCA.pdf',
       plot = g_brcavsbrca, width = 8, height = 5)


## Analysis of reduced signature
genes.importance <- read.csv('~/Data/imp_score_avg.csv')
thres.imp <- 0.7
genes.important <- c(genes.importance$gene_names[genes.importance$HR.proficient > thres.imp],
                     genes.importance$gene_names[genes.importance$HRD > thres.imp])

sig_redux <- sig[genes.important,]

results_HRDredux <- data.frame(
  Patient = rownames(expr.tumor_testing),
  HRD_score = apply(expr.tumor_testing[,rownames(sig_redux)], 1, function(x) cor(x,sig_redux$HRD)) -
    apply(expr.tumor_testing[,rownames(sig_redux)], 1, function(x) cor(x,sig_redux$HR_proficient))
)

results_HRDredux <- merge(x = results_HRDredux, y = ann_tcga_test)
g_redux_hrdvshrd <- ggboxplot(data = results_HRDredux, x = 'HRD', y = 'HRD_score', fill = 'HRD') +
  stat_compare_means() + ylab('HRD score') +
  theme(legend.position = 'top', legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_TCGA_HRDvsHRD_redux.pdf',
       plot = g_redux_hrdvshrd, width = 4, height = 4)

group_comparisons <- list(c('HR-proficient','BRCA2'),
                          c('HR-proficient','BRCA1'),
                          c('HR-proficient','HRD_BRCA+'))
g_redux_hrdvsbrca <- ggboxplot(data = results_HRDredux, x = 'group', y = 'HRD_score', fill = 'group') +
  stat_compare_means(comparisons = group_comparisons) + ylab('HRD score') +
  theme(legend.position = 'top', legend.title = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_TCGA_HRDvsBRCA_redux.pdf',
       plot = g_redux_hrdvsbrca, width = 5, height = 4)

# Barplot of relevant importance values
genes.imp <- genes.importance[genes.importance$gene_names %in% genes.important, ]
genes.imp <- genes.imp[,-1]

genes.imp$Enriched <- sapply(genes.imp$HRD, function(x)
  ifelse(x > 0.7, 'HRD', 'HR-proficient'))
genes.imp$HR.proficient[genes.imp$Enriched == 'HRD'] <- 0
genes.imp$HRD[genes.imp$Enriched == 'HR-proficient'] <- 0
genes.imp$enrich_score <- genes.imp$HRD - genes.imp$HR.proficient
genes.imp <- genes.imp[order(genes.imp$enrich_score), ]
genes.imp$gene_names <- factor(genes.imp$gene_names,
                               levels = genes.imp$gene_names)

g_importance <- ggplot(genes.imp, aes(x = gene_names, y = enrich_score, fill = Enriched)) +
  geom_bar(stat = 'identity') +
  theme_minimal() + theme(axis.title.y = element_blank(),
                          axis.title.x = element_blank(),
                          legend.position = 'top',
                          legend.title = element_blank(),
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure4/Gene_ImportanceRank.pdf',
       plot = g_importance, height = 3)
