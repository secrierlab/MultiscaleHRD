#####
## Compare transcriptional signatures against TCGA-BRCA training cohort
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

# Load complete set of signatures
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.Rdata')
names(signature_alternative.centroid.list) <- paste0('Alternative_',names(signature_alternative.centroid.list))

signature.centroid.list <- c(signature.centroid.list, signature_alternative.centroid.list)
rm(signature_alternative.centroid.list)

# Initialise data frame to collate AUC values
auc.df <- data.frame(
  Model = names(signature.centroid.list),
  BRCA1 = NA, BRCA2 = NA, HRD_BRCApos = NA, HR_BRCA_proficient = NA, 
  HRD = NA
)

# For each model:
#   Calculate HRD and BRCA-specific signature scores
#   Calculate the relevant AUC for each statistic

for (i in 1:length(signature.centroid.list)) {
  
  print(names(signature.centroid.list)[i])
  
  sig.i <- signature.centroid.list[[i]]
  genes.include <- rownames(sig.i)
  
  # Start with BRCA
  sig.group <- strsplit(names(signature.centroid.list)[i],split = '_')[[1]][1]
  
  results_brca.testing <- data.frame(
    group = ann_tcga_test$group,
    BRCA1 = apply(expr.tumor_testing[,genes.include], 1, function(x) cor(x,sig.i$BRCA1)),
    BRCA2 = apply(expr.tumor_testing[,genes.include], 1, function(x) cor(x,sig.i$BRCA2)),
    HRD_BRCApos = apply(expr.tumor_testing[,genes.include], 1, function(x) cor(x,sig.i$HRD_BRCApos)),
    HR_proficient = apply(expr.tumor_testing[,genes.include], 1, function(x) cor(x,sig.i$HR_BRCA_proficient))
  )
  
  auc.df$BRCA1[i] <- roc(group == 'BRCA1' ~ BRCA1, data = results_brca.testing)$auc
  auc.df$BRCA2[i] <- roc(group == 'BRCA2' ~ BRCA2, data = results_brca.testing)$auc
  auc.df$HRD_BRCApos[i] <- roc(group == 'HRD_BRCA+' ~ HRD_BRCApos, data = results_brca.testing)$auc
  auc.df$HR_BRCA_proficient[i] <- roc(group == 'HR-proficient' ~ HR_proficient, data = results_brca.testing)$auc
  
  # Then HRD
  results_hrd.testing <- data.frame(
    HRD_status = ann_tcga_test$HRD,
    HRD = apply(expr.tumor_testing[,genes.include], 1, function(x) cor(x,sig.i$HRD)),
    HR_proficient = apply(expr.tumor_testing[,genes.include], 1, function(x) cor(x,sig.i$HR_proficient))
  )
  results_hrd.testing$HRD_score <- results_hrd.testing$HRD - results_hrd.testing$HR_proficient
  auc.df$HRD[i] <- roc(HRD_status == 'HRD' ~ HRD_score, data = results_hrd.testing)$auc
  
}

# Repeat for relevant gene markers and add to AUC results

gene.markers <- c('BRCA1','BRCA2','POLQ','PARP1')
auc_gene.df <- data.frame(
  Model = paste0('Gene_',gene.markers),
  BRCA1 = NA, BRCA2 = NA, HRD_BRCApos = NA, HR_BRCA_proficient = NA, 
  HRD = NA
)

for (i in 1:length(gene.markers)) {
  
  df.i <- ann_tcga_test
  df.i$expr_i <- expr.tumor_testing[,gene.markers[i]]
  
  auc_gene.df$BRCA1[i] <- roc(group  == 'BRCA1' ~ expr_i, data = df.i)$auc
  auc_gene.df$BRCA2[i] <- roc(group  == 'BRCA2' ~ expr_i, data = df.i)$auc
  auc_gene.df$HRD_BRCApos[i] <- roc(group  == 'HRD_BRCA+' ~ expr_i, data = df.i)$auc
  auc_gene.df$HR_BRCA_proficient[i] <- roc(group  == 'HR-proficient' ~ expr_i, data = df.i)$auc
  
  auc_gene.df$HRD[i] <- roc(HRD == 'HRD' ~ expr_i, data = df.i)$auc
  
}

auc.df <- rbind(auc.df, auc_gene.df)

# For plotting, remove other regression signatures, and plot AUC in descending order of HRD prediction
auc.df_plot <- auc.df[c(1,5:nrow(auc.df)), ]
auc.df_plot$Signature <- sapply(auc.df_plot$Model, function(x) strsplit(x,split='_')[[1]][2])
auc.df_plot$Signature[1] <- 'ElasticNet'
auc.df_plot <- auc.df_plot[order(auc.df_plot$HRD, decreasing = TRUE), ]
auc.df_plot$Signature <- factor(auc.df_plot$Signature, levels = auc.df_plot$Signature)

# Plot HRD AUCs
g_auc_hrd <- ggplot(auc.df_plot, aes(x = Signature, y = HRD, fill = Signature)) + 
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(0.5, 1)) + ylab('AUC') +
  scale_fill_brewer(palette = 'Spectral')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure3/TCGA_AUCbyHRD.pdf',
       plot = g_auc_hrd, height = 4, width = 7)

# Plot BRCA-specific AUCs
auc.df_plot2 <- auc.df_plot[,c(7,2:5)] %>%
  pivot_longer(cols = -Signature, names_to = 'group', values_to = 'AUC')
auc.df_plot2$group <- factor(auc.df_plot2$group,
                             levels = c('BRCA1','BRCA2',
                                        'HRD_BRCApos','HR_BRCA_proficient'))
g_auc_brca <- ggplot(auc.df_plot2, aes(x = Signature, y = AUC, fill = Signature)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_fill_brewer(palette = 'Spectral') +
  facet_wrap(~ group, nrow = 1)
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_TCGA_AUCbyBRCA.pdf',
       plot = g_auc_brca, width = 7, height = 4)
