#####
## Application of Elastic Net signature to TCGA-BRCA testing data
#####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(ggpubr)
library(wesanderson)
library(tidyr)
library(pROC)

# Load in testing data and both HRD and BRCA-defect template
#   and subset data for template genes
load('Data/TCGA_BRCA/TCGA_BRCA_FPKM_filterByExpr_TPregress_testing.Rdata')

load('Results/TranscriptionalSignature/Templates/HRD/template_HRD_ElasticNet.Rdata')
load('Results/TranscriptionalSignature/Templates/BRCA_defect/template_BRCA_ElasticNet.Rdata')

input.x <- input_data.test[ ,rownames(template.brca)]
input.y_hrd <- input_data.test$HRD
input.y_brca <- factor(input_data.test$BRCA_defect,
                       levels = c('HR-proficient', 'HRD_BRCA+', 'BRCA1', 'BRCA2'))

# Apply HRD and BRCA-defect signatures to testing data
testing.HRD <- data.frame(
  HRD_status = input.y_hrd,
  BRCA_status = input.y_brca,
  HRD = apply(input.x, 1, function(x) cor(x, template.hrd$HRD)),
  HR_prof = apply(input.x, 1, function(x) cor(x, template.hrd$HR_proficient))
)
testing.HRD$HRD_score <- testing.HRD$HRD - testing.HRD$HR_prof

save(testing.HRD, file = 'Results/TranscriptionalSignature/SignatureScores_HRD/TCGA_testing_HRDscores.Rdata')

# Plotting
g_HRD_byHRD <- ggboxplot(data = testing.HRD, x = 'HRD_status', y = 'HRD_score',
                         fill = 'HRD_status') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray50') +
  theme_minimal() + theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          legend.position = 'top',
                          legend.title = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = 'Figures/Figure4/TCGA_testing_HRD_byHRD.pdf',
       plot = g_HRD_byHRD)

g_HRD_byBRCA <- ggboxplot(data = testing.HRD, x = 'BRCA_status', y = 'HRD_score',
                          fill = 'BRCA_status') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray50') +
  theme_minimal() + theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          legend.position = 'top',
                          legend.title = element_blank()) +
  scale_fill_manual(values = c(wes_palette('GrandBudapest1')[-2],
                               wes_palette('GrandBudapest2')[1]))
ggsave(filename = 'Figures/Figure4/TCGA_testing_HRD_byBRCA.pdf',
       plot = g_HRD_byBRCA, width = 5, height = 4)


# Repeat, with plotting for BRCA_status by BRCA_status
testing.BRCA <- data.frame(
  HRD_status = input.y_hrd,
  BRCA_status = input.y_brca,
  BRCA1 = apply(input.x, 1, function(x) cor(x, template.brca$BRCA1)),
  BRCA2 = apply(input.x, 1, function(x) cor(x, template.brca$BRCA2)),
  HRD_BRCApos = apply(input.x, 1, function(x) cor(x, template.brca$HRD_BRCApos)),
  HR_proficient = apply(input.x, 1, function(x) cor(x, template.brca$HR_proficient))
)

save(testing.BRCA, file = 'Results/TranscriptionalSignature/SignatureScores_BRCAdefect/TCGA_testing_BRCAscores.Rdata')


testing.BRCA.plot <- testing.BRCA %>%
  pivot_longer(cols = c('BRCA1','BRCA2','HRD_BRCApos','HR_proficient'),
               names_to = 'Signature', values_to = 'Correlation')
testing.BRCA.plot$Signature <- factor(testing.BRCA.plot$Signature,
                                      levels = c('BRCA1','BRCA2',
                                                 'HRD_BRCApos','HR_proficient'))

g_BRCA_byBRCA <- ggboxplot(data = testing.BRCA.plot, x = 'BRCA_status', y = 'Correlation',
                           fill = 'BRCA_status') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray50') +
  theme_minimal() + theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          legend.position = 'top',
                          legend.title = element_blank()) +
  scale_fill_manual(values = c(wes_palette('GrandBudapest1')[-2],
                               wes_palette('GrandBudapest2')[1])) +
  facet_wrap( ~ Signature)
ggsave(filename = 'Figures/Figure4/TCGA_testing_BRCA_byBRCA.pdf',
       width = 8, height = 4)

# Test HRD classification in breast cancer subtypes
#   especially ER+, ER-, TNBC

# Load clinical reference data and match to HRD testing results
load('Data/TCGA_BRCA_HRDHallmarks.Rdata')

testing.HRD$sample <- sapply(rownames(testing.HRD),
                             function(x) substr(x, 1, 12))
testing.HRD_BCsubtype <- merge(
  x = testing.HRD[,c('sample','HRD_status','HRD_score')],
  y = tcga_HRDHallmarks[,c('sample','er_status_by_ihc',
                           'her2_status_by_ihc','pr_status_by_ihc')])
testing.HRD_BCsubtype$TNBC <- testing.HRD_BCsubtype$er_status_by_ihc == 'Negative' &
  testing.HRD_BCsubtype$her2_status_by_ihc == 'Negative' &
  testing.HRD_BCsubtype$pr_status_by_ihc == 'Negative'

# Plot HRD-score vs HRD-status when subsetted on each data type
roc_hrd_erPos <- roc(HRD_status ~ HRD_score, data = testing.HRD_BCsubtype[testing.HRD_BCsubtype$er_status_by_ihc == 'Positive', ])
g_hrd_erPos <- ggboxplot(
  data = testing.HRD_BCsubtype[testing.HRD_BCsubtype$er_status_by_ihc == 'Positive', ],
  x = 'HRD_status', y = 'HRD_score', add = 'jitter', color = 'HRD_status'
) + stat_compare_means() +
  scale_color_manual(values = wes_palette('GrandBudapest1')) +
  theme(legend.title = element_blank(), axis.title.x = element_blank()) +
  ggtitle(paste0('ER-positive: AUC = ', round(roc_hrd_erPos$auc, 2)))

roc_hrd_erNeg <- roc(HRD_status ~ HRD_score, data = testing.HRD_BCsubtype[testing.HRD_BCsubtype$er_status_by_ihc == 'Negative', ])
g_hrd_erNeg <- ggboxplot(
  data = testing.HRD_BCsubtype[testing.HRD_BCsubtype$er_status_by_ihc == 'Negative', ],
  x = 'HRD_status', y = 'HRD_score', add = 'jitter', color = 'HRD_status'
) + stat_compare_means() +
  scale_color_manual(values = wes_palette('GrandBudapest1')) +
  theme(legend.title = element_blank(), axis.title.x = element_blank()) +
  ggtitle(paste0('ER-negative: AUC = ', round(roc_hrd_erNeg$auc, 2)))

roc_hrd_tnbc <- roc(HRD_status ~ HRD_score, data = testing.HRD_BCsubtype[testing.HRD_BCsubtype$TNBC, ])
g_hrd_tnbc <- ggboxplot(
  data = testing.HRD_BCsubtype[testing.HRD_BCsubtype$TNBC, ],
  x = 'HRD_status', y = 'HRD_score', add = 'jitter', color = 'HRD_status'
) + stat_compare_means() +
  scale_color_manual(values = wes_palette('GrandBudapest1')) +
  theme(legend.title = element_blank(), axis.title.x = element_blank()) +
  ggtitle(paste0('TNBC: AUC = ', round(roc_hrd_tnbc$auc, 2)))

g_hrd_BCsubtype <- ggarrange(plotlist = list(
  g_hrd_erPos, g_hrd_erNeg, g_hrd_tnbc), nrow = 1)
ggsave(filename = 'Figures/Supplementary/TCGA_testing_HRD_byHRD_BCsubtype.pdf',
       plot = g_hrd_BCsubtype, width = 9, height = 4)

## Test on reduced signature
genes_important.hrd <- c('SNRPA1','RNU7-171P','IFT22','BRCA1','RAB2A',
                         'ACVR1B','C18orf54','TRAF4','AP2A2','RAB6B',
                         'KCTD2','DZANK1','RAMAC')
genes_important.hrp <- c('IGKV2-30','CYB561D2','GSC','BRCA1','ERGIC3',
                         'MELTF','SNHG21','HNRNPCP1','TFDP1','SEPTIN3')
genes.important <- unique(c(genes_important.hrd, genes_important.hrp))

# Reduce template and testing datasets
template.hrd_redux <- template.hrd[genes.important, ]
input.x_redux <- input.x[,genes.important]

# Calculate HRD scores against reduced template
testing.HRD_redux <- data.frame(
  HRD_status = input.y_hrd,
  BRCA_status = input.y_brca,
  HRD = apply(input.x_redux, 1, function(x) cor(x, template.hrd_redux$HRD)),
  HR_prof = apply(input.x_redux, 1, function(x) cor(x, template.hrd_redux$HR_proficient))
)
testing.HRD_redux$HRD_score <- testing.HRD_redux$HRD - testing.HRD_redux$HR_prof

save(testing.HRD_redux, file = 'Results/TranscriptionalSignature/SignatureScores_HRD/TCGA_testing_reduxSig_HRDscores.Rdata')

# Plotting
g_HRD_byHRD_redux <- ggboxplot(data = testing.HRD_redux, x = 'HRD_status', y = 'HRD_score',
                         fill = 'HRD_status') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray50') +
  theme_minimal() + theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          legend.position = 'top',
                          legend.title = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  stat_compare_means()
ggsave(filename = 'Figures/Supplementary/TCGA_testing_HRD_byHRD_reduxSig.pdf',
       plot = g_HRD_byHRD_redux, width = 3, height = 3)

my_comps <- list(c('HR-proficient','HRD_BRCA+'),
                 c('HR-proficient','BRCA1'),
                 c('HR-proficient','BRCA2'))
g_HRD_byBRCA_redux <- ggboxplot(data = testing.HRD_redux, x = 'BRCA_status', y = 'HRD_score',
                          fill = 'BRCA_status') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray50') +
  theme_minimal() + theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          legend.position = 'top',
                          legend.title = element_blank()) +
  scale_fill_manual(values = c(wes_palette('GrandBudapest1')[-2],
                               wes_palette('GrandBudapest2')[1])) +
  stat_compare_means(comparisons = my_comps)
ggsave(filename = 'Figures/Supplementary/TCGA_testing_HRD_byBRCA_reduxSig.pdf',
       plot = g_HRD_byBRCA_redux, height = 3.5, width = 5)
