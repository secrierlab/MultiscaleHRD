#####
## Plotting hallmarks of HRD
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(wesanderson)

# Load TCGA HRD classifications and hallmarks, and combine
load('Results/TCGA_HRDclassification_BRCAannotation.Rdata')
tcga_HRDHallmarks <- read.table('~/Data/TCGA/HRD_hallmarks.txt', h=T, sep='\t')

resultsHRD <- merge(x = ann_tcga[,c('Patient', 'Phenotype_Assigned', 'HRD')], 
                    y = tcga_HRDHallmarks, all.x = TRUE)
resultsHRD$HRD <- factor(resultsHRD$HRD,
                         levels = c('HR-proficient', 'HRD'))

# # Plot 4 main HRD hallmarks
# resultsHRD.plot <- resultsHRD[,c('Patient','HRD','HRD_index',
#                                  'CX3','log_POLQ_FPKM','ProliferativeCapacity')] %>%
#   pivot_longer(cols = -c(Patient,HRD), names_to = 'Hallmark')
# resultsHRD.plot$Hallmark <- factor(resultsHRD.plot$Hallmark,
#                                    levels = c('HRD_index', 'CX3',
#                                               'log_POLQ_FPKM','ProliferativeCapacity'))
# g_hallmarks <- ggboxplot(data = resultsHRD.plot, x = 'HRD', y = 'value', fill = 'HRD') + 
#   stat_compare_means() + 
#   scale_fill_manual(values = wes_palette('GrandBudapest1')) +
#   theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
#         legend.position = 'top', legend.title = element_blank(),
#         axis.text.x = element_blank()) +
#   facet_wrap(~Hallmark, scales = 'free', nrow = 1)
# ggsave(filename = 'Figures/Figure2/HRDHallmarks.pdf', height = 4)

# MYC amplification
resultsHRD$MYC_status <- factor(resultsHRD$MYC_status, 
                                levels = c('Deletion', 'Normal','Amplification'))
res_MYC <- resultsHRD %>%
  group_by(HRD, MYC_status) %>% summarise(n = n()) %>% drop_na()
res_MYC

g_myc <- ggplot(data = res_MYC, aes(x = HRD, y = n, fill = MYC_status)) +
  geom_bar(stat = 'identity', position = 'fill') +
  theme_minimal() + ylab('% Samples') +
  theme(legend.position = 'top', axis.title.y = element_blank()) + coord_flip() +
  scale_fill_manual(values = c(wes_palette('Zissou1')[3], 'gray90', wes_palette('Zissou1')[1]))
ggsave('Figures/Figure2/HRD_vs_MYCamplification.pdf', g_myc,
       width = 6, height = 3)

# POLQ expression
g_polq <- ggboxplot(data = resultsHRD, x = 'HRD', y = 'log_POLQ_FPKM', fill = 'HRD') + 
  stat_compare_means() + 
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank()) 
ggsave('Figures/Figure2/HRD_vs_POLQexpression.pdf', g_polq,
       width = 4, height = 4)

# HRD index scores
g_hrdIndex <- ggboxplot(data = resultsHRD, x = 'HRD', y = 'HRD_index', fill = 'HRD') + 
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank())
ggsave('Figures/Figure2/HRD_vs_HRDscore.pdf', g_hrdIndex,
       width = 4, height = 4)

# Individual HRD index scores
df_hrdIndex <- resultsHRD[,c('HRD', 'NtAI', 'LST', 'HRD.LOH')] %>%
  pivot_longer(cols = -HRD, names_to = 'HRD_index', values_to = 'score')
g_hrdIndex_individual <- ggboxplot(data = df_hrdIndex, x = 'HRD', y = 'score', fill = 'HRD') +
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  facet_wrap(~HRD_index)
ggsave('Figures/Supp_HRDvsHRDscore_ind.pdf', g_hrdIndex_individual,
       width = 8, height = 4)

# CX3 Copy Number Signature Exposure
g_cx3 <- ggboxplot(data = resultsHRD, x = 'HRD', y = 'CX3', fill = 'HRD') + 
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank())
ggsave('Figures/Figure2/HRD_vs_CX3Exposure.pdf', g_cx3,
       width = 4, height = 4)

# # CX2/5 Copy Number Signature Exposure
# resultsHRD_CXs <- resultsHRD[,c('Patient','HRD','CX2','CX5')] %>%
#   pivot_longer(cols = c(CX2,CX5), names_to = 'CX_Signature', values_to = 'Exposure')
# g_cxs <- ggboxplot(data = resultsHRD_CXs, x = 'HRD', y = 'Exposure', fill = 'HRD') +
#   stat_compare_means() +
#   scale_fill_manual(values = wes_palette('GrandBudapest1')) +
#   theme(axis.title.x = element_blank()) +
#   facet_wrap(~CX_Signature)
# ggsave('Figures/Supplementary/HRD_vs_CXExposure.pdf', g_cxs,
#        width = 5, height = 3.5)

# Quiescence Score Comparison
g_prol <- ggboxplot(data = resultsHRD, x = 'HRD', y = 'ProliferativeCapacity', fill = 'HRD') +
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank())
ggsave('Figures/Figure2/HRD_vs_Quiescence.pdf', g_prol,
       width = 4, height = 4)


# Breast Cancer Subtypes
res.BRCA <- resultsHRD[,c('Patient', 'HRD', 'er_status_by_ihc',
                          'her2_status_by_ihc', 'pr_status_by_ihc')]

res.BRCA$BRCA_subtype <- NA
res.BRCA$BRCA_subtype[res.BRCA$er_status_by_ihc == 'Positive'] <- 'ER+'
res.BRCA$BRCA_subtype[res.BRCA$her2_status_by_ihc == 'Positive'] <- 'HER2+'
res.BRCA$BRCA_subtype[res.BRCA$er_status_by_ihc == 'Positive' &
                        res.BRCA$her2_status_by_ihc == 'Positive'] <- 'ER+HER2+'
res.BRCA$BRCA_subtype[res.BRCA$er_status_by_ihc == 'Negative' &
                        res.BRCA$her2_status_by_ihc == 'Negative' &
                        res.BRCA$pr_status_by_ihc == 'Negative'] <- 'TNBC'

res.BRCA.plot <- res.BRCA %>%
  drop_na(BRCA_subtype) %>%
  group_by(HRD, BRCA_subtype) %>% summarise(n=n())

g_resBRCA <- ggplot(res.BRCA.plot, aes(x = BRCA_subtype, y = n, fill = HRD)) +
  geom_bar(stat = 'identity', position = 'fill') + theme_minimal() +
  xlab('') + ylab('% Samples') + theme(legend.position = 'top') +
  scale_fill_manual(values = wes_palette('GrandBudapest1'))
ggsave(filename = 'Figures/Figure2/HRDvsBRCAsubtype.pdf',
       plot = g_resBRCA, width = 4, height = 4)

g_resBRCA2 <- ggplot(res.BRCA.plot, aes(x = HRD, y = n, fill = BRCA_subtype)) +
  geom_bar(stat = 'identity', position = 'fill') + theme_minimal() +
  xlab('') + ylab('% Samples') + theme(legend.position = 'top') + 
  coord_flip() +
  scale_fill_manual(values = wes_palette('Royal2')[4:1])
ggsave(filename = 'Figures/Figure2/BRCAsubtypevsHRD.pdf',
       plot = g_resBRCA2, width = 6, height = 3)

## Plot supplementary figures

# 1. HRD hallmarks: HR-proficient vs HRD_BRCA+ vs HRD_BRCA-
ann_tcga2 <- ann_tcga
ann_tcga2 <- ann_tcga2[!is.na(ann_tcga$BRCA_status), ]
ann_tcga2$group <- ann_tcga2$HRD
ann_tcga2$group[ann_tcga2$HRD == 'HRD' & 
                  ann_tcga2$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga2$group[ann_tcga2$HRD == 'HRD' & 
                  ann_tcga2$BRCA_status != 'none'] <- 'HRD_BRCA-'
ann_tcga2$group <- factor(ann_tcga2$group,
                          levels = c('HR-proficient','HRD_BRCA+','HRD_BRCA-'))

df_supp1 <- merge(x = ann_tcga2[,c('Patient','group')],
                  y = tcga_HRDHallmarks[,c('Patient','HRD_index','CX3',
                                           'log_POLQ_FPKM','ProliferativeCapacity')])
df_supp1 <- df_supp1 %>%
  pivot_longer(cols = -c(Patient,group), names_to = 'Hallmark', values_to = 'score')

comps_supp1 <- list(c('HR-proficient','HRD_BRCA+'),
                    c('HRD_BRCA+','HRD_BRCA-'),
                    c('HR-proficient','HRD_BRCA-'))

g_brcaPos <- ggboxplot(df_supp1, x = 'group', y = 'score', fill = 'group') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'top') +
  stat_compare_means(comparisons = comps_supp1) +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  facet_wrap(~Hallmark, scales = 'free', nrow=1)

ggsave(filename = 'Figures/Supp_TCGAhallmarksBRCAPositive.pdf',
       plot = g_brcaPos, width = 10, height = 5)


# 1. HRD hallmarks: HR-proficient vs p(HRD) > 0.5 vs p(HRD) > 0.79
ann_tcga3 <- ann_tcga[,c(1,3,4)]
ann_tcga3$group <- sapply(ann_tcga3$HRD_prob, function(x)
  ifelse(x > 0.79, 'p(HRD) > 0.79',
         ifelse(x > 0.5, 'p(HRD) > 0.5', 'HR-proficient')))
ann_tcga3$group <- factor(ann_tcga3$group,
                          levels = c('HR-proficient',
                                     'p(HRD) > 0.5', 'p(HRD) > 0.79'))

df_supp2 <- merge(x = ann_tcga3[,c('Patient','group')],
                  y = tcga_HRDHallmarks[,c('Patient','HRD_index','CX3',
                                           'log_POLQ_FPKM','ProliferativeCapacity')])
df_supp2 <- df_supp2 %>%
  pivot_longer(cols = -c(Patient,group), names_to = 'Hallmark', values_to = 'score')

comps_supp2 <- list(c('HR-proficient','p(HRD) > 0.5'),
                    c('p(HRD) > 0.5','p(HRD) > 0.79'),
                    c('HR-proficient','p(HRD) > 0.79'))

g_hrdThresholds <- ggboxplot(df_supp2, x = 'group', y = 'score', fill = 'group') +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'top') +
  stat_compare_means(comparisons = comps_supp2) +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  facet_wrap(~Hallmark, scales = 'free', nrow=1)

ggsave(filename = 'Figures/Supp_TCGAhallmarksHRDthresholds.pdf',
       plot = g_hrdThresholds, width = 10, height = 5)
