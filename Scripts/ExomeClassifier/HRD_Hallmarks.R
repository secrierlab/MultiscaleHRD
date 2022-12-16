#####
## Plotting hallmarks of HRD
#####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(wesanderson)

# Load TCGA HRD classifications and hallmarks, and combine
load('Results/ExomeClassifier/TCGA_BRCA/TCGA_HRD_resultsSummary.Rdata')
load('Data/TCGA_BRCA_HRDHallmarks.Rdata')

resultsHRD <- merge(x = results.tcga_df[,c('Patient', 'Phenotype_Assigned', 'HRD')], 
                    y = tcga_HRDHallmarks,
                    by.x = 'Patient', by.y = 'sample', all.x = TRUE)
resultsHRD$HRD <- factor(resultsHRD$HRD,
                         levels = c('HR-proficient', 'HRD'))

# MYC amplification
res_MYC <- resultsHRD %>%
  group_by(HRD, MYC) %>% summarise(n = n()) %>% drop_na()
res_MYC

g_myc <- ggplot(data = res_MYC, aes(x = HRD, y = n, fill = MYC)) +
  geom_bar(stat = 'identity', position = 'fill') +
  theme_minimal() + ylab('% Samples') +
  theme(axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = c(wes_palette('Zissou1')[4], 'gray90', wes_palette('Zissou1')[1]))
ggsave('Figures/Figure3/HRD_vs_MYCamplification.pdf', g_myc,
       width = 4, height = 4)

# POLQ expression
g_polq <- ggboxplot(data = resultsHRD, x = 'HRD', y = 'POLQ', fill = 'HRD') + 
  stat_compare_means() + 
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank()) 
ggsave('Figures/Figure3/HRD_vs_POLQexpression.pdf', g_polq,
       width = 4, height = 4)

# HRD index scores
g_hrdIndex <- ggboxplot(data = resultsHRD, x = 'HRD', y = 'HRD_index', fill = 'HRD') + 
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank())
ggsave('Figures/Figure3/HRD_vs_HRDscore.pdf', g_hrdIndex,
       width = 4, height = 4)

# Individual HRD index scores
df_hrdIndex <- resultsHRD[,c('HRD', 'NtAI', 'LST', 'HRD.LOH')] %>%
  pivot_longer(cols = -HRD, names_to = 'HRD_index', values_to = 'score')
g_hrdIndex_individual <- ggboxplot(data = df_hrdIndex, x = 'HRD', y = 'score', fill = 'HRD') +
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank()) +
  facet_wrap(~HRD_index)
ggsave('Figures/Supplementary/HRD_vs_HRDscore_ind.pdf', g_hrdIndex_individual,
       width = 8, height = 4)

# CX3 Copy Number Signature Exposure
g_cx3 <- ggboxplot(data = resultsHRD, x = 'HRD', y = 'CX3', fill = 'HRD') + 
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank())
ggsave('Figures/Figure3/HRD_vs_CX3Exposure.pdf', g_cx3,
       width = 4, height = 4)

# CX2/5 Copy Number Signature Exposure
resultsHRD_CXs <- resultsHRD[,c('Patient','HRD','CX2','CX5')] %>%
  pivot_longer(cols = c(CX2,CX5), names_to = 'CX_Signature', values_to = 'Exposure')
g_cxs <- ggboxplot(data = resultsHRD_CXs, x = 'HRD', y = 'Exposure', fill = 'HRD') +
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank()) +
  facet_wrap(~CX_Signature)
ggsave('Figures/Supplementary/HRD_vs_CXExposure.pdf', g_cxs,
       width = 5, height = 3.5)

# Quiescence Score Comparison
g_qui <- ggboxplot(data = resultsHRD, x = 'HRD', y = 'Quiescence_Score', fill = 'HRD') +
  stat_compare_means() +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank())
ggsave('Figures/Figure3/HRD_vs_Quiescence.pdf', g_qui,
       width = 4, height = 4)

# APOBEC enrichment scores
resultsHRD$Phenotype_HRD_APOBEC <- sapply(resultsHRD$Phenotype_Assigned, function(x)
  strsplit(x, split = '_')[[1]][1])
resultsHRD$Phenotype_HRD_APOBEC[resultsHRD$Phenotype_Assigned == 'HRD_APOBEC'] <- 'HRD_APOBEC'
resultsHRD$Phenotype_HRD_APOBEC <- factor(resultsHRD$Phenotype_HRD_APOBEC,
                                          levels = c('APOBEC','HRD_APOBEC','HRD',
                                                     'SBS5','ID2','ID4'))
g_apobec <- ggboxplot(data = resultsHRD, x = 'Phenotype_HRD_APOBEC', y = 'APOBEC_enrich', fill = 'Phenotype_HRD_APOBEC') +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())
ggsave('Figures/Supplementary/HRD_vs_APOBEC.pdf', g_apobec, 
       width = 6.5, height = 3)

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
ggsave(filename = 'Figures/Figure3/HRD_vs_BRCAsubtype.pdf',
       plot = g_resBRCA, width = 4, height = 4)
