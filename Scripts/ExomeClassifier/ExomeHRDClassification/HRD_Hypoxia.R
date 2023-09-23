setwd('~/Data/TCGA/brca_tcga_pan_can_atlas_2018/')

library(introdataviz)
library(wesanderson)

load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')
hypoxia <- read.table('data_clinical_supp_hypoxia.txt', h=T)

ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status) &
                       !(ann_tcga$BRCA_status %in% c('PALB2','RAD51C')), ]
ann_tcga <- ann_tcga[ann_tcga$ER_status %in% c('Negative','Positive'), ]

df <- merge(x = ann_tcga[,c('Patient','HRD','BRCA_status','ER_status')],
            y = hypoxia[,c('PATIENT_ID','BUFFA_HYPOXIA_SCORE')],
            by.x = 'Patient', by.y = 'PATIENT_ID')
names(df)[5] <- 'Hypoxia_Buffa'
df$ER_status <- sapply(df$ER_status, function(x) 
  ifelse(x == 'Negative', 'ER-', 'ER+'))
df$ER_status <- factor(df$ER_status, levels = c('ER+', 'ER-'))

library(ggpubr)

df$Group <- df$BRCA_status
df$Group[df$BRCA_status == 'none' & df$HRD == 'HRD'] <- 'HRD_BRCA+'
df$Group[df$Group == 'none'] <- 'HR-proficient'
df$Group <- factor(df$Group, levels = c('HR-proficient','HRD_BRCA+',
                                        'BRCA1','BRCA2'))

comps <- list(c('HRD_BRCA+','BRCA2'),c('HR-proficient','HRD_BRCA+'),c('HRD_BRCA+','BRCA1'))

g_hypoxia <- ggplot(df, aes(x = Group, y = Hypoxia_Buffa, fill = ER_status)) +
  geom_split_violin(alpha = 0.4) + geom_boxplot() +
  stat_compare_means(comparisons = comps) +
  scale_fill_manual(values = wes_palette('Royal1')[c(3,2)]) +
  ylab('Hypoxia Score (Buffa)') + theme_minimal() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())
# ggsave(filename = '~/Projects/HRD_MutationalSignature/Figures/Figure2/HRD_Hypoxia.pdf',
#        plot = g_hypoxia, height = 4, width = 7)
ggsave(filename = '~/Projects/Thesis/Chapter 4/HRD_Hypoxia.pdf',
       plot = g_hypoxia)

ggplot(df[df$ER_status == 'ER+', ], aes(x = Group, y = Hypoxia_Buffa)) +
  geom_boxplot() +
  stat_compare_means(comparisons = comps)

lm_hypoxia <- lm(Hypoxia_Buffa ~ ER_status + BRCA_status*HRD, data = df)
anova(lm_hypoxia)

