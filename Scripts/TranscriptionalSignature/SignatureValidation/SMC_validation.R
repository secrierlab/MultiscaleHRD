library(readxl)
library(ggpubr)
library(wesanderson)
library(tidyr)

load('~/Projects/HRD_MutationalSignature/Results/SMC_HRD_resultsSummary.Rdata')
data.brca <- read_excel('~/Data/SMC_BRCA/SMC_BRCA.BrcaStatus.xlsx', skip=2)

results.smc_df$sample_id <- sapply(results.smc_df$Patient, 
                                   function(x) substr(x, 15, nchar(x)))
results.smc_df <- merge(x = results.smc_df, y = data.brca[,c('sample_id', 'gene_symbol')],
                        all.x = TRUE)
results.smc_df$group <- results.smc_df$gene_symbol
results.smc_df$group[results.smc_df$HRD_prob > 0.79 & 
                       is.na(results.smc_df$gene_symbol)] <- 'HRD_BRCA+'
results.smc_df$group[is.na(results.smc_df$group)] <- 'HR-proficient'
results.smc_df$group <- factor(results.smc_df$group,
                               levels = c('HR-proficient','HRD_BRCA+',
                                          'BRCA1','BRCA2'))

# Load validation data and match with HRD classifications
smc_rnaseq <- read.delim('~/Data/SMC_BRCA/data_mrna_seq_tpm.txt', sep='\t')
smc_rnaseq <- smc_rnaseq[!duplicated(smc_rnaseq$Hugo_Symbol), ]
rownames(smc_rnaseq) <- smc_rnaseq$Hugo_Symbol
smc_rnaseq <- smc_rnaseq[,3:ncol(smc_rnaseq)]

samples.intersect <- intersect(colnames(smc_rnaseq), results.smc_df$Patient)
smc_rnaseq <- smc_rnaseq[,samples.intersect]
results.smc_df <- results.smc_df[match(samples.intersect, results.smc_df$Patient), ]

# Load signature centroids
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
sig.interest <- signature.centroid.list$ElasticNet_alpha0.25

# Calculate HRD and BRCA-defect scores
genes.intersect <- intersect(rownames(sig.interest), rownames(smc_rnaseq))
sig.interest <- sig.interest[genes.intersect, ]
smc_rnaseq_hrd <- smc_rnaseq[genes.intersect, ]

df.hrd <- data.frame(
  Patient = colnames(smc_rnaseq_hrd),
  HRD_score = apply(smc_rnaseq_hrd, 2, function(x) cor(x,sig.interest$HRD)) -
    apply(smc_rnaseq_hrd, 2, function(x) cor(x,sig.interest$HR_proficient)),
  BRCA1 = apply(smc_rnaseq_hrd, 2, function(x) cor(x,sig.interest$BRCA1)),
  BRCA2 = apply(smc_rnaseq_hrd, 2, function(x) cor(x,sig.interest$BRCA2)),
  HRD_BRCApos = apply(smc_rnaseq_hrd, 2, function(x) cor(x,sig.interest$HRD_BRCApos)),
  HR_BRCA_proficient = apply(smc_rnaseq_hrd, 2, function(x) cor(x,sig.interest$HR_BRCA_proficient))
)

# Match with HRD classification
df.hrd <- merge(x = df.hrd, y = results.smc_df[,c(2,6,8)])
g_HRDbyHRD <- ggboxplot(df.hrd, x = 'HRD', y = 'HRD_score', fill = 'HRD') +
  stat_compare_means() + scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_HRDvsHRD.pdf',
       plot = g_HRDbyHRD, width = 4, height = 4)

brca_comparisons <- list(c('HR-proficient','HRD_BRCA+'),
                         c('HR-proficient','BRCA1'),
                         c('HR-proficient','BRCA2'))
g_HRDvsBRCA <- ggboxplot(df.hrd, x = 'group', y = 'HRD_score', fill = 'group') +
  stat_compare_means(comparisons = brca_comparisons) + 
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_HRDvsBRCA.pdf',
       plot = g_HRDvsBRCA, width = 6, height = 4)

# Match with BRCA classifications
df.hrd_plot <- df.hrd[,c(1,3:6,8)] %>%
  pivot_longer(cols = -c(Patient,group), names_to = 'Signature', values_to = 'Correlation')
g_BRCAvsBRCA <- ggboxplot(df.hrd_plot, x = 'group', y = 'Correlation', fill = 'group') +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank()) +
  facet_wrap(~Signature, scales = 'free')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_BRCAvsBRCA.pdf',
       plot = g_BRCAvsBRCA, width = 8, height = 5)
