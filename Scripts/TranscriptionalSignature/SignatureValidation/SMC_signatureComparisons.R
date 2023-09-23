library(pROC)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(readxl)

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
smc_rnaseq <- as.data.frame(t(smc_rnaseq))
results.smc_df <- results.smc_df[match(samples.intersect, results.smc_df$Patient), ]

load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.Rdata')
names(signature_alternative.centroid.list) <- paste0('Alternative_',names(signature_alternative.centroid.list))

signature.centroid.list <- c(list(signature.centroid.list[[1]]), signature_alternative.centroid.list)
rm(signature_alternative.centroid.list)
names(signature.centroid.list)[1] <- 'ElasticNet'

auc.df <- data.frame(
  Model = names(signature.centroid.list),
  BRCA1 = NA, BRCA2 = NA, HRD_BRCApos = NA, HR_BRCA_proficient = NA, 
  HRD = NA
)

for (i in 1:length(signature.centroid.list)) {
  
  print(names(signature.centroid.list)[i])
  
  sig.i <- signature.centroid.list[[i]]
  
  genes.include <- intersect(rownames(sig.i), colnames(smc_rnaseq))
  sig.i <- sig.i[genes.include, ]
  
  # Start with BRCA
  sig.group <- strsplit(names(signature.centroid.list)[i],split = '_')[[1]][1]
  
  results_brca.testing <- data.frame(
    group = results.smc_df$group,
    BRCA1 = apply(smc_rnaseq[,genes.include], 1, function(x) cor(x,sig.i$BRCA1)),
    BRCA2 = apply(smc_rnaseq[,genes.include], 1, function(x) cor(x,sig.i$BRCA2)),
    HRD_BRCApos = apply(smc_rnaseq[,genes.include], 1, function(x) cor(x,sig.i$HRD_BRCApos)),
    HR_proficient = apply(smc_rnaseq[,genes.include], 1, function(x) cor(x,sig.i$HR_BRCA_proficient))
  )
  auc.df$BRCA1[i] <- roc(group == 'BRCA1' ~ BRCA1, data = results_brca.testing)$auc
  auc.df$BRCA2[i] <- roc(group == 'BRCA2' ~ BRCA2, data = results_brca.testing)$auc
  auc.df$HRD_BRCApos[i] <- roc(group == 'HRD_BRCA+' ~ HRD_BRCApos, data = results_brca.testing)$auc
  auc.df$HR_BRCA_proficient[i] <- roc(group == 'HR-proficient' ~ HR_proficient, data = results_brca.testing)$auc
  
  # Then HRD
  results_hrd.testing <- data.frame(
    HRD_status = results.smc_df$HRD,
    HRD = apply(smc_rnaseq[,genes.include], 1, function(x) cor(x,sig.i$HRD)),
    HR_proficient = apply(smc_rnaseq[,genes.include], 1, function(x) cor(x,sig.i$HR_proficient))
  )
  results_hrd.testing$HRD_score <- results_hrd.testing$HRD - results_hrd.testing$HR_proficient
  auc.df$HRD[i] <- roc(HRD_status == 'HRD' ~ HRD_score, data = results_hrd.testing)$auc
  
}

gene.markers <- c('BRCA1','BRCA2','POLQ','PARP1')
auc_gene.df <- data.frame(
  Model = paste0('Gene_',gene.markers),
  BRCA1 = NA, BRCA2 = NA, HRD_BRCApos = NA, HR_BRCA_proficient = NA, 
  HRD = NA
)

for (i in 1:length(gene.markers)) {
  
  df.i <- results.smc_df
  df.i$expr_i <- smc_rnaseq[,gene.markers[i]]
  
  auc_gene.df$BRCA1[i] <- roc(group  == 'BRCA1' ~ expr_i, data = df.i)$auc
  auc_gene.df$BRCA2[i] <- roc(group  == 'BRCA2' ~ expr_i, data = df.i)$auc
  auc_gene.df$HRD_BRCApos[i] <- roc(group  == 'HRD_BRCA+' ~ expr_i, data = df.i)$auc
  auc_gene.df$HR_BRCA_proficient[i] <- roc(group  == 'HR-proficient' ~ expr_i, data = df.i)$auc
  
  auc_gene.df$HRD[i] <- roc(HRD == 'HRD' ~ expr_i, data = df.i)$auc
  
}

auc.df <- rbind(auc.df, auc_gene.df)

auc.df <- auc.df[order(auc.df$HRD, decreasing = TRUE), ]
auc.df$Model <- factor(auc.df$Model, levels = auc.df$Model)
# auc.df$group <- sapply(auc.df$Model, function(x) strsplit(x,split='_')[[1]][1])
ggplot(auc.df, aes(x = Model, y = HRD, fill = Model)) + geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'Spectral') +
  coord_cartesian(ylim = c(0.45,1.05)) +
  theme(axis.text.x = element_blank())

auc.df_plot <- auc.df[,-6] %>%
  pivot_longer(cols = -c(Model), names_to = 'Signature', values_to = 'AUC')
ggplot(auc.df_plot, aes(x = Model, y = AUC, fill = Model)) + geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'Spectral') +
  theme(axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(0.45,1.05)) +
  facet_wrap(~Signature, scales = 'free', nrow = 1)

# Plot for figures
auc.df <- auc.df[order(auc.df$HRD, decreasing = TRUE), ]
auc.df$Model <- factor(auc.df$Model, levels = auc.df$Model)

g1 <- ggplot(auc.df, aes(x = Model, y = HRD, fill = Model)) +
  geom_bar(stat = 'identity') + ylim(c(0,1)) + ylab('AUC') +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))

auc.dfb <- auc.df[,1:5] %>%
  pivot_longer(cols = -Model, names_to = 'Signature', values_to = 'AUC')
auc.dfb$Signature <- factor(auc.dfb$Signature,
                            levels = c('BRCA1','BRCA2',
                                       'HRD_BRCApos', 'HR_BRCA_proficient'))
g2 <- ggplot(auc.dfb, aes(x = Model, y = AUC, fill = Model)) +
  geom_bar(stat = 'identity') + ylim(c(0,1)) + 
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18)) +
  facet_wrap(~Signature, nrow = 1)

g_join <- ggarrange(plotlist = list(g1,g2))
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_SMC_AUCresults.pdf',
       plot = g_join, width = 16, height = 6)
