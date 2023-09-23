#####
## Compare HRD exome classifiers
#####

setwd('~/Projects/HRD_MutationalSignature/Results')

# Load libraries
library(tidyr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(pROC)

# 1. My classifier
load('TCGA_HRD_resultsSummary.Rdata')
# results.tcga_df <- data.frame(
#   Tumor = results.tcga_df$Patient,
#   HRD79 = ifelse(results.tcga_df$HRD_prob >= 0.79,'HRD','HR-proficient'),
#   HRD50 = ifelse(results.tcga_df$HRD_prob >= 0.50,'HRD','HR-proficient'),
#   HRD32 = ifelse(results.tcga_df$HRD_prob >= 0.32,'HRD','HR-proficient')
# )
results.tcga_df <- data.frame(
  Tumor = results.tcga_df$Patient,
  p_HRDi = results.tcga_df$HRD_prob,
  HRDi = results.tcga_df$HRD
)

# 2. SigMA (run using webtool)
sigMA <- read.csv('SigMA_TCGA_output.csv')
sigMA <- data.frame(
  Tumor = sapply(sigMA$tumor, function(x) substr(x, 1, 12)),
  p_HRD_SigMA = sigMA$Signature_3_mva,
  SigMA = sigMA$pass_mva
)

# 3. deconstructSigs (run manually)
SBS3 <- read.table('TCGA_BRCA_deconstructSigs.txt', header = TRUE)
SBS3 <- data.frame(
  Tumor = sapply(rownames(SBS3), function(x) substr(x, 1, 12)),
  p_HRD_SBS3 = SBS3$SBS3,
  SBS3_dominant = apply(SBS3, 1, function(x) names(x)[x==max(x)]) == 'SBS3'
)

# 4. CX3 Copy Number Signatures (taken from Drews et al. Nature 2022)
mac.full <- read.delim('~/Data/TCGA/Macintyre_TCGAExposures.txt', sep = '\t')
mac.full <- mac.full[mac.full$Cancer == 'BRCA',c('Sample','Signature','Exposure')] %>% 
  pivot_wider(names_from = Signature, values_from = Exposure)
mac.full <- data.frame(
  Tumor = mac.full$Sample,
  p_HRD_CX3 = mac.full$CX3,
  CX3_dominant = apply(mac.full[,-1], 1, function(x) names(x)[x==max(x)]) == 'CX3'
)
names(mac.full) <- c('Tumor','p_HRD_CX3','CX3')

# 5. HRD index scores (taken from Marquard et al. 2015)
marq <- read.table('~/Data/TCGA/Marquard_HRDScores.txt', h=T)
marq$HRD_index <- apply(marq[,c('NtAI','LST','HRD.LOH')], 1, sum)
marq <- data.frame(
  Tumor = marq$Tumor,
  HRD_index = marq$HRD_index,
  HRDindex_42 = marq$HRD_index >= 42,
  HRDindex_63 = marq$HRD_index > 63
)

# Join all
df.full <- merge(x = results.tcga_df, y = sigMA, all.x = TRUE)
df.full <- merge(x = df.full, y = SBS3, all.x = TRUE)
df.full <- merge(x = df.full, y = mac.full, all.x = TRUE)
df.full <- merge(x = df.full, y = marq, all.x = TRUE)

# Add valieris BRCA defects
valieris <- read_excel('~/Data/TCGA/cancers-12-03687-s001/BRCA.xlsx', sheet = 'class-original')
valieris <- valieris[valieris$BRCA1_somatic_null != 'NA', ]
valieris <- valieris[,c('sample','event.BRCA1','event.BRCA2','event.RAD51C','event.PALB2')]
valieris$BRCA1 <- valieris$event.BRCA1 != 0
valieris$BRCA2 <- !(valieris$event.BRCA2 %in% c(0,'Mono-allelic-inactivation'))
valieris$RAD51C <- valieris$event.RAD51C != 0
valieris$PALB2 <- valieris$event.PALB2 != 0

valieris$BRCA_status <- 'none'
valieris$BRCA_status[valieris$PALB2] <- 'PALB2'
valieris$BRCA_status[valieris$RAD51C] <- 'RAD51C'
valieris$BRCA_status[valieris$BRCA2] <- 'BRCA2'
valieris$BRCA_status[valieris$BRCA1] <- 'BRCA1'

valieris <- data.frame(
  Tumor = valieris$sample,
  BRCA_defective = valieris$BRCA_status != 'none'
)

df.full <- merge(x = df.full, y = valieris, all.x = TRUE)

# Calculate sensitivities, specificities, and F-scores
break_table.func <- function(column) {
  t <- table(df.full[,column], df.full$BRCA_defective)
  
  true_negative = t[1,1]
  false_negative = t[1,2]
  false_positive = t[2,1]
  true_positive = t[2,2]
  
  return(c(true_positive, true_negative, false_positive, false_negative))
}

results <- data.frame(
  Method = c('HRDi','SigMA','SBS3_dominant','CX3','HRDindex_42','HRDindex_63'),
  TP = NA, TN = NA, FP = NA, FN = NA
)

# For each method: calculate total TP, TN, FP, FN 
for (i in 1:length(results$Method)) {
  method = results$Method[i]
  results[i,2:5] <- break_table.func(method)
}

# Calculate sensitivity, specificity, and balanced F-score of each method
weight = 1

results$recall = results$TP/(results$TP + results$FN)
results$precision = results$TP/(results$TP + results$FP)
results$F_score = (1+weight^2)*results$precision*results$recall/((weight^2)*results$precision + results$recall)
results$Method <- factor(results$Method, levels = results$Method[order(results$F_score, decreasing = TRUE)])

# Plot results
g_Fscore <- ggplot(results, aes(x = Method, y = F_score, fill = Method)) + 
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_fill_brewer(palette = 'Paired')
# ggsave(filename = '../Figures/Figure1/TCGA_HRDcomparisonsFscores.pdf', plot = g_Fscore,
#        width = 4.5, height = 4.5)

# Plot sensitivity and specificity
res.plot <- results[,c('Method','recall','precision')] %>%
  pivot_longer(-Method, names_to = 'Measure', values_to = 'value')
g_SensSpec <- ggplot(res.plot, aes(x = Method, y = value, fill = Method)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'top') +
  scale_fill_brewer(palette = 'Paired') +
  facet_wrap(~Measure, ncol = 1)

# ggsave(filename = '../Figures/Supp_TCGAcomparisonsSensSpec.pdf', plot = g_SensSpec,
#        width = 6, height = 3.6)
ggsave(filename = '~/Projects/Thesis/Chapter 3/TCGAcomparisons_SensSpec.pdf', plot = g_SensSpec,
       width = 4.5, height = 4.5)

# Plot comparative AUC curves for respective measures
roc_hrd <- plot(roc(df.full$BRCA_defective, df.full$p_HRDi), print.auc = TRUE)
roc_hrd <- plot(roc(df.full$BRCA_defective, df.full$p_HRD_SigMA), print.auc = TRUE,
                col = 'blue', print.auc.y = .45, add = TRUE)
roc_hrd <- plot(roc(df.full$BRCA_defective, df.full$p_HRD_SBS3), print.auc = TRUE,
                col = 'darkgreen', print.auc.y = .4, add = TRUE)
roc_hrd <- plot(roc(df.full$BRCA_defective, df.full$p_HRD_CX3), print.auc = TRUE,
                col = 'red', print.auc.y = .35, add = TRUE)
roc_hrd <- plot(roc(df.full$BRCA_defective, df.full$HRD_index), print.auc = TRUE,
                col = 'purple', print.auc.y = .3, add = TRUE)
