#####
## Scripts to compare centroid models against PARPi sensitivity in CCLE
#####

# Load libraries
library(data.table)
library(ggplot2)
library(tidyr)
library(ggpubr)

# Load CCLE expression data
setwd('~/Data/CCLE')

expr <- read.csv('Expression_Public_23Q2_subsetted.csv', row.names = 1)

# Load signature centroids (subset for ElasticNet_alpha0.25)
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
signature.centroid.list <- list(signature.centroid.list$ElasticNet_alpha0.25)
names(signature.centroid.list) <- c('ElasticNet_alpha0.25')

load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.Rdata')
names(signature_alternative.centroid.list) <- paste0('Alternative_',names(signature_alternative.centroid.list))

signature.centroid.list <- c(signature.centroid.list, signature_alternative.centroid.list)
rm(signature_alternative.centroid.list)

# Load drug sensitivity data
drugs <- read.csv('CCLE_breast_PRISM_drugSensitivity.csv', row.names = 1)
drugs <- drugs[,c(which(grepl(pattern = 'OLAPARIB', names(drugs))),
                  which(grepl(pattern = 'TALAZOPARIB', names(drugs))),
                  which(grepl(pattern = 'NIRAPARIB', names(drugs))),
                  which(grepl(pattern = 'RUCAPARIB', names(drugs))))]
names(drugs) <- sapply(names(drugs), function(x) strsplit(x,split='[..]')[[1]][1])
drugs$CellLine <- rownames(drugs)

# Initialise results matrix
res_olaparib <- data.frame(Signature = names(signature.centroid.list), Drug = 'olaparib', cor = NA, pVal = NA)
res_talazoparib <- data.frame(Signature = names(signature.centroid.list), Drug = 'talazoparib', cor = NA, pVal = NA)
res_niraparib <- data.frame(Signature = names(signature.centroid.list), Drug = 'niraparib', cor = NA, pVal = NA)
res_rucaparib <- data.frame(Signature = names(signature.centroid.list), Drug = 'rucaparib', cor = NA, pVal = NA)

# Correlate each signature with sensivitiy to each PARP inhibitor
for (i in 1:length(signature.centroid.list)) {
  
  print(names(signature.centroid.list)[i])
  
  sig.i <- signature.centroid.list[[i]] # Extract ith signature
  genes.intersect <- intersect(rownames(sig.i), colnames(expr))
  
  # if (i == 12) genes.intersect <- genes.intersect[genes.intersect != 'FAM170B']
  
  sig.i <- sig.i[genes.intersect, ]
  expr.i <- expr[, genes.intersect]
  
  # Calculate HRD scores
  expr.i.hrdScores <- data.frame(
    CellLine = rownames(expr.i),
    HRD = apply(expr.i, 1, function(x) cor(x,sig.i$HRD)),
    HR_proficient = apply(expr.i, 1, function(x) cor(x,sig.i$HR_proficient))
  )
  expr.i.hrdScores$HRD_score <- expr.i.hrdScores$HRD - expr.i.hrdScores$HR_proficient
  
  # Merge with drug data
  df.i <- merge(x = expr.i.hrdScores, y = drugs)
  
  # Save correlations in relevant results tables
  cor.ol <- cor.test(df.i$HRD_score, df.i$OLAPARIB)
  res_olaparib$cor[i] <- cor.ol$estimate
  res_olaparib$pVal[i] <- cor.ol$p.value
  
  cor.tal <- cor.test(df.i$HRD_score, df.i$TALAZOPARIB)
  res_talazoparib$cor[i] <- cor.tal$estimate
  res_talazoparib$pVal[i] <- cor.tal$p.value
  
  cor.nir <- cor.test(df.i$HRD_score, df.i$NIRAPARIB)
  res_niraparib$cor[i] <- cor.nir$estimate
  res_niraparib$pVal[i] <- cor.nir$p.value
  
  cor.ruc <- cor.test(df.i$HRD_score, df.i$RUCAPARIB)
  res_rucaparib$cor[i] <- cor.ruc$estimate
  res_rucaparib$pVal[i] <- cor.ruc$p.value
  
}

# Organise data table and plot relative significance values
library(data.table)
res <- rbindlist(list(res_olaparib,res_talazoparib,res_niraparib,res_rucaparib))
res$logP <- -log10(res$pVal)
res$group <- sapply(res$Signature, function(x) strsplit(x,split='_')[[1]][1])

library(ggplot2)
# g_pval <- ggplot(res, aes(x = Signature, y = logP, fill = group)) +
#   geom_bar(stat = 'identity') + geom_hline(yintercept = -log10(0.05), col = 'red', linetype = 'dashed') +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = 'top', legend.title = element_blank()) +
#   facet_wrap(~Drug, scales = 'free', nrow = 1)
# ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_CCLEsignificance.pdf',
#        plot = g_pval, height = 4, width = 7)

res$Method <- sapply(res$Signature, function(x) strsplit(x,split='_')[[1]][2])
res$Method[res$Method == 'alpha0.25'] <- 'ElasticNet'
res$Method <- factor(res$Method,
                     levels = c('ElasticNet','Severson','PARPi7','CIN70','Peng'))

g_pval <- ggplot(res, aes(x = Method, y = logP, fill = group)) +
  geom_bar(stat = 'identity') + geom_hline(yintercept = -log10(0.05), col = 'red', linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'top', legend.title = element_blank(), axis.title.x = element_blank()) +
  facet_wrap(~Drug, scales = 'free', nrow = 1) +
  scale_fill_brewer(palette = 'Pastel1')
ggsave(filename = '~/Projects/Thesis/Chapter 5/CCLE_methodComparisons.pdf',
       plot = g_pval, height = 4, width = 8)

# Plot correlate of ElasticNet_alpha0.25 score against PARPi response

sig.interest <- signature.centroid.list$ElasticNet_alpha0.25

genes.intersect <- intersect(rownames(sig.interest), colnames(expr))
sig.interest <- sig.interest[genes.intersect, ]
expr.hrd <- expr[,genes.intersect]

# Calculate HRD scores and match drug sensitivity data
df.hrd <- data.frame(
  CellLine = rownames(expr.hrd),
  HRD_score = apply(expr.hrd, 1, function(x) cor(x,sig.interest$HRD)) -
    apply(expr.hrd, 1, function(x) cor(x,sig.interest$HR_proficient))
)
df.hrd <- merge(x = df.hrd, y = drugs)

# Reorder df.hrd and plot correlations

df.hrd <- df.hrd %>%
  pivot_longer(cols = -c(CellLine, HRD_score),
               names_to = 'Drug', values_to = 'PRISM')
df.hrd$Drug <- factor(df.hrd$Drug,
                      levels = c('RUCAPARIB','NIRAPARIB',
                                 'TALAZOPARIB','OLAPARIB'))
g_ccle <- ggplot(df.hrd, aes(x = HRD_score, y = PRISM)) + 
  geom_point() + geom_smooth(method = 'lm', color = 'darkred') +
  stat_cor() + theme_minimal() +
  facet_wrap(~Drug, nrow = 1)
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure5/CCLE_PARPiResponse.pdf',
       plot = g_ccle, width = 8, height = 4)

