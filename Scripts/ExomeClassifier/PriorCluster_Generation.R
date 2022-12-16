#####
## Generate prior distributions for each signature phenotype cluster
#####

setwd('~/Documents/GitHub/HRD_classification')

# Load libraries
library(tidyr)
library(deconstructSigs) # for context ordering
library(ggplot2)

# 1. Load mutation tallies
load('Data/ICGC_BRCA/BRCA_UKEU_mt_tally.Rdata')
mut_complete <- as.data.frame(cbind(mt_tally.brca_wgs$SBS_96, mt_tally.brca_wgs$ID_83))

# 2. Load annotation (ann)
load('Data/ICGC_BRCA/ICGC_BRCA_PhenotypeAnnotation.Rdata')
ann <- ann[rownames(mut_complete), ]

# Separate mutation counts into clusters
mut_byClust <- list()
for (pheno in unique(ann$Phenotype)) {
  mut_byClust[[pheno]] <- mut_complete[ann$Phenotype == pheno, ]
}

# For each input, generate a complete probability distribution, and collate them
#   NB, whilst summing seemed nice, SigMA used a mean. Let's try both!

collate_function <- function(input, variant_type = 'ALL', collation = 'total') {
  
  sbs.index <- sapply(names(input), function(x) grepl('>',x))
  id.index <- sapply(names(input), function(x) grepl(':',x))
  
  if (variant_type == 'SBS') {
    input_final <- input[,sbs.index]
  } else if (variant_type == 'ID') {
    input_final <- input[,id.index]
  } else {
    input_final <- input
  }
  
  # Add pseudocount
  input_final <- input_final + 1 # original
  
  if (collation == 'mean') {
    dist_temp = t(apply(input_final, 1, function(x) x/sum(x)))
    dist = apply(dist_temp, 2, mean)
  } else {
    if (collation != 'total') print('Set collation to mean or total. Setting to total...')
    dist_temp = apply(input_final, 2, sum)
    dist = dist_temp/sum(dist_temp)
  }
  
  return(dist)
}

mut.dists_total = mut.dists_mean <- data.frame()

for (pheno in names(mut_byClust)) {
  mut.dists_total <- rbind(mut.dists_total, collate_function(input = mut_byClust[[pheno]]))
  mut.dists_mean <- rbind(mut.dists_mean, collate_function(mut_byClust[[pheno]], collation = 'mean'))
}

colnames(mut.dists_total) = colnames(mut.dists_mean) <- colnames(mut_byClust[[1]])
rownames(mut.dists_total) = rownames(mut.dists_mean) <- names(mut_byClust)

# Save prior clusters
save(mut.dists_mean, file = 'Data/PriorClusters/ICGC_Clust19_mclust_meanCont.Rdata')
save(mut.dists_total, file = 'Data/PriorClusters/ICGC_Clust19_mclust_totalCont.Rdata')

# Plot prior clusters
mut.dists_mean.plot <- cbind(mut.dists_mean[,colnames(signatures.cosmic)],
                             mut.dists_mean[,97:ncol(mut.dists_mean)])
mut.dists_mean.plot$Phenotype <- rownames(mut.dists_mean)
mut.dists_mean.plot <- mut.dists_mean.plot %>%
  pivot_longer(cols = -Phenotype, names_to = 'Context', values_to = 'Contribution')
mut.dists_mean.plot$Phenotype <- factor(mut.dists_mean.plot$Phenotype,
                                        levels = sort(rownames(mut.dists_mean)))
mut.dists_mean.plot$Context <- factor(mut.dists_mean.plot$Context,
                                      levels = c(colnames(signatures.cosmic), colnames(mut.dists_mean)[97:ncol(mut.dists_mean)]))
mut.dists_mean.plot$Type <- ifelse(
  grepl(pattern = '>', mut.dists_mean.plot$Context),
  'SBS','indel')

g_meanPlot <- ggplot(mut.dists_mean.plot, aes(x = Context, y = Contribution, fill = Type)) +
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(axis.text.x = element_blank(), legend.position = 'top') +
  facet_wrap(~Phenotype, scales = 'free')
ggsave(filename = 'Figures/Supplementary/PriorDistributions_mean.pdf',
       plot = g_meanPlot, height = 5, width = 10)


# mut.dists_total.plot <- mut.dists_total
# mut.dists_total.plot[,1:96] <- mut.dists_total.plot[,colnames(signatures.cosmic)]
# mut.dists_total.plot$Phenotype <- rownames(mut.dists_total)
# mut.dists_total.plot <- mut.dists_total.plot %>%
#   pivot_longer(cols = -Phenotype, names_to = 'Context', values_to = 'Contribution')
# mut.dists_total.plot$Phenotype <- factor(mut.dists_total.plot$Phenotype,
#                                         levels = sort(rownames(mut.dists_total)))
# mut.dists_total.plot$Context <- factor(mut.dists_total.plot$Context,
#                                       levels = c(colnames(signatures.cosmic), colnames(mut.dists_total)[97:ncol(mut.dists_total)]))
# mut.dists_total.plot$Type <- ifelse(
#   grepl(pattern = '>', mut.dists_total.plot$Context),
#   'SBS','indel')
# 
# g_totalPlot <- ggplot(mut.dists_total.plot, aes(x = Context, y = Contribution, fill = Type)) +
#   geom_bar(stat = 'identity') + theme_minimal() +
#   theme(axis.text.x = element_blank(), legend.position = 'top') +
#   facet_wrap(~Phenotype, scales = 'free')
# ggsave(filename = 'Figures/Supplementary/PriorDistributions_total.pdf',
#        plot = g_totalPlot, height = 5, width = 10)
