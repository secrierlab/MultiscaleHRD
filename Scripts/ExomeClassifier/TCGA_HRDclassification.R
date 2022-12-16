#####
## Apply HRD exome classifier to TCGA-BRCA cohort
#####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Load TCGA data and tally mutation contributions
data <- read.maf(maf = '~/Documents/PhD/Data/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz',
                 isTCGA = TRUE)

mut.maf <- rbind(data@data, data@maf.silent)
mut.maf <- read.maf(mut.maf, 
                    vc_nonSyn = names(table(mut.maf$Variant_Classification)))

mt_tally.tcga <- sig_tally(
  mut.maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg38',
  mode = 'ALL',
  use_syn = TRUE
)

# Check proportion of samples with SBS/ID loads >= 50 and plot distributions
sampleLoads <- data.frame(
  Sample = rownames(mt_tally.tcga$SBS_96),
  SBS = apply(mt_tally.tcga$SBS_96, 1, sum),
  ID = apply(mt_tally.tcga$ID_83, 1, sum),
  log_SBS = log2(apply(mt_tally.tcga$SBS_96, 1, sum)+1),
  log_ID = log2(apply(mt_tally.tcga$ID_83, 1, sum) + 1)
)
g_sbs <- gghistogram(data = sampleLoads, x = 'log_SBS',
                     fill = 'lightgrey') + xlab('log2(#SBS + 1)') +
  geom_vline(xintercept = log2(50+1), col = 'red', linetype = 'dashed')
g_id <- gghistogram(data = sampleLoads, x = 'log_ID',
                    fill = 'lightgrey') + xlab('log2(#ID + 1)')  +
  geom_vline(xintercept = log2(50+1), col = 'red', linetype = 'dashed')
g_sampleLoads <- ggarrange(g_sbs, g_id)
ggsave(filename = 'Figures/Supplementary/TCGA_mutationLoads.pdf',
       plot = g_sampleLoads, width = 6)

table(sampleLoads$SBS >= 50)
table(sampleLoads$ID >= 50)

# Collate SBS_96 and ID_83 contributions
tcga_muts <- as.data.frame(cbind(mt_tally.tcga$SBS_96, mt_tally.tcga$ID_83))

## Load relevant data for classifier
# Prior cluster mean distributions (cluster_distributions = mut.dists_mean)
load('Data/PriorClusters/ICGC_Clust19_mclust_meanCont.Rdata')

# Signature Phenotype Assignment (cluster_assign = pheno_assigned)
load('Data/ICGC_BRCA/ICGC_BRCA_PhenotypeAnnotation.Rdata')
pheno_assigned <- ann$Phenotype

# Likelihood function that aligns a dataset with the designated mean distributions
likelihood_calc <- function(input_data, cluster_distributions, cluster_assign) {
  
  # This function:
  #   Takes a dataset of mutations as its input (rows = samples, cols = mutation types)
  #   For now, this function does not limit mutation types: SBS and indels will be included
  #   Applies a likelihood approach to calculate posterior probabilities of cluster assignment
  #   NB, For now, we ensure that our cluster distributions have the correct mutation types
  #   Output: posterior distributions
  
  # Calculate likelihood: P(mutation spectrum | cluster)
  log_likelihoods <- matrix(NA, nrow = nrow(input_data), ncol = nrow(cluster_distributions))
  rownames(log_likelihoods) <- rownames(input_data); colnames(log_likelihoods) <- rownames(cluster_distributions)
  print('Calculating log-likelihoods...')
  for (i in 1:nrow(input_data)) {
    print(paste0('Calculating log likelihoods for sample ', i, ' of ', nrow(input_data), ': ', rownames(input_data)[i]))
    log_likelihoods[i, ] <- apply(cluster_distributions, 1, 
                                  function(x) sum(log10(x) * input_data[i, ]))
  }
  
  # Set prior: P(cluster)
  marginal.probs <- table(cluster_assign)/length(cluster_assign)
  marginal.probs <- marginal.probs[colnames(log_likelihoods)]
  
  # Calculate log_posteriors and return final posteriors
  log_posteriors <- data.frame(t(apply(log_likelihoods, 1, function(x) 
    x + log10(marginal.probs)
  )))
  
  # Generate final posteriors
  final_probs <- log_posteriors
  for (i in 1:ncol(final_probs)) {
    final_probs[,i] <- 10^(apply(log_posteriors, 1, function(x) x[i] - max(x[1:ncol(final_probs)])))
  }
  
  final_probs <- data.frame(t(apply(final_probs, 1, function(x) x/sum(x))))
  
  return(final_probs)
  
}

# Apply log-likelihood approach
results.tcga_loglik <- likelihood_calc(input_data = tcga_muts, 
                                       cluster_distributions = mut.dists_mean,
                                       cluster_assign = pheno_assigned)
save(results.tcga_loglik,
     file = 'Results/ExomeClassifier/TCGA_BRCA/TCGA_HRD_logLikelihoods.Rdata')

results.tcga_df <- data.frame(
  Patient = rownames(results.tcga_loglik),
  Phenotype_Assigned = apply(results.tcga_loglik, 1,
                             function(x) names(x)[x==max(x)]),
  Phenotype_Assigned.prob = apply(results.tcga_loglik, 1, max),
  HRD_prob = apply(results.tcga_loglik[,grepl(pattern = 'HRD', names(results.tcga_loglik))],
                   1, sum)
)
results.tcga_df$HRD <- sapply(results.tcga_df$HRD_prob,
                              function(x) ifelse(x > .5, 'HRD', 'HR-proficient'))
save(results.tcga_df, file = 'Results/ExomeClassifier/TCGA_BRCA/TCGA_HRD_resultsSummary.Rdata')

## Compare with BRCA defects and clinical features

# Load BRCA-defect data from Valieris et al.
valieris <- read_excel('~/Downloads/cancers-12-03687-s001/BRCA.xlsx',
                       sheet = 'class-original')
valieris <- valieris[valieris$BRCA1_somatic_null != 'NA', ]
valieris <- valieris[,c('sample','event.BRCA1','event.BRCA2')]
valieris$BRCA1 <- valieris$event.BRCA1 != 0
valieris$BRCA2 <- !(valieris$event.BRCA2 %in% c(0,'Mono-allelic-inactivation'))

valieris$BRCA_status <- 'none'
valieris$BRCA_status[valieris$BRCA1] <- 'BRCA1'
valieris$BRCA_status[valieris$BRCA2] <- 'BRCA2'

# Clinical subtypes
load('~/Documents/PhD/Data/TCGA_BRCA_ClinicalSubtypes.Rdata')
dat.patients_brcaSubtype <- dat.patients_brcaSubtype[,1:2]

tcga.clinical <- merge(x = valieris[,c(c('sample','BRCA_status'))], y = dat.patients_brcaSubtype,
                       by.x = 'sample', by.y = 'bcr_patient_barcode')

# Create annotation
ann_tcga <- merge(x = results.tcga_df[,c('Patient','Phenotype_Assigned', 'HRD', 'HRD_prob')],
                  y = tcga.clinical,
                  by.x = 'Patient', by.y = 'sample', all.x = TRUE)
ann_tcga <- ann_tcga[!duplicated(ann_tcga$Patient), ]

rownames(ann_tcga) <- ann_tcga$Patient; ann_tcga <- ann_tcga[,-1]
ann_tcga <- ann_tcga[order(ann_tcga$Phenotype_Assigned), ]

# Plots showing BRCA-defect classification as HRD
ann_tcga.plot <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga.plot$HRD_prob.plot <- ann_tcga.plot$HRD_prob - 0.5
ann_tcga.plot <- ann_tcga.plot[order(ann_tcga.plot$HRD_prob), ]
ann_tcga.plot$index <- 1:nrow(ann_tcga.plot)

g_waterfall <- ggplot(ann_tcga.plot,
                      aes(x = index, y = HRD_prob.plot, fill = BRCA_status)) +
  geom_bar(stat = 'identity') + theme_minimal() + ylab('p(HRD)') +
  scale_fill_manual(values = c('blue','red','lightgray')) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'top')
ggsave(filename = 'Figures/Figure2/TCGA_HRDclassification_BRCAwaterfall.pdf', g_waterfall,
       height = 3, width = 6)

# Barplot of BRCA-defect HRD classification
ann_tcga.plot$BRCA_defective <- ifelse(ann_tcga.plot$BRCA_status != 'none',
                                       'BRCA_defective', 'BRCA+')
ann_tcga.plot_hrd <- ann_tcga.plot %>%
  group_by(HRD, BRCA_defective) %>% summarise(n = n())
g_hrd <- ggplot(ann_tcga.plot_hrd, aes(x = BRCA_defective, y = n, fill = HRD)) +
  geom_bar(stat = 'identity', position = 'fill') +
  theme_minimal() + scale_fill_brewer(palette = 'Paired') +
  theme(axis.title.x = element_blank(),
        legend.position = 'top',
        legend.title = element_blank()) +
  ylab('% Sample')
ggsave(filename = 'Figures/Figure2/TCGA_HRDclassification_BRCAdefective.pdf', plot = g_hrd,
       height = 3, width = 3)

# Barplot of BRCA-defect-specific HRD classification
hrd_brca1type <- c('HRD_ID8', 'HRD_APOBEC', 'HRD_IDmult')
hrd_brca2type <- c('HRD_ID6mid', 'HRD_ID6high')

ann_tcga.plot$HRD_BRCAspecific <- 'HR-proficient'
ann_tcga.plot$HRD_BRCAspecific[ann_tcga.plot$Phenotype_Assigned %in% hrd_brca1type] <- 'BRCA1-type HRD'
ann_tcga.plot$HRD_BRCAspecific[ann_tcga.plot$Phenotype_Assigned %in% hrd_brca2type] <- 'BRCA2-type HRD'

ann_tcga.plot_brca <- ann_tcga.plot %>%
  group_by(HRD_BRCAspecific, BRCA_status) %>% summarise(n = n())
ann_tcga.plot_brca$HRD_BRCAspecific <- factor(ann_tcga.plot_brca$HRD_BRCAspecific,
                                              levels = c('HR-proficient',
                                                         'BRCA2-type HRD',
                                                         'BRCA1-type HRD'))
g_brca <- ggplot(ann_tcga.plot_brca, aes(x = BRCA_status, y = n, fill = HRD_BRCAspecific)) +
  geom_bar(stat = 'identity', position = 'fill') +
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(5,'mm')) +
  scale_fill_manual(values = c('lightgrey','red','blue')) +
  ylab('% Samples')
ggsave(filename = 'Figures/Supplementary/TCGA_HRDclassification_BRCAspecific.pdf', plot = g_brca,
       height = 3, width = 3)

# Continue with heatmap and save annotation
names(ann_tcga)[5] <- 'ER_status'
ann_tcga$Phenotype_Assigned <- as.factor(ann_tcga$Phenotype_Assigned)
ann_tcga$HRD <- as.factor(ann_tcga$HRD)
ann_tcga$BRCA_status <- as.factor(ann_tcga$BRCA_status)
ann_tcga$ER_status <- as.factor(ann_tcga$ER_status)
save(ann_tcga, file = 'Results/ExomeClassifier/TCGA_BRCA/TCGA_HRD_annotation.Rdata')

## Plot heatmap of results

# Reorder results
results.tcga_plot <- as.data.frame(t(results.tcga_loglik[rownames(ann_tcga), order(colnames(results.tcga_loglik))]))

# Sort colours
cols <- colorRampPalette(brewer.pal(9,'Set1'))
cols_pheno <- cols(length(unique(ann_tcga$Phenotype_Assigned)))
names(cols_pheno) <- unique(ann_tcga$Phenotype_Assigned)

ann_tcga_colours <- list(
  Phenotype_Assigned = cols_pheno,
  HRD = c('HRD' = 'black', 'HR-proficient' = 'white'),
  BRCA_status = c('BRCA1' = 'blue', 'BRCA2' = 'red', 'none' = 'white'),
  ER_status = c('Positive' = 'darkgreen', 'Negative' = 'yellow',
                'Indeterminate' = 'navy', '[Not Evaluated]' = 'white')
)

# Use white -> navy scale
cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)

pheatmap(results.tcga_plot, 
         show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
         annotation_col = ann_tcga[,c(1,5,2:4)], 
         annotation_colors = ann_tcga_colours,
         color = cols_scale, fontsize = 8, fontsize_row = 10,
         filename = 'Figures/Figure2/TCGA_HRDclassification_heatmap.pdf')
