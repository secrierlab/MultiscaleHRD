#####
## Apply HRD exome classifier to TCGA-BRCA cohort
#####

setwd('~/Data/TCGA/')

# Load libraries
library(plyr)
library(TCGAbiolinks)
library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# Load TCGA data and tally mutation contributions
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Simple Nucleotide Variation',
  data.type = 'Masked Somatic Mutation'
)
# GDCdownload(query)
tcga_mutations <- GDCprepare(query = query)

# Exclude mutations from non-primary tumours
tcga_mutations$sample_type_code <- sapply(tcga_mutations$Tumor_Sample_Barcode,
                                          function(x) substr(x,14,15))
tcga_mutations <- tcga_mutations[tcga_mutations$sample_type_code == '01', ]

mut.maf <- read.maf(tcga_mutations, isTCGA = TRUE,
                    vc_nonSyn = names(table(tcga_mutations$Variant_Classification)))

mt_tally.tcga <- sig_tally(
  mut.maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg38',
  mode = 'ALL',
  use_syn = TRUE
)

# Check proportion of samples with SBS/ID loads >= 50 and plot distributions
sampleLoads <- data.frame(
  Sample = rownames(mt_tally.tcga$SBS_96),
  SBS = log2(apply(mt_tally.tcga$SBS_96, 1, sum)+1),
  ID = log2(apply(mt_tally.tcga$ID_83, 1, sum) + 1)
)

sampleLoads <- sampleLoads %>%
  pivot_longer(cols = -Sample, names_to = 'MutationType', values_to = 'log_count')
sampleLoads$MutationType <- factor(sampleLoads$MutationType,
                                   levels = c('SBS','ID'))

g_mutLoads <- gghistogram(data = sampleLoads, x = 'log_count',
                           fill = 'lightgrey') +
  geom_vline(xintercept = log2(50+1), col = 'red', linetype = 'dashed') +
  geom_vline(data = ddply(sampleLoads, "MutationType", summarize, wavg = median(log_count)), aes(xintercept=wavg),
             col = 'blue', linetype = 'dashed') +
  facet_wrap(~MutationType, scales = 'free')

ggsave(filename = '~/Projects/HRD_MutationalSignature/Figures/SupplementaryFigures/Supp_TCGAmutationLoads.pdf',
       plot = g_mutLoads, width = 6, height = 3)

table(sampleLoads$log_count[sampleLoads$MutationType == 'SBS'] >= log2(50+1))
table(sampleLoads$log_count[sampleLoads$MutationType == 'ID'] >= log2(50+1))

# Collate SBS_96 and ID_83 contributions and save results
tcga_muts <- as.data.frame(cbind(mt_tally.tcga$SBS_96, mt_tally.tcga$ID_83))

save(tcga_muts, file = 'TCGA_BRCA_mutContributions.Rdata')

load('TCGA_BRCA_mutContributions.Rdata') # if running without data loading

## Load relevant data for classifier

setwd('~/Projects/HRD_MutationalSignature/')

# Prior cluster mean distributions (cluster_distributions = mut.dists_mean)
load('Data/ClusterLikelihoods/ICGC_clust20_mclust_meanCont.Rdata')

# Signature Phenotype Assignment (cluster_assign = pheno_assigned)
load('Results/ICGC_BRCA_IDnormalised_PhenotypeAnnotation_clust20.Rdata')
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

results.tcga_df <- data.frame(
  Patient = rownames(results.tcga_loglik),
  Phenotype_Assigned = apply(results.tcga_loglik, 1,
                             function(x) names(x)[x==max(x)]),
  Phenotype_Assigned.prob = apply(results.tcga_loglik, 1, max),
  HRD_prob = apply(results.tcga_loglik[,grepl(pattern = 'HRD', names(results.tcga_loglik))],
                   1, sum)
)
results.tcga_df$HRD <- sapply(results.tcga_df$HRD_prob,
                              function(x) ifelse(x > .79, 'HRD', 'HR-proficient'))
# results.tcga_df$HRD_lenient <- sapply(results.tcga_df$HRD_prob,
#                                       function(x) ifelse(x >= .27, 'HRD', 'HR-proficient'))
# results.tcga_df$HRD_strict <- sapply(results.tcga_df$HRD_prob,
#                                      function(x) ifelse(x >= .79, 'HRD','HR-proficient'))
save(results.tcga_df, file = 'Results/TCGA_HRD_resultsSummary.Rdata')

# Compare with BRCA defects and clinical features

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

# Clinical subtypes
clin <- read.delim('~/Data/TCGA/TCGA_clinicalStatus.txt', header = TRUE)
tcga.clinical <- merge(x = valieris[,c('sample','BRCA_status')], y = clin[,c('bcr_patient_barcode','er_status_by_ihc')],
                       by.x = 'sample', by.y = 'bcr_patient_barcode')
names(tcga.clinical)[3] <- 'ER_status'

# Create annotation
ann_tcga <- merge(x = results.tcga_df[,c('Patient','Phenotype_Assigned', 'HRD', 'HRD_prob')],
                  y = tcga.clinical,
                  by.x = 'Patient', by.y = 'sample', all.x = TRUE)
ann_tcga <- ann_tcga[!duplicated(ann_tcga$Patient), ]
save(ann_tcga, file = 'Results/TCGA_HRDclassification_BRCAannotation.Rdata')

rownames(ann_tcga) <- ann_tcga$Patient; ann_tcga <- ann_tcga[,-1]
ann_tcga <- ann_tcga[order(ann_tcga$Phenotype_Assigned), ]

# Plots showing BRCA-defect classification as HRD
ann_tcga.plot <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga.plot$HRD_prob.plot <- ann_tcga.plot$HRD_prob+0.005
ann_tcga.plot <- ann_tcga.plot[order(ann_tcga.plot$HRD_prob), ]
ann_tcga.plot$index <- 1:nrow(ann_tcga.plot)

g_waterfall <- ggplot(ann_tcga.plot, aes(x = index, y = HRD_prob.plot, fill = BRCA_status)) +
  geom_bar(stat = 'identity') + theme_minimal() + ylab('p(HRD)') +
  scale_fill_manual(values = c('blue','red','lightgray','gold','darkgreen')) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'top') +
  geom_hline(yintercept = 0.79, col = 'red', linetype = 'dashed')
ggsave(filename = 'Figures/Figure1/TCGA_BRCAWaterfallPlot.pdf',
       plot = g_waterfall, height = 3, width = 7.5)

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
  BRCA_status = c('BRCA1' = 'blue', 'BRCA2' = 'red',
                  'PALB2' = 'gold', 'RAD51C' = 'darkgreen','none' = 'white'),
  ER_status = c('[Not Evaluated]' = 'white', 'Indeterminate' = 'navy',
                'Negative' = 'gold', 'Positive' = 'darkgreen')
)

# Use white -> navy scale
cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)

pheatmap(results.tcga_plot, 
         show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
         annotation_col = ann_tcga[,c(1,5,2:4)], 
         annotation_colors = ann_tcga_colours,
         color = cols_scale, fontsize = 8, fontsize_row = 10,
         filename = 'Figures/Figure1/TCGA_HRDclassificationHeatmapExtended.pdf')

table(ann_tcga$HRD, ann_tcga$BRCA_status != 'none')
table(ann_tcga$HRD, ann_tcga$BRCA_status %in% c('BRCA1','BRCA2'))
table(ann_tcga$HRD, ann_tcga$BRCA_status)

# Plot BRCA-defect/HRD status
ann_tcga$BRCA_status_broad <- sapply(ann_tcga$BRCA_status,
                                     function(x) ifelse(x == 'none', 'BRCA+', 'BRCA-defective'))

ann_tcga.plot <- ann_tcga[!is.na(ann_tcga$BRCA_status_broad), ] %>%
  group_by(HRD, BRCA_status_broad) %>% 
  summarise(n = n())
g_hrdBrca <- ggplot(ann_tcga.plot, aes(x = BRCA_status_broad, y = n, fill = HRD)) +
  geom_bar(stat = 'identity', position = 'fill') +
  theme(legend.position = 'top', axis.title.x = element_blank(),
        legend.title = element_blank()) +
  scale_fill_brewer(palette = 'Paired')
ggsave(filename = 'Figures/Figure1/TCGA_BRCASensitivity.pdf',
       plot = g_hrdBrca, width = 3.5, height = 3.5)

# BRCA type-specific HRD classification
hrd_brca1type <- c('HRD_APOBEC', 'HRD_ID6mid', 'HRD_ID8', 'HRD_SBS8')
hrd_brca2type <- c('HRD_ID6high')

ann_tcga$HRD_BRCAgroup <- ann_tcga$HRD
ann_tcga$HRD_BRCAgroup[ann_tcga$Phenotype_Assigned %in% hrd_brca1type] <- 'BRCA1-type HRD'
ann_tcga$HRD_BRCAgroup[ann_tcga$Phenotype_Assigned %in% hrd_brca2type] <- 'BRCA2-type HRD'
ann_tcga$HRD_BRCAgroup[ann_tcga$HRD_BRCAgroup == 'HRD'] <- 'HRD unassigned'

ann_tcga.plot2 <- ann_tcga[!is.na(ann_tcga$BRCA_status), ] %>%
  group_by(HRD_BRCAgroup, BRCA_status) %>%
  summarise(n=n())
ann_tcga.plot2$BRCA_status <- factor(ann_tcga.plot2$BRCA_status,
                                     levels = c('BRCA1','BRCA2',
                                                'RAD51C','PALB2','none'))
ann_tcga.plot2$HRD_BRCAgroup <- factor(ann_tcga.plot2$HRD_BRCAgroup,
                                       levels = c('HR-proficient','HRD unassigned',
                                                  'BRCA2-type HRD', 'BRCA1-type HRD'))

g_brcaHRDgroup <- ggplot(ann_tcga.plot2, aes(x = BRCA_status, y = n, fill = HRD_BRCAgroup)) +
  geom_bar(stat = 'identity', position = 'fill') +
  theme(legend.position = 'top', axis.title.x = element_blank(),
        legend.title = element_blank(), axis.title.y = element_blank()) +
  scale_fill_manual(values = c('gray90','gray50','red', 'blue'))
ggsave(filename = 'Figures/Supp_BRCAspecificHRD_BRCASensitivity.pdf',
       plot = g_brcaHRDgroup, width = 6, height = 4.5)
