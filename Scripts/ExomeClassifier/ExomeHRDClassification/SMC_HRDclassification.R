#####
## Apply HRD exome classifier to SMC BRCA cohort
#####

setwd('~/Data/SMC_BRCA/')

# Load libraries
library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(readxl)

# Load SMC data and tally mutation contributions
smc_mutations <- read.delim('data_mutations.txt')

mut.maf <- read.maf(smc_mutations, 
                    vc_nonSyn = names(table(smc_mutations$Variant_Classification)))

mt_tally.smc <- sig_tally(
  mut.maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ALL',
  use_syn = TRUE
)

smc_muts <- as.data.frame(cbind(mt_tally.smc$SBS_96, mt_tally.smc$ID_83))

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
results.smc_loglik <- likelihood_calc(input_data = smc_muts, 
                                       cluster_distributions = mut.dists_mean,
                                       cluster_assign = pheno_assigned)

results.smc_df <- data.frame(
  Patient = rownames(results.smc_loglik),
  Phenotype_Assigned = apply(results.smc_loglik, 1,
                             function(x) names(x)[x==max(x)]),
  Phenotype_Assigned.prob = apply(results.smc_loglik, 1, max),
  HRD_prob = apply(results.smc_loglik[,grepl(pattern = 'HRD', names(results.smc_loglik))],
                   1, sum)
)
results.smc_df$HRD <- sapply(results.smc_df$HRD_prob,
                              function(x) ifelse(x >= 0.79, 'HRD', 'HR-proficient'))
save(results.smc_df, file = 'Results/SMC_HRD_resultsSummary.Rdata')

# Match with BRCA status
data.brca <- read_excel('~/Data/SMC_BRCA/SMC_BRCA.BrcaStatus.xlsx', skip=2)

results.smc_df$sample_id <- sapply(results.smc_df$Patient, 
                                   function(x) substr(x, 15, nchar(x)))
results.smc_df <- merge(x = results.smc_df, y = data.brca[,c('sample_id', 'gene_symbol')],
                        all.x = TRUE)

# Match with additional clinical data, specifically BRCA subtype
smc_clinic <- read.delim('~/Data/SMC_BRCA/data_clinical_sample.txt', skip=4)

results.smc_df <- merge(x = results.smc_df, y = smc_clinic[,c('PATIENT_ID', 'SUBTYPE_CONSENSUS')],
                        by.x = 'Patient', by.y = 'PATIENT_ID')
names(results.smc_df)[ncol(results.smc_df)] <- 'Subtype'
names(results.smc_df)[ncol(results.smc_df)-1] <- 'BRCA_defect'

# Plot results
ann_smc <- results.smc_df[,c('BRCA_defect','HRD_prob', 'HRD', 'Phenotype_Assigned', 'Subtype')]
rownames(ann_smc) <- results.smc_df$Patient
ann_smc <- ann_smc[order(ann_smc$Phenotype_Assigned), ]

results.smc_plot <- as.data.frame(t(results.smc_loglik[rownames(ann_smc), order(colnames(results.smc_loglik))]))

# Sort colours
cols <- colorRampPalette(brewer.pal(9,'Set1'))
cols_pheno <- cols(length(unique(ann_smc$Phenotype_Assigned)))
names(cols_pheno) <- unique(ann_smc$Phenotype_Assigned)

ann_smc_colours <- list(
  Phenotype_Assigned = cols_pheno,
  HRD = c('HRD' = 'black', 'HR-proficient' = 'white'),
  BRCA_defect = c('BRCA1' = 'blue', 'BRCA2' = 'red'),
  Subtype = c('ER+' = 'navy', 'HER2+' = 'darkgreen',
              'ER+HER2+' = 'gray', 'TN' = 'yellow')
)
ann_smc <- ann_smc[,ncol(ann_smc):1]

# Use white -> navy scale
cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)

pheatmap(results.smc_plot, 
         show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE,
         annotation_col = ann_smc, 
         annotation_colors = ann_smc_colours,
         color = cols_scale, fontsize = 8, fontsize_row = 10,
         filename = 'Figures/Supp_SMCheatmap.pdf')

# Compare with BRCA subtype
table(ann_smc$HRD_prob >= 0.79, ann_smc$BRCA_defect, useNA = 'always')

# Plotting barplots
ann_smc$BRCA_status <- sapply(ann_smc$BRCA_defect, function(x)
  ifelse(!is.na(x), 'BRCA-defective', 'BRCA+'))
ann_smc$BRCA_status <- factor(ann_smc$BRCA_status,
                              levels = c('BRCA-defective','BRCA+'))

ann_smc$HRDgroup <- 'HR-proficient'
ann_smc$HRDgroup[ann_smc$HRD_prob >= 0.5] <- 'HRD > 0.5'
ann_smc$HRDgroup[ann_smc$HRD_prob >= 0.79] <- 'HRD > 0.79'
ann_smc$HRDgroup <- factor(ann_smc$HRDgroup,
                           levels = c('HR-proficient','HRD > 0.5', 'HRD > 0.79'))

ann_smc.plot1 <- ann_smc %>%
  group_by(BRCA_status, HRDgroup) %>%
  summarise(n = n())

g_brca <- ggplot(ann_smc.plot1, aes(x = BRCA_status, y = n, fill = HRDgroup)) +
  geom_bar(stat = 'identity', position = 'fill') + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ylab('% Samples') +
  scale_fill_brewer(palette = 'Blues')
ggsave(filename = 'Figures/Supp_SMCbrcaClassification.pdf', plot = g_brca,
       height = 4, width = 4)

ann_smc.plot2 <- ann_smc %>%
  group_by(Subtype, HRDgroup) %>%
  summarise(n = n())

g_subtype <- ggplot(ann_smc.plot2, aes(x = Subtype, y = n, fill = HRDgroup)) +
  geom_bar(stat = 'identity', position = 'fill') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ylab('% Samples') +
  scale_fill_brewer(palette = 'Blues')
ggsave(filename = 'Figures/Supp_SMCsubtypeClassification.pdf', plot = g_subtype,
       height = 4, width = 5)

