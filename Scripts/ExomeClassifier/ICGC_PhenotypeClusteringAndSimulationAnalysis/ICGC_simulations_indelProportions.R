#####
## Simulation analysis to determine utility of indels in HRD classification
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(pheatmap)
library(pROC)
library(ggpubr)

# Aim: demonstrate that, in subsampled cohorts, the inclusion of indels improves classification
#   In this case, they improve the probability of reclassifying SBS3-enriched samples as SBS3
#   The ICGC cohort easily clusters into 3 based on SBS signature contributions:
#     SBS2/13 (APOBEC), SBS3 (HRD), SBS5 (Ageing)

## 1. CREATE SBS MUTATIONAL SIGNATURE CLUSTERS

# Signature contributions have already been calculated. Extract SBS signatures
load('Results/ICGC_BRCA_deconstructSigs_Cutoff0.01_SBSandIDnormalised.Rdata')
sigs_complete <- sigs_complete[, grepl(pattern = 'SBS', names(sigs_complete))]

cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)
pheatmap(sigs_complete, cutree_rows = 3,
         show_rownames = FALSE, color = cols_scale,
         file = 'Figures/Supp_ICGCSignaturesSBSHeatmap.pdf')

# Create groups using hierarchical clustering
sigs.dist_mat <- dist(sigs_complete, method = 'euclidean')
sigs.hclust3 <- hclust(sigs.dist_mat, method = 'complete')
sigs.clust3 <- cutree(sigs.hclust3, k = 3)
sigs.clust3 <- data.frame(Cluster.num = sigs.clust3)
sigs.clust3$Phenotype <- sapply(sigs.clust3$Cluster.num, function(x)
  ifelse(x == 1, 'SBS5',
         ifelse(x == 2, 'SBS3', 'APOBEC')))

# Load in ICGC mutation tallies, which have already been processed
load('Data/BRCA_UKEU_mt_tally.Rdata')
mut_complete <- cbind(mt_tally.brca_wgs$SBS_96,
                      mt_tally.brca_wgs$ID_83)

## 2. GENERATE CLUSTER-SPECIFIC MUTATIONAL SPECTRA

# Separate out mut_complete based on cluster assignment
#   Then use collate_function() to generate the SBS+ID spectrum for each

mut_apobec <- mut_complete[sigs.clust3$Phenotype == 'APOBEC', ]
mut_sbs3 <- mut_complete[sigs.clust3$Phenotype == 'SBS3', ]
mut_sbs5 <- mut_complete[sigs.clust3$Phenotype == 'SBS5', ]

collate_function <- function(input, variant_type = 'ALL') {
  
  # Spectra can be generate for SBS/ID-only, or ALL
  sbs.index <- sapply(names(input), function(x) grepl('>',x))
  id.index <- sapply(names(input), function(x) grepl(':',x))
  
  if (variant_type == 'SBS') {
    input_final <- input[,sbs.index]
  } else if (variant_type == 'ID') {
    input_final <- input[,id.index]
  } else {
    input_final <- input
  }
  
  # Return average spectrum
  total = apply(input_final, 2, sum)
  dist = total/sum(total)
  return(dist)
  
}

mut.prob_apobec <- collate_function(mut_apobec)
mut.prob_sbs3 <- collate_function(mut_sbs3)
mut.prob_sbs5 <- collate_function(mut_sbs5)

mut.prob <- rbind(mut.prob_apobec, mut.prob_sbs3, mut.prob_sbs5)
rownames(mut.prob) <- c('APOBEC','SBS3','SBS5')

## 3. LIKELIHOOD FUNCTION

# This function aligns a dataset with the designated mean distribution
likelihood_calc <- function(input_data, cluster_distributions, cluster_assign) {
  
  # This function:
  #   Takes a dataset of mutations as its input (rows = samples, cols = mutation types)
  #   Limits the data to to designated mutation types
  #   Applies a likelihood approach to calculate posterior probabilities of cluster assignment
  #   NB, For now, we ensure that our cluster distributions have the correct mutation types
  #   Output: posterior distributions, note that for WES data a comparison column is not possible
  
  # Calculate likelihood: P(mutation spectrum | cluster)
  likelihoods <- matrix(NA, nrow = nrow(input_data), ncol = nrow(cluster_distributions))
  rownames(likelihoods) <- rownames(input_data); colnames(likelihoods) <- rownames(cluster_distributions)
  for (i in 1:nrow(input_data)) {
    likelihoods[i, ] <- apply(cluster_distributions, 1, 
                              function(x) prod(x ^ input_data[i, ]))
  }
  
  # Set prior: P(cluster)
  marginal.probs <- table(cluster_assign)/length(cluster_assign)
  
  # Calculate posteriors and return
  posteriors <- data.frame(t(apply(likelihoods, 1, function(x) 
    (x * marginal.probs)/sum(x * marginal.probs))))
  return(posteriors)
  
}

## 4. INDEL-BASED SUBSAMPLING

# This function applies the likelihood function to a newly generated subsample
#   This subsample consists of sample_size mutation events extracted with replacement from the original dataset
#   Of these events, indel_prop * sample_size will be indel events, and the rest are SBS (initially indel_prop = 0)

# Apply likelihood fucntion to subsampled data
simulate_likelihood_calc <- function(input_data, cluster_distributions,
                                     cluster_assign, sample_size, indel_prop = 0) {
  
  # Separate data into SBS and ID, for separate sampling
  data_sbs <- input_data[, grepl(pattern = '>', colnames(input_data))]
  data_id  <- input_data[, grepl(pattern = ':', colnames(input_data))]
  
  # Develop simulation matrix (same dimensions as input_data)
  sims <- matrix(0, nrow = nrow(input_data), ncol = ncol(input_data))
  colnames(sims) <- colnames(input_data); rownames(sims) <- rownames(input_data)
  
  for (i in 1:nrow(input_data)) {
    s_sbs <- table(sample(colnames(data_sbs), size = sample_size * (1 - indel_prop),
                          replace = TRUE, prob = data_sbs[i, ]/sum(data_sbs[i, ])))
    sims[i, names(s_sbs)] <- s_sbs
    
    s_id <- table(sample(colnames(data_id), size = sample_size * indel_prop,
                         replace = TRUE, prob = data_id[i, ]/sum(data_id[i, ])))
    sims[i, names(s_id)] <- s_id
  }
  
  # Apply likelihood function to the simulated data
  posteriors <- likelihood_calc(sims, cluster_distributions,
                                cluster_assign)
  
  # As this is simulated data, we can also add the comparisons
  posteriors$PhenoTrue <- cluster_assign
  
  # Calculate AUCs and return
  roc_apobec <- roc(posteriors$PhenoTrue == 'APOBEC', posteriors$APOBEC, quiet = TRUE)
  roc_sbs3 <- roc(posteriors$PhenoTrue == 'SBS3', posteriors$SBS3, quiet = TRUE)
  roc_sbs5 <- roc(posteriors$PhenoTrue == 'SBS5', posteriors$SBS5, quiet = TRUE)
  
  auc.values <- c(roc_apobec$auc, roc_sbs3$auc, roc_sbs5$auc)
  
  return(list(auc.values = auc.values, posteriors = posteriors))
  
}


## 5. COMPLETE SIMULATION RUNNING

# This function runs n simulations, returning a data frame of AUCs and plottable posteriors
run_likelihood_sims <- function(input_data, sample_size, indel_prop,
                                cluster_distributions, cluster_assign,
                                n_simulations) {
  
  # Initialise data frames
  results.mat <- matrix(NA, nrow=0, ncol=3); colnames(results.mat) <- rownames(cluster_distributions)
  posterior.df <- data.frame()
  
  # Run n_simulations of simulate_likelihood_calc()
  for (i in 1:n_simulations) {
    
    set.seed(i)
    
    print(paste0('Sample size = ', sample_size,
                 ' indel_prop = ', indel_prop,
                 ' Running simulation ', i, ' of ', n_simulations, '...'))
    run.i <- simulate_likelihood_calc(
      input_data, cluster_distributions, cluster_assign,
      sample_size, indel_prop)
    results.mat <- rbind(results.mat, run.i$auc.values)
    
    run.i$posteriors$Run <- i
    posterior.df <- rbind(posterior.df, run.i$posteriors)
    
  }
  
  # Further posterior info
  posterior.df$sample_size = sample_size
  posterior.df$indel_prop = indel_prop
  
  results <- data.frame(Pheno = rep(colnames(results.mat), each = nrow(results.mat)),
                        AUC = c(results.mat[,1], results.mat[,2], results.mat[,3]),
                        sample_size = sample_size, indel_prop = indel_prop)
  
  return(list(results = results, posteriors = posterior.df))
  
}

nSim = 100
res.full_25 <- res.full_50 <- res.full_100 <- data.frame()
full_posteriors <- data.frame()

for (indel_props in seq(0, 0.5, by = .05)) {
  
  res.i_25 <- run_likelihood_sims(
    input_data = mut_complete, sample_size = 25, indel_prop = indel_props,
    cluster_distributions = mut.prob, cluster_assign = sigs.clust3$Phenotype,
    n_simulations = nSim
  )
  res.full_25 <- rbind(res.full_25, res.i_25$results)
  full_posteriors <- rbind(full_posteriors, res.i_25$posteriors)
  
  print(paste0('Indel proportion = ', indel_props,
               ' sample size = 50'))
  
  res.i_50 <- run_likelihood_sims(
    input_data = mut_complete, sample_size = 50, indel_prop = indel_props,
    cluster_distributions = mut.prob, cluster_assign = sigs.clust3$Phenotype,
    n_simulations = nSim
  )
  res.full_50 <- rbind(res.full_50, res.i_50$results)
  full_posteriors <- rbind(full_posteriors, res.i_50$posteriors)
  
  print(paste0('Indel proportion = ', indel_props,
               ' sample size = 100'))
  
  res.i_100 <- run_likelihood_sims(
    input_data = mut_complete, sample_size = 100, indel_prop = indel_props,
    cluster_distributions = mut.prob, cluster_assign = sigs.clust3$Phenotype,
    n_simulations = nSim
  )
  res.full_100 <- rbind(res.full_100, res.i_100$results)
  full_posteriors <- rbind(full_posteriors, res.i_100$posteriors)
  
}

myComp <- list(c('0', '0.05'), c('0', '0.1'),
               c('0', '0.15'), c('0', '0.2'), c('0', '0.25'))

res.full_25$best_AUC <- res.full_25$indel_prop == .2
g_sim25 <- ggboxplot(res.full_25[res.full_25$Pheno == 'SBS3', ],
          x = 'indel_prop', y = 'AUC', add = 'jitter', color = 'best_AUC') +
  xlab('Indel Proportions') +
  geom_vline(xintercept = 1+(1/0.05)*0.0675, linetype = 'dashed', color = 'blue') +
  guides(color = 'none') + scale_color_manual(values = c('black','red')) +
  stat_compare_means(comparisons = myComp) + ggtitle('sample size = 25')
# ggsave(filename = 'Figures/Supp_ICGCsimulations_indelProps_sim25.pdf', plot = g_sim25,
#        width = 8, height = 4)
ggsave(filename = '~/Projects/Thesis/Chapter 3/ICGCsimulations_indelProps_sim25.pdf', plot = g_sim25,
       width = 12, height = 4)

g_sim50 <- ggboxplot(res.full_50[res.full_50$Pheno == 'SBS3', ],
          x = 'indel_prop', y = 'AUC', add = 'jitter') +
  geom_vline(xintercept = 1+(1/0.05)*0.0675, linetype = 'dashed', color = 'red') +
  stat_compare_means(comparisons = myComp)
# ggsave(filename = 'Figures/Figure1/ICGCsimulations_indelProps_sim50.pdf', plot = g_sim50,
#        width = 8, height = 4)

res.full_50$best_AUC <- res.full_50$indel_prop == .1
g_sim50 <- ggboxplot(res.full_50[res.full_50$Pheno == 'SBS3', ],
                     x = 'indel_prop', y = 'AUC', add = 'jitter', color = 'best_AUC') +
  xlab('Indel Proportions') +
  geom_vline(xintercept = 1+(1/0.05)*0.0675, linetype = 'dashed', color = 'blue') +
  guides(color = 'none') + scale_color_manual(values = c('black','red')) +
  stat_compare_means(comparisons = myComp) + ggtitle('sample size = 50')
# ggsave(filename = 'Figures/Supp_ICGCsimulations_indelProps_sim50.pdf', plot = g_sim50,
#        width = 8, height = 4)
ggsave(filename = '~/Projects/Thesis/Chapter 3/ICGCsimulations_indelProps_sim50.pdf', plot = g_sim50,
       width = 12, height = 4)

res.full_100$best_AUC <- res.full_100$indel_prop == .05
g_sim100 <- ggboxplot(res.full_100[res.full_100$Pheno == 'SBS3', ],
          x = 'indel_prop', y = 'AUC', add = 'jitter', color = 'best_AUC') +
  xlab('Indel Proportions') +
  geom_vline(xintercept = 1+(1/0.05)*0.0675, linetype = 'dashed', color = 'blue') +
  guides(color = 'none') + scale_color_manual(values = c('black','red')) +
  stat_compare_means(comparisons = myComp) + ggtitle('sample size = 100')
# ggsave(filename = 'Figures/Supp_ICGCsimulations_indelProps_sim100.pdf', plot = g_sim100,
#        width = 8, height = 4)
ggsave(filename = '~/Projects/Thesis/Chapter 3/ICGCsimulations_indelProps_sim100.pdf', plot = g_sim100,
       width = 12, height = 4)
