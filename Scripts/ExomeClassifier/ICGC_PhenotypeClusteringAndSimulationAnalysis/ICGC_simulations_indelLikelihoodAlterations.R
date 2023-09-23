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

pheatmap(sigs_complete, cutree_rows = 3, method = 'average',
         show_rownames = FALSE)

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

## 4. WRITE SIMULATION FUNCTION

# Apply likelihood function to simulated subsampled data
simulate_likelihood_calc <- function(input_data, cluster_distributions,
                                     cluster_assign, sample_size) {
  
  # Develop simulation matrix (same dimensions as input_data)
  sims <- matrix(0, nrow = nrow(input_data), ncol = ncol(input_data))
  colnames(sims) <- colnames(input_data); rownames(sims) <- rownames(input_data)
  
  for (i in 1:nrow(input_data)) {
    s.i <- table(sample(colnames(input_data), size = sample_size,
                        replace = TRUE, prob = input_data[i,]/sum(input_data)))
    sims[i, names(s.i)] <- s.i
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

## 5. Save likelihood distributions with altered indel contributions

# Generate a list of likelihood distributions with indel contribution * indel_alter
#   of which indel_alter = seq(0.2, 2, by=.2)

cluster_distribution.list <- list()
indel_alterations <- c(1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5)
for (i in 1:length(indel_alterations)) {
  mut.prob_sbs <- mut.prob[,grepl(pattern = '>', colnames(mut.prob))]
  mut.prob_id <- mut.prob[, grepl(pattern = ':', colnames(mut.prob))]
  
  indel_alter <- indel_alterations[i]
  
  mut.prob_id <- mut.prob_id * indel_alter
  for (j in 1:3) {
    mut.prob_sbs[j,] <- mut.prob_sbs[j,]*(1-sum(mut.prob_id[j,]))/sum(mut.prob_sbs[j,])
  }
  
  mut.prob_altered <- cbind(mut.prob_sbs, mut.prob_id)
  cluster_distribution.list[[i]] <- mut.prob_altered
  
}

## 6. COMPLETE SIMULATION RUNNING

# This function runs n simulations, returning a data frame of AUCs and plottable posteriors
run_likelihood_sims <- function(input_data, sample_size, 
                                cluster_distributions_index, cluster_assign,
                                n_simulations) {
  
  cluster_distributions <- cluster_distribution.list[[cluster_distributions_index]]
  
  # Initialise data frames
  results.mat <- matrix(NA, nrow=0, ncol=3); colnames(results.mat) <- rownames(cluster_distributions)
  posterior.df <- data.frame()
  
  
  # Run n_simulations of simulate_likelihood_calc()
  for (i in 1:n_simulations) {
    
    set.seed(i)
    
    print(paste0('Sample size = ', sample_size,
                 ' mut.prob_index = ', indel_alterations[cluster_distributions_index],
                 ' Running simulation ', i, ' of ', n_simulations, '...'))
    
    run.i <- simulate_likelihood_calc(
      input_data, cluster_distributions, cluster_assign,
      sample_size)
    results.mat <- rbind(results.mat, run.i$auc.values)
    
    run.i$posteriors$Run <- i
    posterior.df <- rbind(posterior.df, run.i$posteriors)
    
  }
  
  # Further posterior info
  posterior.df$sample_size = sample_size
  posterior.df$indel_alter = indel_alterations[cluster_distributions_index]
  
  results <- data.frame(Pheno = rep(colnames(results.mat), each = nrow(results.mat)),
                        AUC = c(results.mat[,1], results.mat[,2], results.mat[,3]),
                        sample_size = sample_size, 
                        indel_alteration = indel_alterations[cluster_distributions_index])
  
  return(list(results = results, posteriors = posterior.df))
  
}

nSim = 100
sampleSize = 25
res.full_25 <- data.frame()

set.seed(123)
for (i in 1:length(cluster_distribution.list)) {
  
  print(paste0('Indel alteration = ', indel_alterations[i]))
  
  res.i_25 <- run_likelihood_sims(
    input_data = mut_complete, sample_size = sampleSize, 
    cluster_distributions_index = i, 
    cluster_assign = sigs.clust3$Phenotype,
    n_simulations = nSim
  )
  res.full_25 <- rbind(res.full_25, res.i_25$results)

}

# Plot results for HRD Pheno
res_25.sbs3 <- res.full_25[res.full_25$Pheno == 'SBS3', ]

my_comparisons <- list(c('0.2', '1'), c('1', '5'))

g25 <- ggboxplot(res_25.sbs3, x = 'indel_alteration', y = 'AUC', add = 'jitter') +
  stat_compare_means(comparisons = my_comparisons) + ylim(c(0.65,1.05)) +
  ggtitle('sample size = 25') +
  scale_x_discrete(labels = c('1/5', '1/4', '1/3', '1/2', '1',
                              '2', '3', '4', '5'))
# ggsave(filename = 'Figures/Supp_SimulationsIndelLikelihoods_size25.pdf', plot = g25,
#        height = 5, width = 5)
ggsave(filename = '~/Projects/Thesis/Chapter 3/ICGCsimulations_IndelWeighting_sim25.pdf', plot = g25,
       height = 7, width = 4)


# Change sample size
sampleSize = 50
res.full_50 <- data.frame()

set.seed(123)
for (i in 1:length(cluster_distribution.list)) {
  
  print(paste0('Indel alteration = ', indel_alterations[i]))
  
  res.i_50 <- run_likelihood_sims(
    input_data = mut_complete, sample_size = sampleSize, 
    cluster_distributions_index = i, 
    cluster_assign = sigs.clust3$Phenotype,
    n_simulations = nSim
  )
  res.full_50 <- rbind(res.full_50, res.i_50$results)
  
}

# Plot results for HRD Pheno
res_50.sbs3 <- res.full_50[res.full_50$Pheno == 'SBS3', ]

my_comparisons <- list(c('0.2', '1'), c('1', '5'))

g50 <- ggboxplot(res_50.sbs3, x = 'indel_alteration', y = 'AUC', add = 'jitter') +
  stat_compare_means(comparisons = my_comparisons) + ylim(c(0.65, 1.05)) +
  ggtitle('sample size = 50') +
  scale_x_discrete(labels = c('1/5', '1/4', '1/3', '1/2', '1',
                              '2', '3', '4', '5'))
# ggsave(filename = 'Figures/Supp_SimulationsIndelLikelihoods_size50.pdf', plot = g50,
#        height = 5, width = 5)
ggsave(filename = '~/Projects/Thesis/Chapter 3/ICGCsimulations_IndelWeighting_sim50.pdf', plot = g50,
       height = 7, width = 4)

# Change sample size
sampleSize = 100
res.full_100 <- data.frame()

set.seed(123)
for (i in 1:length(cluster_distribution.list)) {
  
  print(paste0('Indel alteration = ', indel_alterations[i]))
  
  res.i_100 <- run_likelihood_sims(
    input_data = mut_complete, sample_size = sampleSize, 
    cluster_distributions_index = i, 
    cluster_assign = sigs.clust3$Phenotype,
    n_simulations = nSim
  )
  res.full_100 <- rbind(res.full_100, res.i_100$results)
  
}

# Plot results for HRD Pheno
res_100.sbs3 <- res.full_100[res.full_100$Pheno == 'SBS3', ]

my_comparisons <- list(c('0.2', '1'), c('1', '5'))

g100 <- ggboxplot(res_100.sbs3, x = 'indel_alteration', y = 'AUC', add = 'jitter') +
  stat_compare_means(comparisons = my_comparisons) + ylim(c(0.65, 1.05)) +
  ggtitle('sample size = 100') +
  scale_x_discrete(labels = c('1/5', '1/4', '1/3', '1/2', '1',
                              '2', '3', '4', '5'))
# ggsave(filename = 'Figures/Supp_SimulationsIndelLikelihoods_size100.pdf', plot = g100,
#        height = 5, width = 5)
ggsave(filename = '~/Projects/Thesis/Chapter 3/ICGCsimulations_IndelWeighting_sim100.pdf', plot = g100,
       height = 7, width = 4)
