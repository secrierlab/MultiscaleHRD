#####
## Simulation analysis to determine overall subsampling reclassification in ICGC-BRCA cohort
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(pROC)
library(ggpubr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# Complete SBS/indel counts, prior probabilities, and cluster likelihood spectra
load('Data/BRCA_UKEU_mt_tally.Rdata')
mut_complete <- as.data.frame(cbind(mt_tally.brca_wgs$SBS_96, mt_tally.brca_wgs$ID_83))

load('Data/ClusterLikelihoods/ICGC_clust20_mclust_meanCont.Rdata')

load('Results/ICGC_BRCA_IDnormalised_PhenotypeAnnotation_clust20.Rdata')
mut_complete <- mut_complete[rownames(ann), ]
pheno_assigned <- ann$Phenotype


# Write a likelihood function that aligns a dataset with the designated mean distribution
likelihood_calc <- function(input_data, cluster_distributions, cluster_assign) {
  
  # This function:
  #   Takes a dataset of mutations as its input (rows = samples, cols = mutation types)
  #   Limits the data to to designated mutation types
  #   Applies a likelihood approach to calculate posterior probabilities of cluster assignment
  #   NB, For now, we ensure that our cluster distributions have the correct mutation types
  #   Output: posterior distributions, note that for WES data a comparison column is not possible
  
  # Calculate likelihood: P(mutation spectrum | cluster)
  log_likelihoods <- matrix(NA, nrow = nrow(input_data), ncol = nrow(cluster_distributions))
  rownames(log_likelihoods) <- rownames(input_data); colnames(log_likelihoods) <- rownames(cluster_distributions)
  for (i in 1:nrow(input_data)) {
    # print(i)
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
  Sys.time(); for (i in 1:ncol(final_probs)) {
    final_probs[,i] <- 10^(apply(log_posteriors, 1, function(x) x[i] - max(x[1:ncol(final_probs)])))
  }
  final_probs <- data.frame(t(apply(final_probs, 1, function(x) x/sum(x)))); Sys.time()
  
  return(final_probs)
  
}



# Apply likelihood function to simulated data
simulate_likelihood_calc <- function(input_data, cluster_distributions, 
                                     cluster_assign, mutation_types = 'ALL', sample_size) {
  
  # input_data: data frame of samples (rows) and 96 trinucloetide contexts (cols)
  # sample_size: number of mutations to be sampled from patient with replacement
  
  # Limit input data to designated mutation types
  index.sbs <- sapply(names(input_data), function(x) grepl(pattern = '>', x))
  index.id <- sapply(names(input_data), function(x) grepl(pattern = ':', x)) 
  
  if (mutation_types == 'SBS') input_data <- input_data[, index.sbs]
  if (mutation_types == 'ID') input_data <- input_data[, index.id]
  
  # Develop simulation matrix (same dimensions as input_data, but with sampled data)
  sims <- matrix(0, nrow = nrow(input_data), ncol = ncol(input_data))
  colnames(sims) <- colnames(input_data); rownames(sims) <- rownames(input_data)
  
  for (i in 1:nrow(input_data)) {
    s <- table(sample(colnames(input_data), size = sample_size, replace = TRUE,
                      prob = input_data[i, ]/sum(input_data[i, ])))
    sims[i,names(s)] <- s
  }
  
  # Apply likelihood function to the simulated data
  posteriors <- likelihood_calc(sims, cluster_distributions, 
                                cluster_assign)
  
  # As this is simulated data, we can add the comparisons
  posteriors$PhenoTrue <- cluster_assign
  
  # Save and return AUC values
  auc.df <- data.frame(Phenotype = unique(posteriors$PhenoTrue),
                       AUC = NA)
  
  for (pheno in posteriors$PhenoTrue) {
    roc.full <- roc(posteriors$PhenoTrue == pheno, posteriors[,pheno], quiet = TRUE)
    auc.val <- roc.full$auc
    auc.df$AUC[auc.df$Phenotype == pheno] <- auc.val
  }
  
  return(list(auc.df = auc.df,
              posteriors = posteriors))
  
}


# Write function to run n simulations and return a data frame of AUC values
run_likelihood_sims <- function(input_data, sample_size, 
                                cluster_distributions, cluster_assign, 
                                mutation_types = 'ALL', n_simulations) {
  
  # All of the same inputs as simulate_likelihood_calc() + n_simulations
  #   Again, this function is about comparing simulated data with true results
  
  # Initialise matrix
  results.mat <- matrix(NA, nrow = 0, ncol = nrow(cluster_distributions)); colnames(results.mat) <- rownames(cluster_distributions)
  
  # Save final posterior distributions in list
  posteriors_list <- list()
  
  # Run n_simulations of simulate_likelihood_calc():
  for (i in 1:n_simulations) {
    print(paste0('Running simulation ', i, ' of ', n_simulations, '...'), quote = FALSE)
    post.i <- simulate_likelihood_calc(input_data = input_data, sample_size = sample_size,
                                       cluster_distributions = cluster_distributions, cluster_assign = cluster_assign,
                                       mutation_types = mutation_types)
    posteriors_list[[i]] <- post.i$posteriors[,1:(ncol(post.i$posteriors)-1)]
    
    auc.i <- post.i$auc.df$AUC
    results.mat <- rbind(results.mat, auc.i)
  }
  
  results <- data.frame(Pheno = rep(colnames(results.mat), each = nrow(results.mat)),
                        AUC = as.vector(results.mat))
  
  # Return collated results and final posterior distribution
  # return(list(results = results, posterior_n = post.i$posteriors))
  return(list(results = results, posteriors_list = posteriors_list))
  
}


# Run likelihood simulations and save/plot output
set.seed(123)
sig.type = 'ALL' # one of c('ALL','SBS','ID)
n_simulations = 100

for (sampleSize in c(25, 50, 100)) {
  
  print(paste0('Running simulations: Sample Size = ', sampleSize, '...'), quote = FALSE)
  
  # Run simulation
  res.simSampleSize <- run_likelihood_sims(input_data = mut_complete, sample_size = sampleSize,
                                           cluster_distributions = mut.dists_mean, cluster_assign = pheno_assigned,
                                           mutation_types = sig.type, n_simulations = n_simulations)
  res.simSampleSize$results$sample_size = sampleSize
  
  # Separate results and plot AUCs
  res_results <- res.simSampleSize$results
  write.table(res_results, file = paste0('Results/ICGC_simulations_AUCs_sims', sampleSize, '.txt'),
              quote = FALSE, sep = '\t', row.names = FALSE)
  
  res_results$Group <- sapply(as.character(res_results$Pheno),
                              function(x) strsplit(x,split='_')[[1]][1])
  res_results$Group[grepl(pattern = 'ID', res_results$Group)] <- 'ID_enriched'
  res_results$Pheno <- factor(res_results$Pheno,
                              levels = names(table(res_results$Pheno))[length(unique(res_results$Pheno)):1])
  
  g_auc <- ggboxplot(data = res_results, x = 'Pheno', y = 'AUC',
                     color = 'Group', orientation = 'horizontal')
  ggsave(filename = paste0('Figures/Supp_ICGCsimulations_PhenoReassign_sim', sampleSize, '_AUCs.pdf'), plot = g_auc)
         
  # Plot overall posterior distributions
  res_totalPosterior <- do.call(rbind, res.simSampleSize$posteriors_list)
  res_totalPosterior$Pheno_Assigned <- apply(res_totalPosterior, 1, function(x) names(x)[x==max(x)])
  res_totalPosterior$Pheno_True <- rep(ann$Phenotype, n_simulations)
  
  res_totalPosterior_summary <- as.matrix(table(res_totalPosterior$Pheno_True, res_totalPosterior$Pheno_Assigned))
  res_totalPosterior_summary <- apply(res_totalPosterior_summary, 1, function(x) x/sum(x))
  
  cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)
  pheatmap(res_totalPosterior_summary, color = cols_scale,
           cluster_rows = FALSE, cluster_cols = FALSE,
           filename = paste0('Figures/Supp_ICGCsimulations_PhenoReassign_sim', sampleSize, '_posteriorHeatmap.pdf'))
  
}

