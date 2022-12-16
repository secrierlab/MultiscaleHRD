#####
## Multinomial Elastic Net Regression to classify BRCA1-/-, BRCA2-/-, HRD, HR-proficient TCGA-BRCA samples
#####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(glmnet)

# Load in training data:
#   TCGA-BRCA, FPKM + log2 normalised
#   filtered away lowly expressed genes, regressed against tumour purity
load('Data/TCGA_BRCA/TCGA_BRCA_FPKM_filterByExpr_TPregress_training.Rdata')

# Separate out expression data from BRCA_defect status
input.x <- as.matrix(input_data.train[,1:(ncol(input_data.train)-2)])
input.y <- input_data.train$BRCA_defect

# Conduct 500 iterations of 10-fold cross validation:
#   (500 iteration ~ 2 days)
#   For each iteration, calculate coefficients using cv$lambda.min
#     (lambda coefficient which gives smallest error)
#   set seed inside loop

iterations <- 1:500
coefs <- list()

for (i in iterations) {
  
  # Set seed with each iteration
  set.seed(123*i)
  
  print(paste0('Iteration ', i, ' of ', iterations[length(iterations)], '...'))
  
  # Conduct grouped multinomial 10-fold cross validation
  cv <- cv.glmnet(x = input.x, y = input.y, family = 'multinomial',
                  alpha = 0.5, type.multinomial = 'grouped')
  coefs[[i]] <- coef(cv$glmnet.fit, s = cv$lambda.min)
  
}

save(coefs, file = 'Results/TranscriptionalSignature/Coefficients/CVcoefficients_ElasticNet_iter1to500_min.Rdata')
