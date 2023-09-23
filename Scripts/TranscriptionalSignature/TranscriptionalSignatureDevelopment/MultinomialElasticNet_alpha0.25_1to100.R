#####
## Multinomial Elastic Net Regression to classify BRCA1-/-, BRCA2-/-, HRD, HR-proficient TCGA-BRCA samples
#####

setwd('~/TranscriptionalSignatures/Revisions')

# Load libraries
library(glmnet)

# Load in training data:
#   TCGA-BRCA, expression deconvolution by BayesPrism, should be followed by low-gene removal and log2-normalisation
#   Genes remain if they have >1 count in >70% samples
#   HRD/BRCA status
load('Data/TCGA_BRCA.ExprDeconvolution_270723_p0.79_training.Rdata')
rownames(Z.tumor_training) <- substr(rownames(Z.tumor_training), 1, 12)
#index.keep <- apply(Z.tumor_training, 2, function(x) sum(x > 1)/length(x) > 0.7)
#Z.tumor_training <- Z.tumor_training[,index.keep]
#Z.tumor_training <- log2(Z.tumor_training + 1)

# Organise output variable
load('Data/TCGA_HRDclassification_BRCAannotation_HRD0.79.Rdata')
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('HRD', 'BRCA_status')]
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('PALB2', 'RAD51C')), ]

# Match output to training data
samples.intersect <- intersect(rownames(Z.tumor_training), rownames(ann_tcga))

ann_tcga <- ann_tcga[samples.intersect, ]
ann_tcga$group <- ann_tcga$BRCA_status
ann_tcga$group[ann_tcga$HRD == 'HRD' & ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga$group[ann_tcga$group == 'none'] <- 'HR-proficient'
ann_tcga$group <- factor(ann_tcga$group, levels = c('HR-proficient', 'HRD_BRCA+',
                                                    'BRCA1', 'BRCA2'))

Z.tumor_training <- Z.tumor_training[samples.intersect, ]


# Separate out expression data from BRCA_defect status
input.x <- Z.tumor_training
input.y <- ann_tcga$group

# Conduct 500 iterations of 10-fold cross validation:
#   (500 iteration ~ 2 days)
#   For each iteration, calculate coefficients using cv$lambda.min
#     (lambda coefficient which gives smallest error)
#   set seed inside loop

iterations <- 1:100
coefs <- list()

for (i in iterations) {
  
  # Set seed with each iteration
  set.seed(123*i)
  
  print(paste0('Iteration ', i, ' of ', iterations[length(iterations)], '...', Sys.time()))
  
  # Conduct grouped multinomial 10-fold cross validation
  cv <- cv.glmnet(x = input.x, y = input.y, family = 'multinomial',
                  alpha = 0.25, type.multinomial = 'grouped')
  coefs[[i]] <- coef(cv$glmnet.fit, s = cv$lambda.min)
  
}

save(coefs, file = 'Results/CVcoefficients_MultiElasticNet_alpha0.25_p0.79_iter1to100.Rdata')
