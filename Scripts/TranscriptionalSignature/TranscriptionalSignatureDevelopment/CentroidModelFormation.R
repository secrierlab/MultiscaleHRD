# Collate centroid models for Runs

# Load deconvoluted training data for reference
load('~/Projects/HRD_TranscriptionalSignature/Data/TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.Rdata')
rownames(Z.tumor_training) <- sapply(rownames(Z.tumor_training), function(x) substr(x,1,12))
Z.tumor_training <- log2(Z.tumor_training + 1)
samples <- rownames(Z.tumor_training)
Z.tumor_training <- apply(Z.tumor_training, 2, scale)
rownames(Z.tumor_training) <- samples

load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')
# ann_tcga$HRD <- sapply(ann_tcga$HRD_prob, function(x) ifelse(x >= 0.32, 'HRD', 'HR-proficient'))
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('HRD', 'BRCA_status')]
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('PALB2', 'RAD51C')), ]

ann_tcga$group <- ann_tcga$BRCA_status
ann_tcga$group[ann_tcga$HRD == 'HRD' & ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga$group[ann_tcga$group == 'none'] <- 'HR-proficient'
# ann_tcga$group[ann_tcga$HRD == 'HR-proficient'] <- 'HR-proficient'
ann_tcga$group <- factor(ann_tcga$group, levels = c('HR-proficient', 'HRD_BRCA+',
                                                    'BRCA1', 'BRCA2'))

samples.intersect <- intersect(rownames(Z.tumor_training), rownames(ann_tcga))

ann_tcga_train <- ann_tcga[samples.intersect, ]
Z.tumor_training <- Z.tumor_training[samples.intersect, ]

# Z.tumor_training <- apply(Z.tumor_training, 2, scale)
# rownames(Z.tumor_training) <- rownames(ann_tcga_train)

# # Apply to non-deconvoluted samples
# library(TCGAbiolinks)
# library(SummarizedExperiment)
# 
# setwd('~/Data/TCGA')
# query <- GDCquery(
#   project = 'TCGA-BRCA',
#   data.category = 'Transcriptome Profiling',
#   data.type = 'Gene Expression Quantification',
#   workflow.type = 'STAR - Counts',
#   barcode = rownames(Z.tumor_training)
# )
# # GDCdownload(query)
# expr.train <- GDCprepare(query = query)
# expr.train <- expr.train[,expr.train$sample_type == 'Primary Tumor']
# expr.train <- expr.train[!duplicated(rowData(expr.train)$gene_name) &
#                          !is.na(rowData(expr.train)$gene_name), ]
# 
# library(SummarizedExperiment)
# expr.tumor_training <- assay(expr.train, 'fpkm_uq_unstrand')
# rownames(expr.tumor_training) <- rowData(expr.train)$gene_name
# colnames(expr.tumor_training) <- sapply(colnames(expr.tumor_training),
#                                        function(x) substr(x,1,12))
# expr.tumor_training <- log2(expr.tumor_training+1)
# 
# expr.tumor_training <- t(expr.tumor_training)
# expr.tumor_training <- expr.tumor_training[rownames(ann_tcga_train), ]

# Move to CV coefficients and generate centroids
setwd('~/Projects/HRD_TranscriptionalSignature/CV_Coefficients/Run6/')
signature.centroid.list <- list()

models <- list.files()
models <- sapply(models, function(x) paste0(strsplit(x,split='_')[[1]][2:3],collapse='_'))
models <- sapply(models, function(x) substr(x,6,nchar(x)))
models <- unique(models)

for (model in models) {
  
  print(model)
  
  files.coefs <- list.files(pattern = model)
  coefs_join <- list()
  for (file in files.coefs) {
    load(file)
    coefs_join <- c(coefs_join, coefs)
    rm(coefs)
  }
  
  coef.mat <- matrix(NA, nrow = nrow(coefs_join[[1]]$`HR-proficient`), ncol = length(coefs_join))
  rownames(coef.mat) <- rownames(coefs_join[[1]]$`HR-proficient`)
  for (i in 1:length(coefs_join)) {
    coef.mat[,i] <- coefs_join[[i]]$`HR-proficient`[,1] != 0
  }
  coef.mat <- coef.mat[-1,]
  genes.include <- rownames(coef.mat)[apply(coef.mat, 1, sum) == 1000]
  
  print(paste0('Number of genes in model ', model, ': ', length(genes.include)))
  
  ## Create centroids
  centroid.model <- data.frame(
    HRD = apply(Z.tumor_training[ann_tcga_train$HRD == 'HRD',genes.include], 2, mean),
    HR_proficient = apply(Z.tumor_training[ann_tcga_train$HRD == 'HR-proficient',genes.include],2,mean),
    BRCA1 = apply(Z.tumor_training[ann_tcga_train$group == 'BRCA1',genes.include],2,mean),
    BRCA2 = apply(Z.tumor_training[ann_tcga_train$group == 'BRCA2',genes.include],2,mean),
    HRD_BRCApos = apply(Z.tumor_training[ann_tcga_train$group == 'HRD_BRCA+',genes.include],2,mean),
    HR_BRCA_proficient = apply(Z.tumor_training[ann_tcga_train$group == 'HR-proficient',genes.include],2,mean)
  )
  
  signature.centroid.list[[model]] <- centroid.model
  
}

save(signature.centroid.list, file = '~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')




