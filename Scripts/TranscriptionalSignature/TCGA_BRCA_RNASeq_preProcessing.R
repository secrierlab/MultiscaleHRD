#####
## TCGA-BRCA RNA-SEQ DATA PRE-PROCESSING
#####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(edgeR)
library(SummarizedExperiment)


# 1. Download data (loaded in by TCGAbiolinks)
load('~/Documents/PhD/Data/TCGA_BRCA_RNASeq_FPKM.Rdata')


# 2. Filter for lowly expressed genes
#   This is achieved using the FilterByExpr() function in edgeR
#   This strategy keeps genes with > k counts in > n samples
#   Here, we keep genes with > 0.2 counts in > 70% (default) samples

index.keep <- filterByExpr(assay(se.brca), min.count = 0.2) # Keep 16,500 features, exclude 39,854
se.brca.filtered <- se.brca # save index.keep until the end



# 3. Regress filtered expression against tumour purity
#   This step aims to remove TME features
#   Tumour purity data obtained from GDC Data Portal
#   Before regression, expression values are log2-normalised

TCGA_purity <- read.table('~/Documents/GitHub/HRD_classification/Data/TCGA_BRCA_purity.txt', header = TRUE)
TCGA_purity$Sample.ID <- sapply(TCGA_purity$sample, function(x)
  paste0(strsplit(x, split='-')[[1]][1:4], collapse='-'))
TCGA_purity <- TCGA_purity[, c('Sample.ID', 'purity')]
TCGA_purity <- TCGA_purity[!is.na(TCGA_purity$purity), ]

# Match with se.brca.filtered
patients.intersect <- intersect(TCGA_purity$Sample.ID, se.brca.filtered$sample)
tp.brca <- TCGA_purity[match(patients.intersect, TCGA_purity$Sample.ID), ]
se.brca.filtered <- se.brca.filtered[, match(patients.intersect, se.brca.filtered$sample)]
assay(se.brca.filtered) <- log2(assay(se.brca.filtered) + 1)

# Regress against tumour purity
se <- assay(se.brca.filtered)
se <- apply(se, 1, function(x)
  lm(x ~ tp.brca$purity)$residuals)

se.brca.purity.filtered <- se.brca.filtered
assay(se.brca.purity.filtered) <- t(se)


# 4. Split into training/testing data
#   Load in results of mutational signature classifier
#   Add HRD/HRP and BRCA-HRD labels from HRD annotations
load('Results/ExomeClassifier/TCGA_BRCA/TCGA_HRD_annotation.Rdata')
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]

ann_tcga$HRD_BRCA_label <- as.character(ann_tcga$BRCA_status)
ann_tcga$HRD_BRCA_label[ann_tcga$HRD_BRCA_label == 'none' &
                          ann_tcga$HRD == 'HRD'] <- 'HRD_BRCA+'
ann_tcga$HRD_BRCA_label[ann_tcga$HRD_BRCA_label == 'none'] <- 'HR-proficient'

# Find overlapping patients and organise RNA-seq data
samples.intersect <- intersect(rownames(ann_tcga),
                               se.brca.purity.filtered$patient)
ann_tcga <- ann_tcga[samples.intersect, ]

se.brca.purity.filtered <- se.brca.purity.filtered[, match(samples.intersect, se.brca.purity.filtered$patient)]
index.keep <- index.keep[!duplicated(rowData(se.brca.purity.filtered)$external_gene_name)]

se.brca.purity.filtered <- se.brca.purity.filtered[!duplicated(rowData(se.brca.purity.filtered)$external_gene_name), ]
se <- assay(se.brca.purity.filtered)
rownames(se) <- rowData(se.brca.purity.filtered)$external_gene_name
se <- t(se)

# Add in BRCA defect HRD status
input_data <- as.data.frame(se)

input_data$HRD <- ann_tcga$HRD
input_data$BRCA_defect <- ann_tcga$HRD_BRCA_label

# 5. Separate into training (2/3) and testing (1/3)

input_data <- input_data[order(rownames(input_data)), ]
set.seed(1234)
training_index <- rbinom(n = nrow(input_data), size = 1, prob = 2/3) == 1

input_data.train_full <- input_data[training_index, ] # for templated formation of alternative signatures
input_data.train <- input_data[training_index, c(index.keep, TRUE, TRUE)]
input_data.test <- input_data[!training_index, ]

# Save training and testing data
save(input_data.train_full,
     file = 'Data/TCGA_BRCA/TCGA_BRCA_FPKM_TPregressONLY_training.Rdata')
save(input_data.train,
     file = 'Data/TCGA_BRCA/TCGA_BRCA_FPKM_filterByExpr_TPregress_training.Rdata')
save(input_data.test,
     file = 'Data/TCGA_BRCA/TCGA_BRCA_FPKM_filterByExpr_TPregress_testing.Rdata')
