#####
## Establishing training data and running BayesPrism for Expression Deconvolution of TCGA-BRCA
#####

setwd('~/Daniel/TCGAExpressionDeconvolution/')

# Load libraries
library(Seurat)
library(dplyr)
library(SummarizedExperiment)
library(caret)
library(BayesPrism)

# 1. Load Qian et al. 2020 data and create matrix

bc.data <- Read10X(data.dir = 'Data/Qian2020/BC_counts/')
bc_sc <- as.matrix(bc.data)
bc_sc <- t(bc_sc)

# 2. Match cell.type.labels from Qian et al. metadata
bc_met <- read.csv('Data/Qian2020/2103-Breastcancer_metadata.csv.gz')
bc_cellMeta <- bc_met$CellType

# 3a. Load TCGA bulk data obtained from TCGAbiolinks

load('Data/TCGA_BRCA.counts.SE.Rdata')
tcga_brca.rnaseq <- tcga_brca.rnaseq[!duplicated(rowData(tcga_brca.rnaseq)$gene_name), ]
tcga_brca.rnaseq <- tcga_brca.rnaseq[,!duplicated(tcga_brca.rnaseq$patient)]

# 3b. Match to HRD/BRCA-status output and remove PALB2/RAD51C defects
load('Data/TCGA_HRDclassification_BRCAannotation.Rdata')
# ann_tcga$patient <- rownames(ann_tcga)
# ann_tcga$HRD <- sapply(ann_tcga$HRD_prob, function(x) ifelse(x >= 0.79, 'HRD', 'HR-proficient'))
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga$HRD_BRCAstatus <- ann_tcga$BRCA_status
ann_tcga$HRD_BRCAstatus[ann_tcga$HRD == 'HRD' &
                          ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga$HRD_BRCAstatus[ann_tcga$HRD_BRCAstatus == 'none'] <- 'HR-proficient'

ann_tcga <- ann_tcga[!(ann_tcga$HRD_BRCAstatus %in% c('PALB2', 'RAD51C')), ]

# 3c. Separate into training and testing (preserving class proportions), and save both
patients.intersect <- intersect(ann_tcga$Patient, tcga_brca.rnaseq$patient)
ann_tcga <- ann_tcga[match(patients.intersect, ann_tcga$patient), ]
tcga_brca.rnaseq <- tcga_brca.rnaseq[, match(patients.intersect, tcga_brca.rnaseq$patient)]

set.seed(1234)
inTrain <- createDataPartition(y = ann_tcga$HRD_BRCAstatus, p = 2/3, list = FALSE)
inTrain <- 1:ncol(tcga_brca.rnaseq) %in% inTrain
tcga_brca.rnaseq_train <- tcga_brca.rnaseq[, inTrain]
tcga_brca.rnaseq_test <- tcga_brca.rnaseq[, !inTrain]

tcga_brca.rnaseq_train$HRD_BRCAstatus <- ann_tcga$HRD_BRCAstatus[inTrain]
tcga_brca.rnaseq_test$HRD_BRCAstatus <- ann_tcga$HRD_BRCAstatus[!inTrain]

save(tcga_brca.rnaseq_train, file = 'Data/TCGA_BRCA.counts.SE_050823_training.Rdata')
save(tcga_brca.rnaseq_test, file = 'Data/TCGA_BRCA.counts.SE_050823_testing.Rdata')

bc_bulk <- t(assay(tcga_brca.rnaseq_train))
colnames(bc_bulk) <- rowData(tcga_brca.rnaseq_train)$gene_name

# 4. Preprocessing of Qian et al.

# Plot single cell outliers
sc.stat <- plot.scRNA.outlier(
  input = bc_sc,
  cell.type.labels = bc_met$CellType,
  species = 'hs',
  return.raw = FALSE,
  pdf.prefix = 'Qian2020_outlierPlot'
)

# Filter genes expressed in <2% cancer cells
bc_sc_cancer <- bc_sc[bc_cellMeta == 'Cancer', ]
index.gene_propCancer <- apply(bc_sc_cancer, 2, function(x) mean(x>0))
bc_sc <- bc_sc[,index.gene_propCancer > 0.02]

# Filter outlier genes from scRNA-seq data
bc_sc.filtered <- cleanup.genes(
  input = bc_sc,
  input.type = 'count.matrix',
  species = 'hs',
  gene.group = c('Rb', 'Mrp', 'other_Rb', 'chrM', 'MALAT1', 'chrX', 'chrY')
)
dim(bc_sc.filtered)

# Check concordance of varying gene types between scRNA-seq and bulk
plot.bulk.vs.sc(sc.input = bc_sc.filtered,
                bulk.input = bc_bulk)

# Subset protein coding genes, since these are most concordant and 
#   computation can be sped up
bc_sc.filtered.pc <- select.gene.type(
  input = bc_sc.filtered,
  gene.type = 'protein_coding'
)

# Construct a prism object
myPrism <- new.prism(
  reference = bc_sc.filtered.pc,
  mixture = bc_bulk,
  input.type = 'count.matrix',
  cell.type.labels = bc_cellMeta,
  cell.state.labels = bc_cellMeta,
  key = 'Cancer',
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

# Run BayesPrism
Sys.time(); tcga_brca.bayesPrism <- run.prism(prism = myPrism, n.cores = 20); Sys.time()

save(tcga_brca.bayesPrism, file = 'TCGA_BRCA.BayesPrism_050823_p0.79_training.Rdata')

# Extract cancer-specific expression
Z.tumor_training <- get.exp(
  bp = tcga_brca.bayesPrism, 
  state.or.type = 'type', cell.name = 'Cancer'
)
save(Z.tumor_training, file = 'TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.Rdata')

