#####
## Script to generate centroid templates for alternative transcriptional signatures
#####

# Libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(readxl)

# Load TCGA transcriptional data
#   Templates will be formed from FPKM-normalised training data
#   Therefore, we must also load the training data to obtain sample IDs

load('~/Projects/HRD_TranscriptionalSignature/Data/TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.Rdata')
barcodes.training <- rownames(Z.tumor_training)

setwd('~/Data/TCGA')
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts',
  barcode = barcodes.training
)
# GDCdownload(query)
expr.train <- GDCprepare(query = query)
# expr.train <- expr.train[,expr.train$sample_type == 'Primary Tumor']
expr.train <- expr.train[!duplicated(rowData(expr.train)$gene_name) &
                         !is.na(rowData(expr.train)$gene_name), ]
genes.tcga <- rowData(expr.train)$gene_name

# Get signatures
signatures_alternative <- list()

parpi7 <- c('BRCA1', 'MRE11', 'NBN', 'TDG', 'XPA', 'CHEK2', 'MAPKAPK2')
cin70 <- c('TPX2','PRC1','FOXM1','CDK1','TGIF2','MCM2','H2AZ1','TOP2A','PCNA','UBE2C',
           'MELK','TRIP13','CEP250','MCM7','RNASEH2A','RAD51AP1','KIF20A','CDC45','MAD2L1','ESPL1',
           'CCNB2','FEN1','TTK','CCT5','RFC4','ATAD2','CKAP5','NUP205','CDC20','CKS2',
           'RRM2','ELAVL1','CCNB1','RRM1','AURKB','MSH6','EZH2','CTPS1','DKC1','OIP5',
           'CDCA8','PTTG1','CEP55','H2AX','CMAS','NCAPH','MCM10','LSM4','NCAPG2','ASF1B',
           'ZWINT','PBK','ZWILCH','CDCA3','ECT2','CDC6','UNG','MTCH2','RAD21','ACTL6A',
           'PDCD2L','SRSF2','HDGF','NXT1','NEK2','DHCR7','AURKA','NDUFAB1','NEMP1','KIF4A')

signatures_alternative[['PARPi7']] <- parpi7
signatures_alternative[['CIN70']] <- cin70

# Severson
sev_init <- read_excel('../../../Downloads/BRCA1nessSignature_Severson2017.xlsx')
sev_init <- as.data.frame(sev_init)
sev <- sev_init[,1]
sev[which(!(sev %in% genes.tcga))] <- c('JAML','FAM241B','PIMREG',
                                        'HIF1A','PLAAT1','IRAG2')
signatures_alternative[['Severson']] <- sev

# Peng 2014
peng <- read_excel('../../../Downloads/HRDSignature_Peng2014.xlsx', skip=1)
peng <- as.data.frame(peng)
peng <- peng$`Gene Symbol`
peng[which(!(peng %in% genes.tcga))] <- c('MCMBP','FAM170B','DDIAS','SKA3','CEP128','TICRR',
                                          'TEDC2','METTL22','ATAD5','KIZ','ISM1','SMIM14',
                                          'SNHG32','DSCC1','DEFB1','DDX39A','HJURP','DLGAP5',
                                          'DNA2','RETREG1','H1-2','H2BC5','H2AC18','H2BC21',
                                          'HSP90AA2P','CREBRF','LOC554223','LOC649679','LOC729843','LOC91431',
                                          'VWA5A','ETFRF1','CENPU','MTARC1','BEX3','LRR1',
                                          'SRSF2','EPB41L4A-AS1','SLC35G1','TUBB4B','TUBB7P','WRAP53')
peng <- peng[peng %in% genes.tcga]
signatures_alternative[['Peng']] <- peng

save(signatures_alternative, file = '~/Projects/HRD_TranscriptionalSignature/AlternativeSignatures.Rdata')

# Collate centroids for alternative signatures
signature_alternative.centroid.list <- list()

# # Add Severson signature immediately
# rownames(sev_init) <- sev
# sev_init <- sev_init[,-1]
# colnames(sev_init) <- c('HRD','HR_proficient')
# signature_alternative.centroid.list[['Severson']] <- sev_init

# Prepare FPKM-normalised TCGA training data
expr.train_fpkm <- assay(expr.train, 'fpkm_unstrand')
rownames(expr.train_fpkm) <- rowData(expr.train)$gene_name
expr.train_fpkm <- log2(expr.train_fpkm + 1)
expr.train_fpkm <- apply(expr.train_fpkm, 1, scale)
rownames(expr.train_fpkm) <- sapply(expr.train$barcode, function(x) substr(x,1,12))
expr.train_fpkm <- as.data.frame(expr.train_fpkm)

# Add mutational signature/BRCA defect classification
load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('HRD','BRCA_status')]
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('PALB2', 'RAD51C')), ]

ann_tcga$group <- ann_tcga$BRCA_status
ann_tcga$group[ann_tcga$HRD == 'HRD' & ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
ann_tcga$group[ann_tcga$group == 'none'] <- 'HR-proficient'
# ann_tcga$group[ann_tcga$HRD == 'HR-proficient'] <- 'HR-proficient'
ann_tcga$group <- factor(ann_tcga$group, levels = c('HR-proficient', 'HRD_BRCA+',
                                                    'BRCA1', 'BRCA2'))

samples.intersect <- intersect(rownames(expr.train_fpkm), rownames(ann_tcga))

ann_tcga <- ann_tcga[samples.intersect, ]
expr.train_fpkm <- expr.train_fpkm[samples.intersect, ]

# For the remaining three signatures:
#   Generate average centroids and add to signature list
for (signature in c('Severson','PARPi7','CIN70','Peng')) {
  sig.genes <- signatures_alternative[[signature]]
  
  centroid.signature <- data.frame(
    HRD = apply(expr.train_fpkm[ann_tcga$HRD == 'HRD',sig.genes], 2, mean),
    HR_proficient = apply(expr.train_fpkm[ann_tcga$HRD == 'HR-proficient',sig.genes],2,mean),
    BRCA1 = apply(expr.train_fpkm[ann_tcga$group == 'BRCA1',sig.genes],2,mean),
    BRCA2 = apply(expr.train_fpkm[ann_tcga$group == 'BRCA2',sig.genes],2,mean),
    HRD_BRCApos = apply(expr.train_fpkm[ann_tcga$group == 'HRD_BRCA+',sig.genes],2,mean),
    HR_BRCA_proficient = apply(expr.train_fpkm[ann_tcga$group == 'HR-proficient',sig.genes],2,mean)
  )
  
  if (signature == 'Severson') {
    centroid.signature$HRD = sev_init$`BRCA1ness template Pearson correlations`
    centroid.signature$HR_proficient = sev_init$`non-BRCAness template Pearson correlations`
  }
  
  signature_alternative.centroid.list[[signature]] <- centroid.signature
  
}

save(signature_alternative.centroid.list, file = '~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives_zNorm.Rdata')
