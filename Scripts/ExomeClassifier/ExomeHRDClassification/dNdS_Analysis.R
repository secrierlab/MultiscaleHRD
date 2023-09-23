#####
## dN/dS analysis to identify mutations under positive selection in HRD/BRCA groups
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(TCGAbiolinks)
library(maftools)
library(dndscv)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(wesanderson)

# Load data and dN/dS references
load('~/Data/RefCDS_human_GRCh38.p12_dNdScv.0.1.0.rda')

setwd('~/Data/TCGA')
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Simple Nucleotide Variation',
  data.type = 'Masked Somatic Mutation'
)
# GDCdownload(query)
mutations <- GDCprepare(query = query)

data <- read.maf(maf = mutations, isTCGA = TRUE)
data <- rbind(data@data, data@maf.silent)
dat_small <- data[,c('Tumor_Sample_Barcode', 'Chromosome',
                     'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2')]
names(dat_small) <- c('sampleID', 'chr', 'pos', 'ref', 'mut')

# Load HRD/BRCA groups and arrange data accordingly (note, samples can overlap)
setwd('~/Projects/HRD_MutationalSignature/')

load('Results/TCGA_HRDclassification_BRCAannotation.Rdata')
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ] # remove samples with missing BRCA_status
# ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('RAD51C','PALB2')), ] # remove samples with deficiencies in other HR genes
dat_small <- merge(x = dat_small, y = ann_tcga[,c('Patient', 'HRD', 'BRCA_status')],
                   by.x = 'sampleID', by.y = 'Patient')

# HRD vs HR-proficient
data_HRD <- dat_small[dat_small$HRD == 'HRD', 1:5]
data_HRproficient <- dat_small[dat_small$HRD == 'HR-proficient', 1:5]

# BRCA-defect categories
data_BRCA1 <- dat_small[dat_small$BRCA_status == 'BRCA1', 1:5]
data_BRCA2 <- dat_small[dat_small$BRCA_status == 'BRCA2', 1:5]
data_RAD51C <- dat_small[dat_small$BRCA_status == 'RAD51C', 1:5]
data_HRDBRCApos <- dat_small[dat_small$HRD == 'HRD' &
                               dat_small$BRCA_status == 'none', 1:5]

## Run dN/dS analysis on each group separately

# BRCA1
dndsoutBRCA1 <- dndscv(data_BRCA1, cv = NULL,
                       refdb = RefCDS)
sel_cvBRCA1 <- dndsoutBRCA1$sel_cv

# BRCA2
dndsoutBRCA2 <- dndscv(data_BRCA2, cv = NULL,
                       refdb = RefCDS)
sel_cvBRCA2 <- dndsoutBRCA2$sel_cv

# RAD51C
dndsoutRAD51C <- dndscv(data_RAD51C, cv = NULL,
                        refdb = RefCDS)
sel_cvRAD51C <- dndsoutRAD51C$sel_cv

# HRD_BRCA+
dndsoutHRDBRCApos <- dndscv(data_HRDBRCApos, cv = NULL,
                            refdb = RefCDS)
sel_cvHRDBRCApos <- dndsoutHRDBRCApos$sel_cv

# HRD full
dndsoutHRD <- dndscv(data_HRD, cv = NULL,
                     refdb = RefCDS)
sel_cvHRD <- dndsoutHRD$sel_cv

# HR-proficient
dndsoutHRprof <- dndscv(data_HRproficient, cv = NULL,
                        refdb = RefCDS)
sel_cvHRprof <- dndsoutHRprof$sel_cv


## Plot HRD vs HR-proficient gene selections

# Find genes under positive selection in at least one group
sig_genes <- unique(c(sel_cvHRD$gene_name[sel_cvHRD$qind_cv < .1],
                      sel_cvHRprof$gene_name[sel_cvHRprof$qind_cv < .1]))

sel_cvSigGenes <- merge(
  x = sel_cvHRD[sel_cvHRD$gene_name %in% sig_genes, c('gene_name', 'wind_cv', 'qind_cv')],
  y = sel_cvHRprof[sel_cvHRprof$gene_name %in% sig_genes, c('gene_name', 'wind_cv', 'qind_cv')],
  by = 'gene_name'
)
names(sel_cvSigGenes)[-1] <- c('dNdS_HRD', 'FDR_HRD',
                               'dNdS_HRprof', 'FDR_HRprof')

# Prepare labels for plotting
sel_cvSigGenes$Significant <- 'Both'
sel_cvSigGenes$Significant[sel_cvSigGenes$FDR_HRD > .1] <- 'HR-proficient ONLY'
sel_cvSigGenes$Significant[sel_cvSigGenes$FDR_HRprof > .1] <- 'HRD ONLY'

sel_cvSigGenes$dNdS_HRD <- log2(sel_cvSigGenes$dNdS_HRD + 1)
sel_cvSigGenes$dNdS_HRprof <- log2(sel_cvSigGenes$dNdS_HRprof + 1)

# save(sel_cvSigGenes, file = 'Results/ExomeClassifier/TCGA_BRCA/dNdS_HRDvsHRprof.Rdata')

ggplot(sel_cvSigGenes, aes(x = dNdS_HRD, y = dNdS_HRprof, 
                           color = Significant, label = gene_name)) +
  geom_point() + theme_minimal() + geom_text_repel()

  
# Plot dN/dS of indsense variants across all groups
sel_cvFull <- rbind(sel_cvHRprof, sel_cvHRD, sel_cvBRCA1,
                    sel_cvBRCA2, sel_cvRAD51C, sel_cvHRDBRCApos)
sel_cvFull <- sel_cvFull[,c('gene_name','wmis_cv','qmis_cv')]
sel_cvFull$group <- rep(c('HR-proficient','HRD full','BRCA1',
                          'BRCA2','RAD51C','HRD BRCA+'),
                        each = nrow(sel_cvFull)/6)

sel_cvFull <- sel_cvFull[sel_cvFull$qmis_cv < 0.05, ]

sel_cvFull$log_dNdS <- log2(sel_cvFull$wmis_cv)
sel_cvFull$log_FDR <- -log10(sel_cvFull$qmis_cv)
sel_cvFull$log_FDR[sel_cvFull$log_FDR == 'Inf' | sel_cvFull$log_FDR >= 9] <- 9

sel_cvFull$group <- factor(sel_cvFull$group,
                           levels = c('HR-proficient','HRD full','BRCA1',
                                      'BRCA2','RAD51C','HRD BRCA+'))

g_dnds <- ggballoonplot(sel_cvFull, x = 'gene_name', y = 'group',
              size = 'log_dNdS', fill = 'log_FDR') +
  scale_fill_continuous(low = 'lightblue', high = 'navy') +
  geom_hline(yintercept = 2.5, col = 'red', linetype = 'dashed') +
  theme(legend.position = 'top')
# ggsave(filename = 'Figures/Figure2/HRD_dNdSballoonPlot.pdf',
#        plot = g_dnds, width = 4, height = 7)
ggsave(filename = '~/Projects/Thesis/Chapter 4/HRD_dNdSballoonPlot.pdf',
       plot = g_dnds, width = 7, height = 3)
