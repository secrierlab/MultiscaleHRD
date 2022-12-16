##### 
## Creation of signature phenotypes in 614 ICGC-BRCA samples
#####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(mclust)
library(pheatmap)
library(readxl)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

# Load ICGC deconstructSigs data
load('Data/ICGC_BRCA/ICGC_BRCA_deconstructSigs_Cutoff05.Rdata')

# Run mixture modelling using mclust
sigs.BIC <- mclustBIC(sigs_complete, G = 2:20)
summary(sigs.BIC) # Top clusters: VEI, 19, 20, 16

pdf(file = 'Figures/Supplementary/ICGC_Signatures_MixtureModelling.pdf',
    width = 5, height = 5)
plot(sigs.BIC) # save this plot
dev.off()

# Apply optimal clustering in sigs.BIC
mod_sigs.BIC <- Mclust(sigs_complete, x = sigs.BIC)
table(mod_sigs.BIC$classification) # classification

## Form annotation document for visualisation

# CHORD
chord <- read_excel('~/Downloads/41467_2020_19406_MOESM4_ESM.xlsx', sheet = 'CHORD')
chord <- chord[!grepl(pattern = 'HMF', chord$group), ]
chord <- chord[,c('sample','response','hr_status','hrd_type')]
names(chord)[3:4] <- c('CHORD', 'CHORD_type')

# HRDetect
hrdetect.samples <- read_excel('~/Downloads/41591_2017_BFnm4292_MOESM9_ESM.xlsx', sheet = 'Data', skip = 2)
hrdetect.pred <- read_excel('~/Downloads/41591_2017_BFnm4292_MOESM11_ESM.xlsx', sheet = 'b.Predictor', skip = 2)
hrdetect.pred <- hrdetect.pred[match(hrdetect.samples$Sample, hrdetect.pred$sample), ]

hrdetect <- merge(x = hrdetect.samples[,c(1,3:8)],
                  y = hrdetect.pred[,c('sample','predictorProb')],
                  by.x = 'Sample', by.y = 'sample')
hrdetect$Gene[hrdetect$isBrcaMonoallelic] <- NA
hrdetect <- hrdetect[,c('Sample','ER status','Gene','predictorProb')]
names(hrdetect) <- c('sample','ER_status','BRCA_defect', 'HRDetect')
hrdetect$sample <- sapply(hrdetect$sample,
                          function(x) substr(x, 1, nchar(x)-1)) # remove last letter

# HRDetect validation: sensitiity = 76/77 = 98.7%
table(hrdetect$BRCA_defect, hrdetect$HRDetect > .7,
      dnn = c('BRCA_defect', 'score > .7'), useNA = 'always')

# Load ICGC sample data (for sampleID matching)
samples.eu <- read.table('~/Downloads/sample.BRCA-EU.tsv', h=T, sep='\t')
samples.uk <- read.table('~/Downloads/sample.BRCA-UK.tsv', h=T, sep='\t')
samples.icgc <- rbind(samples.eu, samples.uk)
samples.icgc <- samples.icgc[,c('project_code','submitted_sample_id','icgc_donor_id')]
samples.icgc <- samples.icgc[!duplicated(samples.icgc$icgc_donor_id), ]

samples.icgc$final_letter <- sapply(samples.icgc$submitted_sample_id,
                                    function(x) substr(x,nchar(x),nchar(x)))
samples.icgc <- samples.icgc[samples.icgc$final_letter %in% c('a','b'),]
samples.icgc$sample_id <- sapply(samples.icgc$submitted_sample_id,
                                 function(x) substr(x,1,nchar(x)-1))
samples.icgc <- samples.icgc[!duplicated(samples.icgc$sample_id), ]

# Combine sample data with HRDetect and CHORD
ann.icgc <- merge(x = samples.icgc, y = hrdetect,
                  by.x = 'sample_id', by.y = 'sample', all.x = TRUE)
ann.icgc <- merge(x = ann.icgc, y = chord,
                  by.x = 'sample_id', by.y = 'sample', all.x = TRUE)
rownames(ann.icgc) <- ann.icgc$icgc_donor_id
ann.icgc <- ann.icgc[,c('ER_status', 'BRCA_defect',
                        'HRDetect', 'CHORD', 'CHORD_type')]

# Create annotation with finite mixture model clusters
ann <- data.frame(BIC_clust = factor(mod_sigs.BIC$classification,
                                     levels = 1:length(unique(mod_sigs.BIC$classification))))
ann <- merge(x = ann, y = ann.icgc,
             by = 0, all.x = TRUE)
rownames(ann) <- ann$Row.names; ann <- ann[,-1]
ann <- ann[order(ann$BIC_clust), ]

# Order samples by classification
sigs_order <- as.data.frame(t(sigs_complete[rownames(ann), ]))

# Set colours
cols <- colorRampPalette(brewer.pal(9,'Set1'))
cols_BIC <- cols(length(unique(ann$BIC_clust)))
names(cols_BIC) <- unique(ann$BIC_clust)

ann_colors = list(
  BIC_clust = cols_BIC,
  ER_status = c(positive = 'darkgreen', negative = 'yellow'),
  BRCA_defect = c(BRCA1 = 'blue', BRCA2 = 'red'),
  CHORD = c(cannot_be_determined = 'grey', HR_deficient = 'black', HR_proficient = 'white'),
  CHORD_type = c(cannot_be_determined = 'grey', BRCA1_type = 'blue', BRCA2_type = 'red', none = 'white')
)

# Use white -> navy scale
cols_scale <- colorRampPalette(colors = c('white','navy'))(1000)

pheatmap(sigs_order, cluster_cols = FALSE, show_colnames = FALSE,
         annotation_col = ann, annotation_colors = ann_colors,
         color = cols_scale)

# Naming clusters as signature phenotypes (by sight based on above heatmap)
pheno <- c('SBS5_1', 'HRD_ID8', 'HRD_ID6high', 'HRD_APOBEC',
           'SBS5_2', 'SBS5_SBS18', 'SBS5_3', 'APOBEC_ID9',
           'HRD_IDmult', 'SBS5_4', 'SBS5_5', 'APOBEC_SBS2',
           'SBS5_6', 'SBS5_7', 'APOBEC_SBS13', 'HRD_ID9',
           'HRD_ID6mid', 'ID4', 'ID2')

ann$Phenotype <- factor(pheno[as.numeric(ann$BIC_clust)],
                        levels = sort(pheno))
save(ann, file = 'Data/ICGC_BRCA/ICGC_BRCA_PhenotypeAnnotation.Rdata')

ann <- ann[order(ann$Phenotype), ]

sigs_order <- sigs_order[, rownames(ann)]

# Redo heatmap with phenotype labels
cols_Pheno <- cols(length(unique(ann$Phenotype)))
names(cols_Pheno) <- unique(ann$Phenotype)

ann_colors = list(
  Phenotype = cols_Pheno,
  ER_status = c(positive = 'darkgreen', negative = 'yellow'),
  BRCA_defect = c(BRCA1 = 'blue', BRCA2 = 'red'),
  CHORD = c(cannot_be_determined = 'grey', HR_deficient = 'black', HR_proficient = 'white'),
  CHORD_type = c(cannot_be_determined = 'grey', BRCA1_type = 'blue', BRCA2_type = 'red', none = 'white')
)

pheatmap(sigs_order, cluster_cols = FALSE, show_colnames = FALSE,
         annotation_col = ann[,c('Phenotype','BRCA_defect', 'ER_status',
                                 'HRDetect', 'CHORD', 'CHORD_type')], 
         annotation_colors = ann_colors,
         color = cols_scale, fontsize = 8, fontsize_row = 10,
         filename = 'Figures/Figure1/ICGC_Signatures_Heatmap.pdf')

## Investigate BRCA-defect distribution across HRD classifications
##   using the 'ann' data frame

ann.hrd <- ann
ann.hrd$HRD_cluster <- sapply(ann.hrd$Phenotype,
                              function(x) ifelse(grepl(pattern = 'HRD', x), 
                                                 'HRD', 'HR-proficient'))

# HRD clusters: sensitivity = 98.7%, specificity = 43.8%
table(ann.hrd$BRCA_defect, ann.hrd$HRD_cluster, 
      dnn = c('BRCA_defect', 'HRD cluster'), useNA = 'always')

ann.hrd$BRCA_defective <- ifelse(is.na(ann.hrd$BRCA_defect), 'BRCA+', 'BRCA_defective')
ann.hrd$HRD <- ifelse(grepl(pattern = 'HRD', ann.hrd$Phenotype), 'HRD', 'HR-proficient')
ann.hrd_summary <- ann.hrd %>%
  group_by(BRCA_defective, HRD) %>%
  summarise(n = n())

g_annHRDsummary <- ggplot(ann.hrd_summary, aes(x = BRCA_defective, y = n, fill = HRD)) +
  geom_bar(stat = 'identity', position = 'fill') + theme_minimal() +
  scale_fill_brewer(palette="Paired") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'top',
        legend.title = element_blank())
ggsave(filename = 'Figures/Figure1/ICGC_Signatures_HRDclassify.pdf', 
       plot = g_annHRDsummary, width = 3, height = 3)

# HRDetect > 0.7: sensitivity = 98.7%, specificity = 61.7%
table(ann.hrd$BRCA_defect, ann.hrd$HRDetect > 0.7,
      dnn = c('BRCA_defect', 'HRDetect > 0.7'), useNA = 'always')

# CHORD: sensitivity = 85.3%, specificity = 66.7%
table(ann.hrd$BRCA_defect, ann.hrd$CHORD,
      dnn = c('BRCA_defect', 'CHORD'), useNA = 'always')

## BRCA-type specific HRD classification

hrd.brca1_type <- c('HRD_APOBEC', 'HRD_ID8', 'HRD_IDmult')
hrd.brca2_type <- c('HRD_ID6mid', 'HRD_ID6high')

ann.hrd$HRD_BRCA_cluster <- 'HR-proficient'
ann.hrd$HRD_BRCA_cluster[ann.hrd$Phenotype %in% hrd.brca1_type] <- 'HRD_BRCA1type'
ann.hrd$HRD_BRCA_cluster[ann.hrd$Phenotype %in% hrd.brca2_type] <- 'HRD_BRCA2type'

# HRD signature phenotypes:
#   BRCA1: sensitivity = 88.9%, specificity = 47.1%
#   BRCA2: sensitivity = 86.7%, specificity = 54.2%
table(ann.hrd$BRCA_defect, ann.hrd$HRD_BRCA_cluster,
      dnn = c('BRCA_defect', 'HRD BRCA cluster'), useNA = 'always')

hrd_brca1type <- c('HRD_ID8', 'HRD_APOBEC', 'HRD_IDmult')
hrd_brca2type <- c('HRD_ID6high', 'HRD_ID6mid')

ann.hrd$BRCAtype_HRD <- 'HR-proficient'
ann.hrd$BRCAtype_HRD[ann.hrd$Phenotype %in% hrd_brca1type] <- 'BRCA1-type HRD'
ann.hrd$BRCAtype_HRD[ann.hrd$Phenotype %in% hrd_brca2type] <- 'BRCA2-type HRD'

ann.hrd$BRCA_defect_label <- ann.hrd$BRCA_defect
ann.hrd$BRCA_defect_label[is.na(ann.hrd$BRCA_defect_label)] <- 'BRCA+'

ann.brca_summary <- ann.hrd %>%
  group_by(BRCA_defect_label, BRCAtype_HRD) %>%
  summarise(n = n())
ann.brca_summary$BRCA_defect_label <- factor(ann.brca_summary$BRCA_defect_label,
                                             levels = c('BRCA1','BRCA2','BRCA+'))
ann.brca_summary$BRCAtype_HRD <- factor(ann.brca_summary$BRCAtype_HRD,
                                        levels = c('HR-proficient', 'BRCA2-type HRD', 'BRCA1-type HRD'))

g_annBRCAsummary <- ggplot(ann.brca_summary, aes(x = BRCA_defect_label, y = n, fill = BRCAtype_HRD)) +
  geom_bar(stat = 'identity', position = 'fill') + theme_minimal() +
  scale_fill_manual(values = c('grey','red','blue')) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(5,'mm'))
ggsave(filename = 'Figures/Figure1/ICGC_Signatures_BRCAclassify.pdf', 
       plot = g_annBRCAsummary, width = 3.3, height = 3)

# CHORD:
#   BRCA1: sensitivity = 73.3%, specificity = 55.0%
#   BRCA2: sensitivity = 93.3%, specificity = 77.8%
table(ann.hrd$BRCA_defect, ann.hrd$CHORD_type,
      dnn = c('BRCA_defect', 'CHORD'), useNA = 'always')

