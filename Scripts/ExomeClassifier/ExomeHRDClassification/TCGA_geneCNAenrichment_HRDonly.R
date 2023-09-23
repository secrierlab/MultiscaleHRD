#####
## Fisher's tests to determine differential chromosome Gene alterations in HRD vs HR-proficient samples
#####

setwd('~/Projects/HRD_MutationalSignature/')

# Load libraries
library(ggplot2)
library(wesanderson)
library(ggrepel)

# Load Gene-level CNA data, and organise to remove metadata
cna <- read.delim('~/Data/TCGA/brca_tcga_pan_can_atlas_2018/data_cna.txt')
cna <- cna[!duplicated(cna$Hugo_Symbol), ]

# Extract cancer genes with alterations in >5% cases
cgc_data <- read.csv('~/Data/Census_allMon Jul  3 15_56_40 2023.csv')
cgc_data.genes <- as.character(cgc_data$Gene.Symbol)

cna <- cna[cna$Hugo_Symbol %in% cgc_data.genes, ]

# Sort dataframe
rownames(cna) <- cna$Hugo_Symbol
cna <- cna[,-c(1,2)]
cna <- as.data.frame(t(cna))
cna$Patient <- sapply(rownames(cna), function(x)
  paste(strsplit(x,split='[.]')[[1]][1:3], collapse='-'))

# Load HRD/HR-proficiency labels and merge with CNA data
load('Results/TCGA_HRDclassification_BRCAannotation.Rdata')
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[ann_tcga$HRD == 'HRD', ]
ann_tcga$group <- sapply(ann_tcga$BRCA_status,
                         function(x) ifelse(x == 'none', 'BRCA+', 'BRCA-defective'))

df.ann <- data.frame(
  Patient = ann_tcga$Patient,
  group = ann_tcga$group
)

df <- merge(x = df.ann, y = cna)

# Initialise data frames to track gain/loss enrichments
fishers.gain = fishers.loss <- data.frame(
  Gene = names(df)[-c(1,2)], Estimate = NA, pVal = NA
)

# For each chromosome Gene:
#   Conduct a Fisher's exact test for HRD and HR-proficient samples against:
#     1) Gain vs not-Gain
#     2) Loss vs not-Loss
#   Save odds ratios and p-values in initialised datasets
for (i in 1:nrow(fishers.gain)) {
  Gene.i <- fishers.gain$Gene[i]
  
  df.i <- data.frame(
    group = df$group,
    Gain = df[,Gene.i] > 0,
    Loss = df[,Gene.i] < 0
  )
  df.i <- df.i[!is.na(df.i[,'Gain']), ]
  
  df.i$group <- factor(df.i$group, levels = c('BRCA-defective', 'BRCA+'))
  df.i$Gain <- factor(df.i$Gain, levels = c(FALSE,TRUE))
  df.i$Loss <- factor(df.i$Loss, levels = c(FALSE,TRUE))
  
  # Gains
  table.gain.i <- table(df.i$group, df.i$Gain)
  fishers.gain.i <- fisher.test(table.gain.i)
  
  fishers.gain$Estimate[i] <- fishers.gain.i$estimate
  fishers.gain$pVal[i] <- fishers.gain.i$p.value
  
  # Losses
  table.loss.i <- table(df.i$group, df.i$Loss)
  fishers.loss.i <- fisher.test(table.loss.i)
  
  fishers.loss$Estimate[i] <- fishers.loss.i$estimate
  fishers.loss$pVal[i] <- fishers.loss.i$p.value
  
}

## Add information to gain/loss results:
#   - Adjust p-values
#   - If Estimate == 'Inf', set to maximum estimate
#   - If Estimate == 0, set to minimum estimate
#   - log2-normalise estimates and -log10(p-adjust)
#   - Add labels for Genes with enrichment significance < .01
#   - 'Volcano plot' displaying results

# Gains
fishers.gain$padj <- p.adjust(fishers.gain$pVal)
# fishers.gain$padj <- fishers.gain$pVal
fishers.gain$Estimate[fishers.gain$Estimate == 'Inf'] <- max(fishers.gain$Estimate)
fishers.gain$Estimate[fishers.gain$Estimate == 0] <- min(fishers.gain$Estimate)
fishers.gain$l2fc <- log2(fishers.gain$Estimate)
fishers.gain$logp <- -log10(fishers.gain$padj)
fishers.gain$label <- fishers.gain$Gene
fishers.gain$label[fishers.gain$padj > .05] <- NA

g_gains <- ggplot(fishers.gain, aes(x = l2fc, y = logp, label = label)) + 
  geom_point() + theme_minimal() + geom_text_repel() +
  ggtitle('Gain vs !Gain')
# ggsave(filename = 'Figures/Supplementary/ChrGeneEnrich_GainVsNotGain.pdf',
#        plot = g_gains)

# Losses
fishers.loss$padj <- p.adjust(fishers.loss$pVal)
# fishers.loss$padj <- fishers.loss$pVal
fishers.loss$Estimate[fishers.loss$Estimate == 'Inf'] <- max(fishers.loss$Estimate)
fishers.loss$Estimate[fishers.loss$Estimate == 0] <- min(fishers.loss$Estimate)
fishers.loss$l2fc <- log2(fishers.loss$Estimate)
fishers.loss$logp <- -log10(fishers.loss$padj)
fishers.loss$label <- fishers.loss$Gene
fishers.loss$label[fishers.loss$padj > .05] <- NA

g_losses <- ggplot(fishers.loss, aes(x = l2fc, y = logp, label = label)) + 
  geom_point() + theme_minimal() + geom_text_repel() +
  ggtitle('Loss vs !Loss')
# ggsave(filename = 'Figures/Supplementary/ChrGeneEnrich_LossVsNotLoss.pdf',
#        plot = g_losses)

# Merge log-fold changes (+ labels, showing significant Genes for each test)
df.fishers <- merge(x = fishers.gain[,c('Gene','l2fc','label')],
                    y = fishers.loss[,c('Gene','l2fc','label')], by = 'Gene')
names(df.fishers) <- c('Gene','l2fc_GAIN', 'label_GAIN','l2fc_LOSS','label_LOSS')

# Add GAIN/LOSS/GAIN+LOSS labels for which tests are significant
df.fishers$Label <- 'none'
df.fishers$Label[!is.na(df.fishers$label_GAIN)] <- 'GAIN'
df.fishers$Label[!is.na(df.fishers$label_LOSS)] <- 'LOSS'
df.fishers$Label[!is.na(df.fishers$label_GAIN) & !is.na(df.fishers$label_LOSS)] <- 'GAIN+LOSS'
df.fishers$Gene_Label <- NA
df.fishers$Gene_Label[df.fishers$Label != 'none'] <- df.fishers$Gene[df.fishers$Label != 'none']

# df.fishers$Label <- factor(df.fishers$Label, levels = c('GAIN','GAIN+LOSS','LOSS','none'))

g_GainvsLoss <- ggplot(df.fishers, aes(x = l2fc_GAIN, y = l2fc_LOSS, color = Label, label = Gene_Label)) +
  geom_point() + geom_text_repel() + theme_minimal() +
  scale_color_manual(values = c('gray90')) +
  xlab('log2(Fold Change) - Gain') + ylab('log2(Fold Change) - Loss') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'top',
        legend.title = element_blank())
ggsave(filename = 'Figures/Supp_GeneCnaEnrich_GainVsLoss_HRDonly.pdf',
       plot = g_GainvsLoss, height = 4, width = 4)
