#####
## Try different HRD thresholds to maximise F-score
#####

# Load libaries
library(tidyr)
library(ggplot2)
library(pROC)

# Function to calculate sensitivity/specificity/F-score

f_score <- function(t1, weight = 1) {
  recall <- t1['HRD','defective']/sum(t1[,'defective'])
  precision <- t1['HRD','defective']/sum(t1['HRD',])
  f_score <- (1 + weight^2)*(precision*recall)/((weight^2)*precision+recall)
  return(c(recall,precision,f_score))
}

hrd_thresholds <- data.frame()

# Load TCGA assignments with BRCA-defect and ER status labels
load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')

# In incremenets of 0.01, calculate F-score and add to the hrd_thresholds data frame
for (i in seq(0,1,.01)) {
  tcga_results <- ann_tcga
  tcga_results$HR_defect <- factor(ifelse(tcga_results$BRCA_status != 'none',
                                   'defective','proficient'), levels = c('defective','proficient'))
  tcga_results$HRD <- factor(ifelse(tcga_results$HRD_prob > i,
                                    'HRD','HR-proficient'), levels = c('HRD','HR-proficient'))
  
  t.i <- table(tcga_results$HRD, tcga_results$HR_defect)
  f_scores.i <- c(i, f_score(t.i,weight = 1))
  
  hrd_thresholds <- rbind(hrd_thresholds, f_scores.i)
}

names(hrd_thresholds) <- c('p_HRD', 'recall', 'precision', 'F-score')

# Plot results
hrd_thresholds.plot <- hrd_thresholds %>%
  pivot_longer(cols = -p_HRD, values_to = 'value', names_to = 'measure')
g_thres <- ggplot(hrd_thresholds.plot, aes(x = p_HRD, y = value, colour = measure)) +
  geom_line() + theme_minimal() +
  geom_vline(xintercept = hrd_thresholds$p_HRD[hrd_thresholds$`F-score` == max(hrd_thresholds$`F-score`, na.rm = TRUE)][1], 
             linetype = 'dashed') +
  xlab('p(HRD)') + theme(axis.title.y = element_blank())
ggsave(filename = '~/Projects/HRD_MutationalSignature/Figures/Supp_TCGAoptimalHRDthreshold.pdf',
       plot = g_thres, width = 5, height = 4)

# Plot AUC curve here
ann_tcga.auc <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga.auc$HR_geneDefective <- ann_tcga.auc$BRCA_status != 'none'
pdf('~/Projects/HRD_MutationalSignature/Figures/Supp_TCGAHRDthresholdAUC.pdf',
    width = 5, height = 4)
pROC_HRD <- roc(HR_geneDefective ~ HRD_prob, data = ann_tcga.auc,
                plot = TRUE, print.auc = TRUE)
dev.off()

# Repeat above analysis with subsetting for Positive and Negative ER status
er.sub <- 'Positive'
ann_tcga.er <- ann_tcga[ann_tcga$ER_status == er.sub, ]

hrd_thresholds.er <- data.frame()
for (i in seq(0,1,.01)) {
  tcga_results <- ann_tcga.er
  tcga_results$HR_defect <- factor(ifelse(tcga_results$BRCA_status != 'none',
                                          'defective','proficient'), levels = c('defective','proficient'))
  tcga_results$HRD <- factor(ifelse(tcga_results$HRD_prob > i,
                                    'HRD','HR-proficient'), levels = c('HRD','HR-proficient'))
  
  t.i <- table(tcga_results$HRD, tcga_results$HR_defect)
  f_scores.i <- c(i, f_score(t.i,weight = 1))
  
  hrd_thresholds.er <- rbind(hrd_thresholds.er, f_scores.i)
}

names(hrd_thresholds.er) <- c('p_HRD', 'recall', 'precision', 'F-score')

hrd_thresholds.er.plot <- hrd_thresholds.er %>%
  pivot_longer(cols = -p_HRD, values_to = 'value', names_to = 'measure')
g_thres.er <- ggplot(hrd_thresholds.er.plot, aes(x = p_HRD, y = value, colour = measure)) +
  geom_line() + theme_minimal() +
  geom_vline(xintercept = hrd_thresholds$p_HRD[hrd_thresholds$`F-score` == max(hrd_thresholds$`F-score`, na.rm = TRUE)][1], 
             linetype = 'dashed') +
  xlab('p(HRD)') + theme(axis.title.y = element_blank()) +
  ggtitle(paste0('ER-',er.sub, ' samples only'))
ggsave(filename = paste0('~/Projects/HRD_MutationalSignature/Figures/Supp_TCGAoptimalHRDthreshold_',er.sub,'.pdf'),
       plot = g_thres.er, width = 5, height = 4)
