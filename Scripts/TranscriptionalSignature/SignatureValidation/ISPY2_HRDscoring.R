#####
## Correlate HRD transcriptional score with pCR status in I-SPY2 trial
#####

# Load libraries
library(ggpubr)

## Organise Data

# Load I-SPY2 expression and clinical data from Puzstai et al. 2021
ispy2_expr <- read.table('~/Data/ClinicalData/ISPY2_Puzstai2021_expression.txt',h=T,row.names = 1)

ispy2_response <- read.csv('~/Data/ClinicalData/GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv')
ispy2_response <- ispy2_response[ispy2_response$Arm == 'durvalumab/olaparib', ]
ispy2_response$pCR.status[ispy2_response$pCR.status == -1] <- 0 # present in control arm

# Extract clinical arm from expression data
ispy2_response$ResearchID <- paste0('X',ispy2_response$ResearchID)
ispy2_expr <- ispy2_expr[,ispy2_response$ResearchID]

# Load signature centroids, extract ElasticNet_alpha0.25, and organise it with expression matrix
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
sig <- signature.centroid.list$ElasticNet_alpha0.25

genes.intersect <- intersect(rownames(sig), rownames(ispy2_expr))
sig <- sig[genes.intersect, ]
ispy2_expr <- ispy2_expr[genes.intersect, ]

## Calculate HRD scores and match with relevant clinical data
ispy2_hrd <- data.frame(
  ResearchID = colnames(ispy2_expr),
  HRD_score = apply(ispy2_expr, 2, function(x) cor(x,sig$HRD)) -
    apply(ispy2_expr, 2, function(x) cor(x,sig$HR_proficient))
)

ispy2_hrd <- merge(x = ispy2_hrd, y = ispy2_response[,c('ResearchID','pCR.status', 'PARPi7_sig.')])

# Format pCR status into responders vs non-responders
ispy2_hrd$Response <- sapply(ispy2_hrd$pCR.status, function(x)
  ifelse(x == 1, 'Responder', 'Non-responder'))
ispy2_hrd$Response <- factor(ispy2_hrd$Response,
                             levels = c('Non-responder','Responder'))

# Plot HRD score against pCR status
g_ispy2 <- ggboxplot(data = ispy2_hrd, x = 'Response', y = 'HRD_score',
          add = 'jitter', color = 'Response') +
  stat_compare_means() + scale_color_brewer(palette = 'Paired') +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank()) +
  ylab('HRD score')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure5/ISPY2_Response.pdf',
       plot = g_ispy2, width = 5, height = 5)
# ggsave(filename = '~/Projects/Thesis/Chapter 5/ISPY2_controlArm.pdf',
#        plot = g_ispy2, width = 5, height = 5)

# Plot HRD score against PARPi7 score
g_hrdvsparpi7 <- ggplot(data = ispy2_hrd, aes(x = PARPi7_sig., y = HRD_score)) +
  geom_point(aes(color = Response)) + 
  geom_smooth(method = 'lm', color = 'gray60') + stat_cor() +
  xlab('PARPi7 Signature Score') + ylab('HRD score') +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = 'top') +
  scale_color_brewer(palette = 'Paired')
# ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure5/ISPY2_HRDvsPARPi7.pdf',
#        plot = g_hrdvsparpi7, width = 5, height = 5)

# Plot PARPi7 score against pCR status
g_parpi7 <- ggboxplot(data = ispy2_hrd, x = 'Response', y = 'PARPi7_sig.',
                     add = 'jitter', color = 'Response') +
  stat_compare_means() + scale_color_brewer(palette = 'Paired') +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank()) +
  ylab('PARPi7 Signature Score')
# ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_ISPY2_PARPi7.pdf',
#        plot = g_parpi7, width = 5, height = 5)
