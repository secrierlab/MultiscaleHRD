#####
## Check gene dropout in Qian et al 2020 BC dataset
#####

setwd('~/Data/scRNASeq/Qian2020/')

# Load libraries
library(dplyr)
library(ggplot2)
library(anndata)

# Load expression and metadata and subset for cancer cells
load('exprData_Qian2020.Rdata')

meta.qian <- read.csv('2103-Breastcancer_metadata.csv', header = TRUE)
meta.qian <- meta.qian[match(colnames(expr.data_qian2020), meta.qian$Cell), ]
expr.cancer <- expr.data_qian2020[,meta.qian$CellType == 'Cancer']

# Load signature centroids
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')

# load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_alternatives.Rdata')
# names(signature_alternative.centroid.list) <- paste0('Alternative_Sig_',names(signature_alternative.centroid.list))
# 
# signature.centroid.list <- c(signature.centroid.list, signature_alternative.centroid.list)
# rm(signature_alternative.centroid.list)

stats.df <- data.frame(
  Model = names(signature.centroid.list),
  prop_nonZeroCells = NA, 
  median_GenesExpressedPerCell = NA, mean_GenesExpressedPerCell = NA
)

for (i in 1:length(signature.centroid.list)) {
  
  print(names(signature.centroid.list)[i])
  
  sig.i <- signature.centroid.list[[i]]
  
  genes.intersect <- intersect(rownames(sig.i), rownames(expr.cancer))
  
  sig.i <- sig.i[genes.intersect, ]
  expr.cancer.i <- expr.cancer[genes.intersect, ]
  
  # Collate relevant results
  stats.df$prop_nonZeroCells[i] <- mean(apply(expr.cancer.i, 2, sum) > 0)
  stats.df$median_GenesExpressedPerCell[i] <- median(apply(expr.cancer.i, 2, function(x) sum(x>0)))
  stats.df$mean_GenesExpressedPerCell[i] <- mean(apply(expr.cancer.i, 2, function(x) sum(x>0)))
  
  
}

# stats.df$group <- sapply(stats.df$Model, function(x) strsplit(x,split='_')[[1]][1])

library(ggplot2)
library(tidyr)

# g1 <- ggplot(stats.df, aes(x = Model, y = prop_nonZeroCells, fill = Model)) + geom_bar(stat ='identity') + theme(axis.text.x = element_blank())
# g2 <- ggplot(stats.df, aes(x = Model, y = median_GenesExpressedPerCell, fill = Model)) + geom_bar(stat ='identity') + theme(axis.text.x = element_blank())
# g3 <- ggplot(stats.df, aes(x = Model, y = mean_GenesExpressedPerCell, fill = Model)) + geom_bar(stat ='identity') + theme(axis.text.x = element_blank())
# 
# ggarrange(plotlist = list(g1,g2,g3))

stats.df_plot <- stats.df %>%
  pivot_longer(cols = -Model, names_to = 'Measure')
g_dropout <- ggplot(stats.df_plot, aes(x = Model, y = value, fill = Model)) + 
  geom_bar(stat = 'identity') + theme_minimal() +
  theme(legend.position = 'top',
        axis.text.x = element_blank()) +
  facet_wrap(~Measure, scales = 'free')
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_QianSignatureComparison.pdf',
       plot = g_dropout, width = 8, height = 4)
