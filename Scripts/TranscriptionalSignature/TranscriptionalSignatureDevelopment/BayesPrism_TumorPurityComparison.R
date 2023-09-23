library(TCGAbiolinks)
library(ggpubr)

load('~/Data/TCGA/TCGA_BRCA.BayesPrism.theta.Rdata')
df.theta <- data.frame(
  Sample.ID = sapply(rownames(theta), function(x) substr(x, 1, 16)),
  CancerCellFraction = theta[,'Cancer']
)

df.theta <- merge(x = df.theta, y = Tumor.purity[,c('Sample.ID','ESTIMATE')])

# Function for converting ESTIMATE column in Tumor.purity data.frame
comma.function <- function(x) {
  x.tmp = as.character(x)
  x.tmp = strsplit(x.tmp,split=',')[[1]]
  x.final = as.numeric(paste(x.tmp,collapse='.'))
  return(x.final)
}

df.theta$TumorPurity <- sapply(df.theta$ESTIMATE, comma.function)

# Compare estimates
g_tumorPurity <- ggplot(df.theta, aes(x = TumorPurity, y = CancerCellFraction)) +
  geom_point(alpha = 0.3) + geom_smooth(method = 'lm') + stat_cor() +
  theme_minimal()
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Supp_BayesPrismEstimates.pdf',
       plot = g_tumorPurity)
