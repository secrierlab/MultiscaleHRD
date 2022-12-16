#####
## Create HRD and BRCA-defect templates from signatures and training data
#####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(pheatmap)
library(RColorBrewer)
library(wesanderson)

# Template formation:
#   - Load in coefficients from 1000 iterations of regularised regression
#   - Extract genes which appear in all 1000 iterations
#   - Subset training data for HRD/BRCA-defect
#   - Create a representative sample for each category

# Write a function which collates non-zero coefficients across all iterations
#   coefficients appear in a list
coef_process <- function(coef_list) {
  
  # Initialise data frame
  df <- matrix(NA, nrow = length(coef_list[[1]]$BRCA1),
               ncol = length(coef_list))
  
  for (i in 1:length(coef_list)) {
    df[,i] <- coef_list[[i]]$BRCA1[,1] != 0
  }
  df <- as.data.frame(df)
  rownames(df) <- rownames(coef_list[[1]]$BRCA1) # gene list
  
  return(df)
  
}

# Load in coefficients and process
load('Results/TranscriptionalSignature/Coefficients/CVcoefficients_ElasticNet_iter1to500_min.Rdata')
coef1 <- coef_process(coefs)

load('Results/TranscriptionalSignature/Coefficients/CVcoefficients_ElasticNet_iter501to1000_min.Rdata')
coef2 <- coef_process(coefs)

coef.df <- cbind(coef1, coef2); rm(coefs, coef1, coef2)

# Count number of inclusions for each gene
coef.counts <- data.frame(
  Gene = rownames(coef.df),
  Inclusions = apply(coef.df, 1, sum)
)
coef.counts <- coef.counts[-1, ] # Remove (Intercept)

genes.include <- coef.counts$Gene[coef.counts$Inclusions == 1000]

# Create template by loading training data and
#   calculating median expression for each gene in the signature
#   across the samples in each category

load('Data/TCGA_BRCA/TCGA_BRCA_FPKM_filterByExpr_TPregress_training.Rdata')

input.x <- input_data.train[, genes.include]
input.y_hrd <- input_data.train$HRD
input.y_brca <- input_data.train$BRCA_defect

# BRCA-defect template
template.brca <- data.frame(
  BRCA1 = apply(input.x[input.y_brca == 'BRCA1', ], 2, median),
  BRCA2 = apply(input.x[input.y_brca == 'BRCA2', ], 2, median),
  HRD_BRCApos = apply(input.x[input.y_brca == 'HRD_BRCA+', ], 2, median),
  HR_proficient = apply(input.x[input.y_brca == 'HR-proficient', ], 2, median)
)
save(template.brca, 
     file = 'Results/TranscriptionalSignature/Templates/BRCA_defect/template_BRCA_ElasticNet.Rdata')

# HRD template
template.hrd <- data.frame(
  HRD = apply(input.x[input.y_hrd == 'HRD', ], 2, median),
  HR_proficient = apply(input.x[input.y_hrd == 'HR-proficient', ], 2, median)
)
save(template.hrd, 
     file = 'Results/TranscriptionalSignature/Templates/HRD/template_HRD_ElasticNet.Rdata')


## Plot training data as a heatmap

# Create annotation for BRCA defect and HRD labels
ann <- data.frame(BRCA_defect = input.y_brca,
                  HRD = input.y_hrd)
rownames(ann) <- rownames(input.x)

ann_cols <- list(
  BRCA_defect = c('BRCA1' = 'blue', 'BRCA2' = 'red', 'HRD_BRCA+' = 'white', 'HR-proficient' = 'white'),
  HRD = c('HRD' = wes_palette('GrandBudapest1')[2], 'HR-proficient' = wes_palette('GrandBudapest1')[1])
)

# Set colour range
paletteLength <- 100
myColor <- colorRampPalette(c('navy', 'darkblue', 'white', 'red', 'darkred'))(paletteLength)
myBreaks <- c(seq(min(input.x), -2, length.out = ceiling(paletteLength/4)+1),
              seq(-2, 0, length.out = floor(paletteLength/4))[-1],
              seq(0, 2, length.out = floor(paletteLength/4))[-1],
              seq(2, max(input.x), length.out = floor(paletteLength/4))[-1])

# Transpose matrix so that columns are samples
input.x_plot <- t(input.x)

pheatmap(input.x_plot, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = ann, annotation_colors = ann_cols,
         color = myColor, breaks = myBreaks,
         filename = 'Figures/Figure4/TCGA_training_SignatureHeatmap.pdf')
