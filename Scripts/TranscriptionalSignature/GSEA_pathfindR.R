#####
## Run Gene Set Enrichment Analysis on the 130-gene signature using the pathfindR package
#####

setwd('~/Documents/GitHub/HRD_classification/')

# Load libraries
library(pathfindR)

# Load in training data and 130-gene signature, and extract gene signature from data
load('Data/TCGA_BRCA/TCGA_BRCA_FPKM_filterByExpr_TPregress_training.Rdata')
load('Results/TranscriptionalSignature/Templates/HRD/template_HRD_ElasticNet.Rdata')

input.x <- input_data.train[, rownames(template.hrd)]
input.y <- input_data.train$BRCA_defect

# For each gene, run an ANOVA against the four HRD/BRCA-defect groups
#   Save the p-values and adjust accordingly
df.res <- data.frame(
  Gene.symbol = colnames(input.x),
  pVal = apply(input.x, 2, function(x) 
    summary(aov(x ~ input.y))[[1]]$`Pr(>F)`[1])
)
df.res$adj.P.Val = p.adjust(df.res$pVal)
df.res <- df.res[,c(1,3)]

# Run pathfindR and save results
df.res_output <- run_pathfindR(df.res, p_val_threshold = 1,
                               gene_sets = 'Reactome',
                               min_gset_size = 10,
                               output_dir = 'Results/TranscriptionalSignature/pathfindR')

# Hierarchical clustering of enriched terms and extract representative clusters
res_clustered <- cluster_enriched_terms(
  df.res_output, plot_dend = FALSE, plot_clusters_graph = FALSE)
res_clustered <- res_clustered[res_clustered$Status == 'Representative', ]

pdf('Figures/Supplementary/pathfindR_RepresentativeTerms.pdf',
    width = 9, height = 4)
enrichment_chart(res_clustered)
dev.off()
