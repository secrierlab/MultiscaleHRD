#####
## Run Gene Set Enrichment Analysis on the 130-gene signature using the pathfindR package
#####

library(pathfindR)

load('~/Projects/HRD_TranscriptionalSignature/Data/TCGA_BRCA.ExprDeconvolution_050823_p0.79_training.Rdata')
rownames(Z.tumor_training) <- substr(rownames(Z.tumor_training), 1, 12)
Z.tumor_training <- log2(Z.tumor_training + 1)

load('~/Projects/HRD_MutationalSignature/Results/TCGA_HRDclassification_BRCAannotation.Rdata')
rownames(ann_tcga) <- ann_tcga$Patient
ann_tcga <- ann_tcga[,c('HRD', 'BRCA_status')]
ann_tcga <- ann_tcga[!is.na(ann_tcga$BRCA_status), ]
ann_tcga <- ann_tcga[!(ann_tcga$BRCA_status %in% c('PALB2', 'RAD51C')), ]

ann_tcga$group <- ann_tcga$BRCA_status
ann_tcga$group[ann_tcga$HRD == 'HRD' & ann_tcga$BRCA_status == 'none'] <- 'HRD_BRCA+'
# ann_tcga$group[ann_tcga$group == 'none'] <- 'HR-proficient'
ann_tcga$group[ann_tcga$HRD == 'HR-proficient'] <- 'HR-proficient'
ann_tcga$group <- factor(ann_tcga$group, levels = c('HR-proficient', 'HRD_BRCA+',
                                                    'BRCA1', 'BRCA2'))

samples.intersect <- intersect(rownames(Z.tumor_training), rownames(ann_tcga))

ann_tcga <- ann_tcga[samples.intersect, ]
Z.tumor_training <- Z.tumor_training[samples.intersect, ]

# Extract relevant signature
load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
centroid.brca <- signature.centroid.list$ElasticNet_alpha0.25

input.x = Z.tumor_training[,rownames(centroid.brca)]

# For each gene, run an ANOVA against the four HRD/BRCA-defect groups
#   Save the p-values and adjust accordingly
df.res <- data.frame(
  Gene.symbol = colnames(input.x),
  pVal = apply(input.x, 2, function(x)
    summary(aov(x ~ ann_tcga$BRCA_status))[[1]]$`Pr(>F)`[1])
)
df.res$adj.P.Val = p.adjust(df.res$pVal, method = 'BH')
df.res <- df.res[,c(1,3)]

# Run pathfindR and save results
df.res_output <- run_pathfindR(df.res, p_val_threshold = 1,
                               gene_sets = 'GO-BP',
                               min_gset_size = 10,
                               output_dir = '~/Projects/HRD_TranscriptionalSignature/Results/pathfindR')

# Hierarchical clustering of enriched terms and extract representative clusters
res_clustered <- cluster_enriched_terms(
  df.res_output, plot_dend = FALSE, plot_clusters_graph = FALSE)
res_clustered <- res_clustered[res_clustered$Status == 'Representative', ]

pdf('~/Projects/HRD_TranscriptionalSignature/Figures/Supp_GSEApathfindR.pdf')
enrichment_chart(res_clustered)
dev.off()
