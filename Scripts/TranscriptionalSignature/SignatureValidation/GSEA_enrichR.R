library(enrichR)

dbs <- listEnrichrDbs()
dbs <- dbs$libraryName[c(212:214,173)]

load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
enriched <- enrichr(rownames(signature.centroid.list$ElasticNet_alpha0.25),
                    dbs)

pdf('~/Projects/HRD_TranscriptionalSignature/Figures/Supp_GSEAenrichR.pdf')
plotEnrich(enriched$KEGG_2021_Human, showTerms = 15)
dev.off()
