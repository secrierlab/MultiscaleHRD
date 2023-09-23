expr <- read.csv('../Downloads/TcgaBrca_exprHRDsig.csv', row.names = 1)
expr <- as.data.frame(t(expr))

# Load important genes
imp.df <- read.csv('Data/imp_score_avg.csv')

genes.imp_hrd <- imp.df$gene_names[imp.df$HRD > 0.7]
genes.imp_hrprof <- imp.df$gene_names[imp.df$HR.proficient > 0.7]

genes.imp <- c(genes.imp_hrd, genes.imp_hrprof)
genes.imp_group <- c(rep('HRD', length(genes.imp_hrd)),
                     rep('HR-prof', length(genes.imp_hrprof)))

expr.redux <- expr[,genes.imp]

library(WGCNA)
library(tidyr)

adj <- as.data.frame(adjacency(expr.redux))
diag(adj) <- NA
adj$source <- rownames(adj)
adj$group <- genes.imp_group
adj[adj < 0.0002] <- NA

adj <- adj %>%
  pivot_longer(cols = -c(source, group), names_to = 'target', values_to = 'int')
adj <- adj[!is.na(adj$int), ]

adj$node1 <- apply(adj[,c('source','target')], 1, function(x) sort(x)[1])
adj$node2 <- apply(adj[,c('source','target')], 1, function(x) sort(x)[2])
adj$link <- apply(adj[,c('node1','node2')], 1, function(x) paste(x,collapse='_'))
adj <- adj[!duplicated(adj$link), ]

write.csv(adj[,c(1,3,4)], file = '~/Projects/HRD_TranscriptionalSignature/Results/adjacency_reduxSig.csv')

adj.df <- data.frame(
  Gene = genes.imp,
  Group = genes.imp_group
)
# adj.df <- adj.df[adj.df$name %in% c(adj$source,adj$target), ]
write.csv(adj.df, file = '~/Projects/HRD_TranscriptionalSignature/Results/adjInfo_reduxSig.csv')
