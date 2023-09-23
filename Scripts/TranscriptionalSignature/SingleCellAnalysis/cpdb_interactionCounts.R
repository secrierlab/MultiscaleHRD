setwd('~/Projects/HRD_TranscriptionalSignature/CellphoneDB/Results/')

library(tidyr)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(networkD3)
library(wesanderson)

qian <- read.delim('Bassez2021/significant_means.txt')
names(qian) <- gsub('HR.proficient','HR_proficient',names(qian))

# target_cancer.index <- which(sapply(names(qian), function(x) strsplit(x,split='[.]')[[1]][2])
#       %in% c('CancerHR_proficient','CancerHRD'))
qian <- qian[,c(2,13:ncol(qian))]

qian <- qian %>%
  pivot_longer(cols = -interacting_pair, names_to = 'cell_int', values_to = 'interaction')
qian$SOURCE <- sapply(qian$cell_int, function(x) strsplit(x,split='[.]')[[1]][1])
qian$TARGET <- sapply(qian$cell_int, function(x) strsplit(x,split='[.]')[[1]][2])

qian <- qian %>%
  group_by(SOURCE,TARGET) %>%
  summarise(interactions = sum(!is.na(interaction)))

qian_heatmap <- qian %>%
  pivot_wider(names_from = TARGET, values_from = interactions)
qian_heatmap <- as.data.frame(qian_heatmap)
rownames(qian_heatmap) <- qian_heatmap$SOURCE; qian_heatmap <- qian_heatmap[,-1]

col1 = 'dodgerblue4'; col2 = 'peachpuff'; col3 = 'deeppink4'
col.heatmap <- colorRampPalette(c(col1, col2, col3))(1000)
pheatmap(t(qian_heatmap), color = col.heatmap)

# cellType_order <- c('Fibroblast','EC','Myeloid','DC','Mast',
#                     'T_cell','B_cell','CancerHR_proficient','CancerHRD')
cellType_order <- c('Fibroblast','Endothelial_cell','Myeloid_cell','pDC','Mast_cell',
                    'T_cell','B_cell','Cancer_cellHR_proficient','Cancer_cellHRD')

qian_heatmap <- qian_heatmap[cellType_order, cellType_order]
pheatmap(t(qian_heatmap), color = col.heatmap,
         cluster_rows = FALSE, cluster_cols = FALSE,
         filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure7/Bassez_Heatmap.pdf',
         width = 5, height = 5)

qian_sankey <- qian[!grepl(pattern = 'Cancer',qian$SOURCE) &
                      grepl(pattern = 'Cancer',qian$TARGET), ]
qian_sankey.Nodes <- data.frame(
  name = c(as.character(qian_sankey$SOURCE),
           as.character(qian_sankey$TARGET)) %>% unique()
)
qian_sankey$IDsource <- match(qian_sankey$SOURCE, qian_sankey.Nodes$name)-1
qian_sankey$IDtarget <- match(qian_sankey$TARGET, qian_sankey.Nodes$name)-1

# sankeyNetwork(Links = qian_sankey, Nodes = qian_sankey.Nodes,
#               Source = 'IDsource', Target = 'IDtarget',
#               Value = 'interactions', NodeID = 'name', fontSize = 20)

g_counts <- ggplot(qian_sankey, aes(x = TARGET, y = interactions, fill = TARGET)) + 
  geom_bar(stat = 'identity') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  scale_fill_manual(values = wes_palette('GrandBudapest1')) +
  facet_wrap(~SOURCE, scales = 'free', nrow=2)
ggsave(filename = '~/Projects/HRD_TranscriptionalSignature/Figures/Figure7/Bassez_countBarplot.pdf',
       plot = g_counts, width = 6, height = 3.5)
