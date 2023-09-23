setwd('~/Projects/HRD_TranscriptionalSignature/CellphoneDB/Results/')

library(dplyr)
library(tidyr)
library(ggpubr)
library(wesanderson)

cellType.of.interest <- 'Myeloid_cell'

# Process Qian
qian <- read.delim('Qian2020/significant_means.txt')
qian <- qian[,c('interacting_pair',
                paste(cellType.of.interest,'CancerHR.proficient',sep='.'),
                paste(cellType.of.interest,'CancerHRD',sep='.'))]
qian$status <- NA
qian$status[!is.na(qian[,2])] <- 'ONLY HR-proficient'
qian$status[!is.na(qian[,3])] <- 'ONLY HRD'
qian$status[!is.na(qian[,2]) &
              !is.na(qian[,3])] <- 'Both'
qian <- qian[!is.na(qian$status), ]
qian$Dataset <- 'Qian2020'

# Process Bassez
bassez <- read.delim('Bassez2021/significant_means.txt')
bassez <- bassez[,c('interacting_pair',
                    paste(cellType.of.interest,'Cancer_cellHR.proficient',sep='.'),
                    paste(cellType.of.interest,'Cancer_cellHRD',sep='.'))]
bassez$status <- NA
bassez$status[!is.na(bassez[,2])] <- 'ONLY HR-proficient'
bassez$status[!is.na(bassez[,3])] <- 'ONLY HRD'
bassez$status[!is.na(bassez[,2]) &
              !is.na(bassez[,3])] <- 'Both'
bassez <- bassez[!is.na(bassez$status), ]
bassez$Dataset <- 'Bassez2021'

names(bassez) <- names(qian)

# Merge
df <- rbind(qian,bassez)
names(df) <- c('interacting_pair','int_HRproficient','int_HRD','status','Dataset')
df 
df <- df %>%
  pivot_longer(cols = c(int_HRproficient, int_HRD),
               names_to = 'CancerStatus', values_to = 'Interaction')
df$group <- paste(df$Dataset, df$CancerStatus, sep = '_')
df$group <- factor(df$group,
                   levels = c('Qian2020_int_HRD','Bassez2021_int_HRD',
                              'Qian2020_int_HRproficient','Bassez2021_int_HRproficient'))

g_ints <- ggballoonplot(df, x = 'interacting_pair', y = 'group',
              size = 'Interaction', fill = 'status') +
  scale_fill_manual(values = c('grey90',wes_palette('GrandBudapest1')[1:2])) +
  ggtitle(paste0('Cell type: ',cellType.of.interest)) +
  geom_hline(yintercept = 2.5, color = 'red', linetype = 'dashed') +
  theme(legend.position = 'top')
ggsave(filename = paste0('~/Projects/HRD_TranscriptionalSignature/Figures/Supp_CPDB_Interactions_',cellType.of.interest,'.pdf'),
       plot = g_ints, height = 4, width = 10)
