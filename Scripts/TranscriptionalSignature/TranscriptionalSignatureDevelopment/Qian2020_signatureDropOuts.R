setwd('~/Data/scRNASeq/Qian2020/')

load('exprData_Qian2020.Rdata')

meta.qian <- read.csv('2103-Breastcancer_metadata.csv', header = TRUE)
meta.qian <- meta.qian[match(colnames(expr.data_qian2020), meta.qian$Cell), ]
expr.cancer <- expr.data_qian2020[,meta.qian$CellType == 'Cancer']

# Collate genes lists
genes.g0 <- read.csv('~/Data/QuiescenceBiomarkers.csv')
genes.g0 <- intersect(genes.g0$Genes, rownames(expr.cancer))

load('~/Projects/HRD_TranscriptionalSignature/Results/centroids_complete_Run6_HRD0.79_zNorm.Rdata')
genes.hrd <- rownames(signature.centroid.list$ElasticNet_alpha0.25)

library(readxl)
genes.ddr_table <- readxl::read_excel('../../../../Downloads/Pearl_DDRgenes.xlsx', skip=1)
genes.ddr <- genes.ddr_table$`Gene ID`
genes.ddr <- genes.ddr[!is.na(genes.ddr)]
genes.ddr <- genes.ddr[!duplicated(genes.ddr)]
genes.ddr <- intersect(genes.ddr, rownames(expr.cancer))

# Summarise
summary(apply(expr.cancer[genes.g0,], 1, function(x) mean(x>0)))
hist(apply(expr.cancer[genes.g0,], 1, function(x) mean(x>0)), breaks = 50, 
     xlim = c(0,1), main = 'G0 genes expression')

summary(apply(expr.cancer[genes.ddr, ], 1, function(x) mean(x>0)))
hist(apply(expr.cancer[genes.ddr, ], 1, function(x) mean(x>0)), breaks = 50, 
     xlim = c(0,1), main = 'DDR genes expression')

summary(apply(expr.cancer[genes.hrd,], 1, function(x) mean(x>0)))
hist(apply(expr.cancer[genes.hrd,], 1, function(x) mean(x>0)), breaks = 50, 
     xlim = c(0,1), main = 'HRD genes expression')

summary(apply(expr.cancer[genes.hrd,], 2, function(x) sum(x>0)))

# 
# summary(apply(expr.cancer, 1, function(x) mean(x>0)))
# hist(apply(expr.cancer, 1, function(x) mean(x>0)), breaks = 50, 
#      xlim = c(0,1), main = 'All genes expression')

