# Extract SBS and ID signatures in the ICGC-BRCA cohort

setwd('~/Data/ICGC/')

library(tidyr)

# Load SBS data
sbs_pcawg <- read.csv('PCAWG_sigProfiler_SBS_signatures_in_samples.csv', row.names = 2)
sbs_nopcawg <- read.csv('nonPCAWG_WGS_sigProfiler_SBS_signatures_in_samples_2018_04_13.csv', row.names = 2)

sbs <- rbind(sbs_pcawg, sbs_nopcawg)
sbs_brca <- sbs[grepl(pattern = 'Breast', sbs$Cancer.Types), -c(1,2)]

sbs_brca_props <- data.frame(
  Sigs = colnames(sbs_brca),
  Prop_present = apply(sbs_brca, 2, function(x) sum(x > 0)/length(x))
)

# Load ID data
id <- read.csv('PCAWG_SigProfiler_ID_signatures_in_samples.csv', row.names = 2)
id_brca <- id[grepl(pattern = 'Breast', id$Cancer.Types), -c(1,2)]

id_brca_props <- data.frame(
  Sigs = colnames(id_brca),
  Prop_present = apply(id_brca, 2, function(x) sum(x > 0)/length(x))
)

# Collate and save
brca_sigs <- rbind(sbs_brca_props, id_brca_props)
write.table(brca_sigs, file = 'ICGC_BRCA_sigProfilerCont.txt',
            row.names = FALSE, col.names = TRUE)
