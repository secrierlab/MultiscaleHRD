setwd('~/Data/ICGC/')
#####
## Processing of ICGC WGS breast cancer samples to infer SBS and indel mutation type contributions
####

# Load libraries
library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)

## Load UK and EU BRCA data (US is exome, FR contains no indels)
# Extract WGS strategies, remove duplicated mutations

brca_uk <- read.table('~/Data/ICGC/simple_somatic_mutation.open.BRCA-UK.tsv.gz',
                      sep = '\t', header = TRUE)
brca_uk_wgs <- brca_uk[brca_uk$sequencing_strategy == 'WGS', ]
brca_uk_wgs <- brca_uk_wgs[!duplicated(brca_uk_wgs$icgc_mutation_id), ]

brca_eu <- read.table('~/Data/ICGC/simple_somatic_mutation.open.BRCA-EU.tsv.gz',
                      sep = '\t', header = TRUE)
brca_eu_wgs <- brca_eu[brca_eu$sequencing_strategy == 'WGS', ]
brca_eu_wgs <- brca_eu_wgs[!duplicated(brca_eu_wgs$icgc_mutation_id), ]

# Collate ICGC data and extract mutations in intron or intergenic regions

brca_wgs <- rbind(brca_uk_wgs, brca_eu_wgs)
brca_wgs_input <- data.frame(
  Tumor_Sample_Barcode = brca_wgs$icgc_donor_id,
  Hugo_Symbol = NA,
  Chromosome = brca_wgs$chromosome,
  Start_position = brca_wgs$chromosome_start,
  End_position = brca_wgs$chromosome_end,
  Variant_Classification = brca_wgs$consequence_type,
  Variant_Type = sapply(brca_wgs$mutation_type,
                        function(x) ifelse(x == 'single base substitution', 'SNP',
                                           ifelse(x == 'insertion of <=200bp', 'INS',
                                                  ifelse(x == 'deletion of <=200bp', 'DEL',NA)))),
  Reference_Allele = brca_wgs$reference_genome_allele,
  Tumor_Seq_Allele2 = brca_wgs$mutated_to_allele
)

# brca_exome_input <- brca_wgs_input[!(brca_wgs_input$Variant_Classification %in%
#                                        c('intergenic_region', 'intron_variant')), ]
brca_exome_input <- brca_wgs_input[brca_wgs_input$Variant_Classification != 'intergenic_region',]

brca_wgs_maf <- read.maf(maf = brca_wgs_input,
                         vc_nonSyn = names(table(brca_wgs_input$Variant_Classification)))

brca_exome_maf <- read.maf(maf = brca_exome_input,
                           vc_nonSyn = names(table(brca_exome_input$Variant_Classification)))

# Calculate ID-83 counts for each MAF file
mt_tally.brca_wgs <- sig_tally(
  object = brca_wgs_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ID',
  useSyn = TRUE
)

mt_tally.brca_exome <- sig_tally(
  object = brca_exome_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ID',
  useSyn = TRUE
)

# Count total mutations of each ID-83 type
icgc_indelCounts <- data.frame(
  wgs = apply(mt_tally.brca_wgs$all_matrices$ID_83, 2, sum),
  exome = apply(mt_tally.brca_exome$all_matrices$ID_83, 2, sum)
)

write.table(icgc_indelCounts, file = 'ICGC_BRCA_indelCounts.txt')

icgc_indelCounts$ex_gen <- icgc_indelCounts$exome/icgc_indelCounts$wgs
hist(icgc_indelCounts$ex_gen)

library(ggplot2)
ggplot(icgc_indelCounts, aes(x = log10(wgs), y = log10(exome))) +
  geom_point() + geom_smooth(method = 'lm')


