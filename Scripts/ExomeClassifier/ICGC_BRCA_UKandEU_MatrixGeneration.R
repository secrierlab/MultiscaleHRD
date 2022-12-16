#####
## Processing of ICGC WGS breast cancer samples to infer SBS and indel mutation type contributions
####

# Load libraries
library(maftools)
library(sigminer)
library(BSgenome.Hsapiens.UCSC.hg19)

## Load UK and EU BRCA data (US is exome, FR contains no indels)
# Extract WGS strategies, remove duplicated mutations

brca_uk <- read.table('~/Downloads/simple_somatic_mutation.open.BRCA-UK.tsv',
                      sep = '\t', header = TRUE)
brca_uk_wgs <- brca_uk[brca_uk$sequencing_strategy == 'WGS', ]
brca_uk_wgs <- brca_uk_wgs[!duplicated(brca_uk_wgs$icgc_mutation_id), ]

brca_eu <- read.table('~/Downloads/simple_somatic_mutation.open.BRCA-EU.tsv',
                      sep = '\t', header = TRUE)
brca_eu_wgs <- brca_eu[brca_eu$sequencing_strategy == 'WGS', ]
brca_eu_wgs <- brca_eu_wgs[!duplicated(brca_eu_wgs$icgc_mutation_id), ]

# Collate ICGC data and organise into MAF-readable format

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

brca_maf <- read.maf(maf = brca_wgs_input,
                     vc_nonSyn = names(table(brca_wgs_input$Variant_Classification)))

# Run mt_tally() from sigminer package to collate mutation type contributions

mt_tally.brca_wgs <- sig_tally(
  object = brca_maf,
  ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  mode = 'ALL',
  useSyn = TRUE
)

save(mt_tally.brca_wgs, 
     file = '~/Documents/GitHub/HRD_classification/Data/ICGC_BRCA/BRCA_UKEU_mt_tally.Rdata')

