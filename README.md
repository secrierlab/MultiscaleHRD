# Multiscale HRD classification

### Author: Daniel H Jacobson, UCL

This repository contains the workflow for classification of homologous recombination deficient breast cancers via mutational signature classification and transcriptomic signature approaches. Scripts are separated into these two groups, and are presented here in the order in which they must be executed.

![alt text](TCGA_HRDclassificationHeatmap.jpg)

# Scripts

## Exome Classifier

- **ICGC_BRCA_UKandEU_MatrixGeneration.R:** Collates mutation data from the BRCA-EU and BRCA-UK projects and uses the sigminer R package to collate the mutational spectra of each sample according to the SBS96 and ID83 signatures.
- **ICGC_deconstructSigs_genome2exome.R:** Calculates contributions of breast cancer-associated SBS and ID mutational signatures to each of the ICGC samples, establishing the signature profiles of each sample, alongside correction for genome-to-exome normalisation.
- **ICGC_PhenotypeClusterDevelopment_normalised.R:** Applies finite mixture modelling to cluster the ICGC samples depending on their signature profiles, and enables name assignment to each cluster generated based on their most prevalent signature contributions.
- **LikelihoodCluster_Generation.R:** Collates the signature profiles for all samples within a cluster and generates a mean mutational spectrum representative of them, therefore generating the probability distributions for each cluster.
- **TCGA_HRDclassification.R:** Using the prior probabilities and likelihoods generated using ICGC, we collate the mutational profiles of 986 exome sequenced breast cancers from TCGA and calculate the posterior probabilities of assignment to each of the ICGC-generated clusters. The sum of probabilities of assignment to the HRD-associated clusters equals to the probability of HRD.
- **TCGA_HRDhallmarks.R:** Comparison of HRD and HR-proficiency assignments in TCGA to HRD-associated features (Myriad HRD score, CX3 copy number signature contribution, POLQ expression, proliferation capacity)

## TranscriptionalSignature

- **TCGA_BRCA.RNAseq_prep.R:** Pre-processing of TCGA-BRCA expression data including removal of lowly expressed genes, cancer cell expression deconvolution using BayesPrism, and separation into training (~2/3) and testing (~1/3) cohorts. This includes both HRD/HR-proficiency assignment according to the exome classifier, as well as BRCA1/BRCA2/HRD_BRCA+/HR-proficiency classifications which are used for signature development.
- **MultinomialElasticNet_alpha0.25_1to100.R:** Performs 100 iterations of 10-fold cross-validated multinomial elastic net regression. On each iteration, the gene parameter coefficients are saved. To run all 1000 iterations, an addition nine scripts were run in parallel with the seed set to the iteration value. Additionally, these analyses were repeated with alpha = 0.5 to generate an alternative signature.
- **CentroidModelFormation.R:** Collates the coefficients from the 1000 iterations of elastic net regression and extracts the 228 genes which appear as non-zero in all of them. The median expression of each gene is calculated across te HRD/HR-proficiency and HRD/BRCA-defect groups to generate the templates.
- **TCGA_testScoring.R:** Correlate the TCGA-BRCA testing cohort against each of the templates, saving the Pearson's correlation coefficient which represents the associated 'score'. The 'HRD score' is calculated by subtracting the correlation with the HR-proficiency template against the correlation with the HRD template.
- **GSEA_pathfindR.R:** Runs gene set enrichment analysis using the pathfindR tool. To enable this, for each gene an ANOVA is run of correlation against the HRD/BRCA-defect group, with the significance saved and adjusted. 
- **CCLE_jointComparisons:** Analysis of associations between HRD scores and PARP inhibitor sensitivity in breast cancer cell lines obtained from the Cancer Cell Line Encyclopaedia, and comparison against alternative HRD signatures.
- **ISPY2_HRDscoring.R:** Analysis of HRD scores across breast cancer patients treated with olaparib and durvalumab as part of the I-SPY2 trial, and comparison with the PARPi7 score
- **Chung2017_analysis.R:** Separates the bulk and single cell expression profiles from the Chung 2017 single cell breast cancer atlas. The tumour cells are extracted from the single cell data. HRD scores are calculated for each sample and tumour cell, and the sample-wide HRD scores are compared with the mean scores generated in the individual tumour cells.
- **Qian2020_preprocessing.R:** Preprocessing of the Qian 2020 breast cancer cohort, including removal of cells with high cell stress and unreasonable gene counts, and normalisation of expression scores.
- **Qian2020_analysis.R:** Analyses of distributions of HRD scores across the Qian 2020 cohort, and file preparation for CellphoneDB analysis
- **Qian2020_HRDprofiling.R** Analyses of distributions of HRD scores across cancer cells and the tumour microenvironment across the Qian 2020 cohort and UMAP plotting.


# Copyright

This code is free and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the GNU General Public License for more details.
