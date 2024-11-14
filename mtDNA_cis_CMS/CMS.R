
#!/usr/bin/env Rscript

#' -----------------------------------------------------------------------------
#' Script Title: Perform GTEx cis-eQTL on updated 147 SNPs: CMS
#' Description: This script performs cis-eQTL analysis on mitochondrial genes across multiple tissues using CMS.
#' Author: Xenofon Giannoulis <xenofon.giannoulis@helmholtz-muenchen.de>
#' Date: Sat Feb 5, 2022
#' -----------------------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(foreign)
library(qvalue)

# Source utility functions
source("~/scripts/utilities.R")
source("/path/to/CMS/CMS_v1.0.R")
source("/path/to/scripts/rowR.R")

# 1. Prepare the results dataframe
initialize_results_df <- function() {
  data.frame(
    tissue = NA, gene = NA, snp = NA, beta.MA = NA, SD.MA = NA,
    pval.MA = NA, beta.MC = NA, SD.MC = NA, pval.MC = NA,
    Yrsq = NA, Ncov = NA, Lcov = NA, pMUL = NA,
    PcutoffBias = NA, Lcov_names = NA, NPEER = NA, KeyID = NA
  )
}

# 2. Load data and perform CMS analysis
perform_cms_analysis <- function(features, tissues, my_predictor, PCA_dir, mtgenedir, load) {
  res_df <- initialize_results_df()
  
  for (feature in features) {
    for (tissue in tissues) {
      cat("##############[ ITERATION ", feature, "-", tissue, " ]##############\n")
      
      # Load genotype data
      geno_file <- file.path("/path/to/genotypes/", paste0(tissue, '.v8.eur.maf01.recoded.raw'))
      mtgeno <- read.delim(geno_file, header = TRUE, sep = ' ') %>%
        select(-IID, -PAT, -MAT, -SEX, -PHENOTYPE) %>%
        mutate(across(everything(), ~ replace(., . == "2", 1)))
      
      snps <- colnames(mtgeno)[-1]
      
      for (snp in snps) {
        res_df <- bind_rows(res_df, run_cms_for_snp(feature, tissue, snp, mtgeno, PCA_dir, mtgenedir, load))
      }
    }
  }
  return(res_df)
}

# 3. Run CMS for a specific SNP
run_cms_for_snp <- function(feature, tissue, snp, mtgeno, PCA_dir, mtgenedir, load) {
  # Load expression data
  expression_file <- file.path(mtgenedir, paste0(tissue, "_mtgenes.txt"))
  expression <- read.table(expression_file, header = TRUE, as.is = TRUE)
  expr <- expression[expression$feature_id == feature, ]
  
  exp_run <- t(expr[-1])
  exp_run <- as.data.frame(exp_run)
  colnames(exp_run) <- expr$feature_id
  exp_run <- rownames_to_column(exp_run, "FID") %>%
    mutate(FID = gsub(".", "-", FID, perl = FALSE, fixed = TRUE))
  
  # Prepare genotype data
  geno_run <- mtgeno %>%
    select(FID = 1, snp)
  
  # Merge data
  sub_total_run <- geno_run %>%
    inner_join(exp_run, by = "FID") %>%
    inner_join(load_covariates(tissue, load), by = "FID") %>%
    inner_join(load_pca_data(PCA_dir), by = "FID") %>%
    inner_join(load_peer_factors(tissue, load), by = "FID")
  
  # Prepare for CMS
  total_run <- select(sub_total_run, -FID)
  DATAMAT_X <- as.matrix(total_run)
  
  # Run CMS
  test <- MC(
    DATAMAT = DATAMAT_X,
    myOUTCOME = feature,
    myPREDICTOR = snp,
    myFIXCOV = c("AGE", "SEX", paste0("PCA", 1:10)),
    listCOVARIATES = colnames(total_run)[15:ncol(total_run)],
    optionStand = 1,
    optionImpute = 1,
    minTotSample = 60,
    sigmaMax = 2,
    dweigh = 8,
    Tmul = 0.05,
    verbose = FALSE
  )
  
  # Construct result row
  Lcov_names <- ifelse(test[[9]] == "0" | is.na(test[[9]]), "ZERO", 
                       paste(colnames(total_run)[as.numeric(unlist(strsplit(test[[9]], ",")))], collapse = ","))
  
  data.frame(
    tissue = tissue, 
    gene = feature, 
    snp = snp, 
    beta.MA = test[[1]], 
    SD.MA = test[[2]], 
    pval.MA = test[[3]], 
    beta.MC = test[[4]], 
    SD.MC = test[[5]],
    pval.MC = test[[6]], 
    Yrsq = test[[7]], 
    Ncov = ncol(total_run), 
    Lcov = ncol(sub_total_run) - ncol(total_run), 
    pMUL = test[[8]], 
    PcutoffBias = test[[10]], 
    Lcov_names = Lcov_names, 
    NPEER = ncol(load_peer_factors(tissue, load)) - 1, 
    KeyID = snp
  )
}

# 4. Load covariates
load_covariates <- function(tissue, load) {
  cov_file <- file.path(load, "covariates/limix_in/", paste0(tissue, "_sexage.txt"))
  read.table(cov_file, header = TRUE, as.is = TRUE) %>%
    rename(FID = 1)
}

# 5. Load PCA data
load_pca_data <- function(PCA_dir) {
  pca_file <- file.path(PCA_dir, "kinship.all.top10.vect")
  read.table(pca_file, header = TRUE, as.is = TRUE) %>%
    select(-1) %>%
    rename(FID = 1) %>%
    setNames(paste0("PCA", 1:10))
}

# 6. Load PEER factors
load_peer_factors <- function(tissue, load) {
  peer_file <- file.path(load, "covariates/limix_in/", paste0(tissue, "_allcovs.txt"))
  read.table(peer_file, as.is = TRUE, header = TRUE) %>%
    rename(FID = 1) %>%
    select(-c(SEX, AGE))
}

# 7. Main execution
features <- c("ENSG00000198695", "ENSG00000198712", "ENSG00000228253", 
              "ENSG00000212907", "ENSG00000198938", "ENSG00000198899", 
              "ENSG00000198888", "ENSG00000198886", "ENSG00000198840", 
              "ENSG00000198804", "ENSG00000198786", "ENSG00000198763", 
              "ENSG00000198727")

tissues <- c("Adipose_Subcutaneous", "Artery_Aorta", "Brain_Substantia_nigra", "Uterus", 
             "Cells_EBV-transformed_lymphocytes", "Brain_Spinal_cord_cervical_c-1", 
             "Minor_Salivary_Gland", "Brain_Amygdala", "Vagina", 
             "Brain_Anterior_cingulate_cortex_BA24", "Ovary", "Small_Intestine_Terminal_Ileum", 
             "Brain_Hippocampus", "Brain_Putamen_basal_ganglia", "Brain_Cerebellar_Hemisphere", 
             "Brain_Hypothalamus", "Brain_Frontal_Cortex_BA9", 
             "Brain_Caudate_basal_ganglia", "Artery_Coronary", "Liver", "Spleen", 
             "Brain_Nucleus_accumbens_basal_ganglia", "Prostate", "Brain_Cortex", 
             "Brain_Cerebellum", "Adrenal_Gland", "Pituitary", "Pancreas", 
             "Stomach", "Colon_Sigmoid", "Testis", "Esophagus_Gastroesophageal_Junction", 
             "Colon_Transverse", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", 
             "Breast_Mammary_Tissue", "Esophagus_Muscularis", "Adipose_Visceral_Omentum", 
             "Cells_Cultured_fibroblasts", "Esophagus_Mucosa", "Skin_Not_Sun_Exposed_Suprapubic", 
             "Lung", "Nerve_Tibial", "Artery_Tibial", "Thyroid", "Skin_Sun_Exposed_Lower_leg", 
             "Whole_Blood", "Muscle_Skeletal")

PCA_dir <- "/path/to/PCA/"
mtgenedir <- "/path/to/mtgenes/"
load <- "/path/to/load/"

# Execute analysis
results <- perform_cms_analysis(features, tissues, my_predictor, PCA_dir, mtgenedir, load)

# Save results
output_file <- "/path/to/output/results_cis_eQTL.csv"
write.csv(results, file = output_file, row.names = FALSE)
