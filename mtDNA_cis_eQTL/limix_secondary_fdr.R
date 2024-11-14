#!/usr/bin/env Rscript

#' -----------------------------------------------------------------------------
#' Script Title: Secondary Multiple Testing for mtDNA cis eQTL
#' Description: This script performs secondary multiple testing on cis-eQTL data
#'              for various tissues. It calculates q-values and identifies
#'              significant results across multiple tissues.
#' Author: Xenofon Giannoulis
#' -----------------------------------------------------------------------------

library(qvalue)  # Load qvalue package for multiple testing correction
setwd("/workingdir") #set working directory 

# Function to read and process tissue data
process_tissue_data <- function(tissue, POST_TFB, bim) {
  r_t <- read.table(paste0(POST_TFB, tissue, ".top_qtl_results_all.txt"), header = TRUE, as.is = TRUE)
  r_t$tissue <- noquote(tissue)
  r_t$tissueqvalue <- qvalue(r_t$empirical_feature_p_value, lambda = 0)$qvalue
  r_t$A1 <- bim$V5[match(r_t$snp_id, bim$V2)]
  r_t$A0 <- bim$V6[match(r_t$snp_id, bim$V2)]
  r_t$beta <- 0 - r_t$beta  # Reverse the beta values, because tested allele in limix is major allele, different from all other methods, manually correct here 
  return(r_t)
}

# Function to process significant results
process_significant_results <- function(sig_TFB, POST_SFO) {
  SFO <- NULL
  for (t in unique(sig_TFB$tissue)) {
    r_t <- read.table(paste0(POST_SFO, t, ".qtl_results_all.txt"), header = TRUE, as.is = TRUE)
    r_t$tissue <- noquote(t)
    
    # Apply threshold from permuted empirical p-values
    threshold <- max(sig_TFB$empirical_feature_p_value[sig_TFB$tissue == t])
    pthreshold <- sig_TFB$p_value[sig_TFB$empirical_feature_p_value == threshold]
    r_t$adjusted_threshold <- as.numeric(pthreshold)
    
    # Classify results as significant or not
    r_t$sig <- ifelse(r_t$p_value <= pthreshold, "Yes", "No")
    SFO <- rbind(SFO, r_t)
  }
  return(SFO)
}

# Main execution
tissues_tested <- c("Artery_Aorta", "Brain_Spinal_cord_cervical_c-1", "Minor_Salivary_Gland", 
                     "Vagina", "Brain_Putamen_basal_ganglia", "Brain_Cerebellar_Hemisphere", 
                     "Brain_Hypothalamus", "Brain_Frontal_Cortex_BA9", "Brain_Caudate_basal_ganglia", 
                     "Artery_Coronary", "Spleen", "Brain_Nucleus_accumbens_basal_ganglia", 
                     "Prostate", "Brain_Cortex", "Brain_Cerebellum", "Pituitary", "Pancreas", 
                     "Stomach", "Colon_Sigmoid", "Testis", "Esophagus_Gastroesophageal_Junction", 
                     "Colon_Transverse", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", 
                     "Breast_Mammary_Tissue", "Esophagus_Muscularis", "Adipose_Visceral_Omentum", 
                     "Cells_Cultured_fibroblasts", "Esophagus_Mucosa", 
                     "Skin_Not_Sun_Exposed_Suprapubic", "Lung", "Nerve_Tibial", 
                     "Artery_Tibial", "Adipose_Subcutaneous", "Thyroid", 
                     "Skin_Sun_Exposed_Lower_leg", "Whole_Blood", "Muscle_Skeletal")

cases <- c("3_logTPM_scale")
SAVE <- getwd()
TFB <- NULL
SFO <- NULL

## read in genotype bim file for matching 
bim <- read.table("allv8.eur.maf01.bim", as.is = TRUE)

for (case in cases) {
  print(case)
  POST_TFB <- paste0(SAVE, case, '/TFB/')
  POST_SFO <- paste0(SAVE, case, '/SFO/')

  for (t in tissues_tested) {
    TFB <- rbind(TFB, process_tissue_data(t, POST_TFB, bim))
    print(t)
  }
}

TFB$qvalue <- qvalue(TFB$empirical_feature_p_value)$qvalue
TFB$tissues_tested <- length(tissues_tested)
TFB$round <- "secondary"

sig_TFB <- TFB[TFB$tissueqvalue < 0.05, ]
SFO <- process_significant_results(sig_TFB, POST_SFO)

# Further processing for SFO
SFO$A1 <- bim$V5[match(SFO$snp_id, bim$V2)]
SFO$A0 <- bim$V6[match(SFO$snp_id, bim$V2)]
SFO$beta <- 0 - SFO$beta 
SFO$round <- "secondary"
SFO$tissues_tested <- length(tissues_tested)

# Save results
write.table(sig_TFB, 'GTEX_cis_eQTL_TFB_sig_results_secondary.txt', 
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(SFO[SFO$sig == "Yes", ], 'GTEX_cis_eQTL_SFO_sig_results_secondary.txt', 
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(TFB, 'GTEX_TFB_secondary_all.txt', 
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(SFO, 'GTEX_SFO_secondary_all.txt', 
            row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
