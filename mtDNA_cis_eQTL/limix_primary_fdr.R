#!/usr/bin/env Rscript

#' -----------------------------------------------------------------------------
#' Script Title: cis-eQTL Analysis
#' Description: This script processes cis-eQTL results from various tissues,
#' applies multiple testing correction, and saves significant results.
#' Author: Xenofon Giannoulis
#' -----------------------------------------------------------------------------

# Load required libraries
library(qvalue)

# Define global parameters
cases <- c("3_logTPM_scale")
SAVE <- "/path/to/save/directory/"
tissues <- c("Artery_Aorta", "Brain_Substantia_nigra", "Uterus", 
             "Cells_EBV-transformed_lymphocytes", "Brain_Spinal_cord_cervical_c-1", 
             "Minor_Salivary_Gland", "Brain_Amygdala", "Vagina", 
             "Brain_Anterior_cingulate_cortex_BA24", "Ovary", 
             "Small_Intestine_Terminal_Ileum", "Brain_Hippocampus", 
             "Brain_Putamen_basal_ganglia", "Brain_Cerebellar_Hemisphere", 
             "Brain_Hypothalamus", "Brain_Frontal_Cortex_BA9", 
             "Brain_Caudate_basal_ganglia", "Artery_Coronary", 
             "Liver", "Spleen", "Brain_Nucleus_accumbens_basal_ganglia", 
             "Prostate", "Brain_Cortex", "Brain_Cerebellum", 
             "Adrenal_Gland", "Pituitary", "Pancreas", "Stomach", 
             "Colon_Sigmoid", "Testis", "Esophagus_Gastroesophageal_Junction", 
             "Colon_Transverse", "Heart_Atrial_Appendage", 
             "Heart_Left_Ventricle", "Breast_Mammary_Tissue",  
             "Esophagus_Muscularis", "Adipose_Visceral_Omentum", 
             "Cells_Cultured_fibroblasts", "Esophagus_Mucosa", 
             "Skin_Not_Sun_Exposed_Suprapubic", "Lung", 
             "Nerve_Tibial", "Artery_Tibial", "Adipose_Subcutaneous", 
             "Thyroid", "Skin_Sun_Exposed_Lower_leg", "Whole_Blood", 
             "Muscle_Skeletal")

bim <- read.table("/path/to/bim/file.bim", as.is = TRUE)

# Function to process TFB results
process_tfb_results <- function(case, tissues, SAVE) {
    TFB <- NULL
    POST_TFB <- paste0(SAVE, case, '/TFB/')

    for (t in tissues) {
        r_t <- read.table(paste0(POST_TFB, t, ".top_qtl_results_all.txt"), header = TRUE, as.is = TRUE)
        r_t$tissue <- noquote(t)
        r_t$tissueqvalue <- qvalue(r_t$empirical_feature_p_value, lambda = 0)$qvalue
        TFB <- rbind(TFB, r_t)
        print(t)
    }

    TFB$qvalue <- qvalue(TFB$empirical_feature_p_value)$qvalue
    TFB$A1 <- bim$V5[match(TFB$snp_id, bim$V2)]
    TFB$A0 <- bim$V6[match(TFB$snp_id, bim$V2)]
    TFB$beta <- -TFB$beta  # reverse beta for consistency
    TFB$tissues_tested <- length(tissues)
    TFB$round <- "primary"

    return(TFB)
}

# Process TFB results
for (case in cases) {
    print(case)
    TFB <- process_tfb_results(case, tissues, SAVE)
}

# Function to process SFO results
process_sfo_results <- function(sig_TFB, tissues, SAVE) {
    SFO <- NULL
    POST_SFO <- paste0(SAVE, case, '/SFO/')

    for (t in unique(sig_TFB$tissue)) {
        r_t <- read.table(paste0(POST_SFO, t, ".qtl_results_all.txt"), header = TRUE, as.is = TRUE)
        r_t$tissue <- noquote(t)
        
        threshold <- max(sig_TFB$empirical_feature_p_value[sig_TFB$tissue == t])
        pthreshold <- sig_TFB$p_value[sig_TFB$empirical_feature_p_value == threshold]
        r_t$adjusted_threshold <- as.numeric(pthreshold)
        
        r_t$sig <- ifelse(r_t$p_value <= pthreshold, "Yes", "No")
        SFO <- rbind(SFO, r_t)
        print(t)
    }

    return(SFO)
}

# Generate significant results for SFO
sig_TFB <- TFB[TFB$tissueqvalue < 0.05, ]
SFO <- process_sfo_results(sig_TFB, tissues, SAVE)

# Additional processing and output
SFO$key <- paste(SFO$tissue, SFO$feature_id, SFO$snp_id, sep = ".")
SFO$A1 <- bim$V5[match(SFO$snp_id, bim$V2)]
SFO$A0 <- bim$V6[match(SFO$snp_id, bim$V2)]
SFO$beta <- -SFO$beta
SFO$round <- "primary"

# Write results to files
write.table(sig_TFB, '/path/to/results/GTEX_cis_eQTL_TFB_sig_results_primary.txt', row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(SFO, '/path/to/results/GTEX_SFO_primary_all.txt', row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

# Setup for secondary analysis
out <- "/path/to/secondary/3_logTPM_scale/fvc/"
tissues_tested <- unique(sig_TFB$tissue)

for (t in tissues_tested) {
    r_save <- sig_TFB[sig_TFB$tissue == t, c("snp_id", "feature_id")]
    names(r_save)[2] <- "feature"
    write.table(r_save, paste0(out, t, "_fvc_input.txt"), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    
    # Create directories for round 2
    dir.create(paste0("/path/to/secondary/3_logTPM_scale/", t))
    print(t)
}
