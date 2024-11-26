#!/usr/bin/env Rscript
#' -----------------------------------------------------------------------------
#' Script Title: Partial Correlation of Expression Residuals for nucDNA SNPs
#' Description: Partial correlation analysis between nuclear SNPs and genes (nucDNA cis + trans)
#' Analyzing colocalized nucSNP-nuGenes-mtGenes
#' -----------------------------------------------------------------------------

# Load libraries
library(data.table)
library(dplyr)
library(tidyverse)

# Function to perform partial correlation analysis for colocalized eGenes
partial_correlation_eQTL <- function(line_mapper, nucDNA_cis_eQTL_results, trans_eQTL_results, pca, peer_cov_run) {

  # Initialize an empty data frame to store results
  results <- data.frame()

  # Iterate over each colocalized SNP-gene pair
  for (r_key in unique(line_mapper$key2)) {
    print(paste("Processing SNP-Gene pair:", r_key))

    # Filter the current SNP-gene pair from line_mapper and cis-eQTL results
    r_line_mapper <- line_mapper %>% filter(key2 == r_key)
    r_lines_mapper_tfb <- nucDNA_cis_eQTL_results %>% filter(key2 == r_key)
    
    # Get the list of colocalized nuclear eGenes for the current SNP
    colocalized_nucgenes <- unique(r_lines_mapper_tfb$gene)
    
    # Iterate through each colocalized nuclear eGene
    for (nucgene in colocalized_nucgenes) {
      print(paste("Processing eGene:", nucgene))
      
      # Filter the gene expression data for the current nuclear eGene
      # Load the gene expression for the current tissue and gene
      gene_expr_data <- read.delim(paste0("/path/to/gene_expression_data/", nucgene, ".txt"), header = TRUE)

      # Merge with covariates (PCA, sex-age, PEER residuals)
      covariates <- merge(peer_cov_run, pca, by = "FID")
      gene_expr_data <- merge(gene_expr_data, covariates, by = "FID")
      
      # Fit the linear model for partial correlation
      model_formula <- paste0(nucgene, " ~ ", paste0(colnames(covariates), collapse = " + "))
      lm_model <- lm(model_formula, data = gene_expr_data)
      
      # Extract residuals from the linear model
      residuals <- lm_model$residuals

      # Store the results (SNP-gene pair, nuclear eGene, and partial correlation results)
      results <- rbind(results, data.frame(
        SNP = r_key,
        eGene = nucgene,
        partial_residuals = residuals
      ))

    }  # End loop for each colocalized eGene
    
  }  # End loop for each SNP-gene pair

  # Return the final results
  return(results)
}

# Example of how to call the function
line_mapper <- read.delim('path_to_line_mapper.txt')
nucDNA_cis_eQTL_results <- read.delim('path_to_nucDNA_cis_eQTL_results.txt')
trans_eQTL_results <- read.delim('path_to_trans_eQTL_results.txt')
pca <- read.table('path_to_pca_file.txt', header = TRUE)
peer_cov_run <- read.table('path_to_peer_covariates.txt', header = TRUE)

# Call the function and store results
partial_corr_results <- partial_correlation_eQTL(line_mapper, nucDNA_cis_eQTL_results, trans_eQTL_results, pca, peer_cov_run)

# Write the results to a file (optional)
write.table(partial_corr_results, 'partial_correlation_results.txt', row.names = FALSE, quote = FALSE, sep = "\t")
