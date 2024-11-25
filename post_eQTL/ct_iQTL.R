#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(magrittr)

# Define the interaction eQTL analysis function
interaction_eqtl_analysis <- function(
  tissue, snp, gene_expression_data, cell_counts_data, genotype_data, covariates_data, peer_factors, pca_data
) {
  # Input:
  # - tissue: name of the tissue being analyzed
  # - snp: SNP identifier
  # - gene_expression_data: data frame with gene expression values
  # - cell_counts_data: data frame with normalized cell types
  # - genotype_data: data frame with genotype values
  # - covariates_data: data frame with additional covariates (age, sex)
  # - peer_factors: data frame with PEER factors
  # - pca_data: data frame with PCA components

  # Merge all data sources on a common identifier (e.g., SUBID)
  merged_data <- reduce(
    list(gene_expression_data, genotype_data, cell_counts_data, covariates_data, peer_factors, pca_data),
    function(x, y) full_join(x, y, by = "SUBID")
  )
  
  # Ensure numeric columns are properly formatted
  merged_data <- merged_data %>%
    mutate(across(where(is.character), as.numeric))
  
  # Define covariates
  covariate_columns <- names(merged_data)[which(names(merged_data) == "SEX"):which(names(merged_data) == "PC10")]
  
  # Linear model for main effect of SNP
  main_formula <- as.formula(paste("Gene_Expression ~", snp, "+", paste(covariate_columns, collapse = "+")))
  main_model <- lm(main_formula, data = merged_data)
  main_results <- summary(main_model)
  
  # Interaction testing with cell types
  interaction_results <- list()
  cell_types <- setdiff(names(cell_counts_data), "SUBID")  # Exclude identifier column
  
  for (cell_type in cell_types) {
    interaction_formula <- as.formula(
      paste("Gene_Expression ~", snp, "+", cell_type, "+", paste(covariate_columns, collapse = "+"), 
            "+", paste0(snp, ":", cell_type))
    )
    
    interaction_model <- lm(interaction_formula, data = merged_data)
    interaction_summary <- summary(interaction_model)
    
    # Extract interaction term
    interaction_term <- paste0(snp, ":", cell_type)
    if (interaction_term %in% rownames(interaction_summary$coefficients)) {
      interaction_results[[cell_type]] <- list(
        Estimate = interaction_summary$coefficients[interaction_term, "Estimate"],
        SE = interaction_summary$coefficients[interaction_term, "Std. Error"],
        PValue = interaction_summary$coefficients[interaction_term, "Pr(>|t|)"]
      )
    } else {
      interaction_results[[cell_type]] <- list(Estimate = NA, SE = NA, PValue = NA)
    }
  }
  
  # Return results
  return(list(
    SNP_PValue = snp_pvalue,
    Interaction_Results = interaction_results
  ))
}
