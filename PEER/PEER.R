#!/usr/bin/env Rscript
#' -----------------------------------------------------------------------------
#' Run PEER on nuclear DNA gene expression data in 48 tissues in GTEx 
#' @author Xenofon Giannoulis <xenofon.giannoulis@helmholtz-muenchen.de>
#' @date Sat Feb  5 00:02:53 CET 2022
#' -----------------------------------------------------------------------------
'%!in%' <- function(x,y)!('%in%'(x,y))
library(peer); library(dplyr); library(data.table); library(ggplot2)

run_peer <- function(tissue, peer_in, ge_outdir) {
  print(paste0('% processing case % .... ', tissue))
  
  r_gene_expr <- paste0(peer_in, tissue, '.forpeer.tab')
  r_samples <- paste0(peer_in, tissue, '.forpeer.sampleids')
  r_genes <- paste0(peer_in, tissue, '.forpeer.geneids')
  dataanno <- paste0(peer_in, tissue, '.forpeer.anno.txt')
  r_peer_out <- paste0(ge_outdir, tissue, "/")
  
  peer_expr <- fread(r_gene_expr, sep='\t', header=F)
  idnames <- read.delim(r_samples, header=F)$V1
  mygenes <- read.delim(r_genes, header=F)$V1
  
  # Determine the number of factors based on sample size
  num_factors <- ifelse(dim(peer_expr)[1] < 150, 15,
                   ifelse(dim(peer_expr)[1] < 250, 30,
                   ifelse(dim(peer_expr)[1] < 350, 45, 60)))
  
  Nmax_iterations <- 2000
  model <- PEER()
  print(paste0("Starting unsupervised run with ", num_factors, " PEER factors for ", tissue, 
               " ... and maximum limit of iterations ", Nmax_iterations))

  # Set data and parameters
  PEER_setNk(model, num_factors)
  PEER_setPhenoMean(model, as.matrix(peer_expr))
  PEER_setPriorAlpha(model, 0.001, 0.1)
  PEER_setPriorEps(model, 0.1, 10)
  PEER_setNmax_iterations(model, Nmax_iterations)
  
  # Perform inference
  PEER_update(model)
  
  # Investigating results
  X <- PEER_getX(model) # factors
  W <- PEER_getW(model) # weights
  Alpha <- PEER_getAlpha(model) # ARD parameters
  Yc <- PEER_getResiduals(model) # get corrected dataset: Residuals
  
  my_peer_names <- paste0("PEER", seq(1:num_factors))
  rownames(Alpha) <- my_peer_names
  
  # Writing Results
  residuals <- scale(data.frame(Yc))
  colnames(residuals) <- mygenes
  rownames(residuals) <- idnames
  
  Alpha <- as.data.frame(Alpha)
  names(Alpha)[1] <- 'Alpha_value'
  
  factors <- data.frame(X)
  names(factors) <- my_peer_names
  factors$sample_id <- idnames    
  factors <- factors %>% select(sample_id, everything()) # limix ready only factors

  weights <- data.frame(W)
  colnames(weights) <- my_peer_names
  weights$feature_id <- mygenes
  weights <- weights %>% select(feature_id, everything())
  
  # Summary PEER plots
  peer_plots <- "/no-backup/xenofon/Emperor/22_GTEX/master/peer/plots/"
  tiff(filename=paste0(peer_plots, tissue, "_Factor_Relevance.PEER", num_factors, ".tiff"), 
       bg="white", units="in", width=12, height=12, res=80, type = c("cairo"))
  plot(1.0 / Alpha$Alpha_value, xlab="Factors", ylab="Factor relevance", main=paste0(tissue))
  dev.off()
  
  tiff(filename=paste0(peer_plots, tissue, "PlotModel.PEER", num_factors, ".tiff"), 
       bg="white", units="in", width=12, height=12, res=80, type = c("cairo"))
  PEER_plotModel(model)
  dev.off()
  
  r_corr <- subset(factors, select=-c(sample_id))
  r_cor_mat <- cor(as.matrix(r_corr))
  
  tiff(filename=paste0(peer_plots, tissue, "Correlation.PEER", num_factors, ".tiff"), 
       bg="white", units="in", width=12, height=12, res=80, type = c("cairo"))
  heatmap(x = r_cor_mat, symm = TRUE, main = paste0(tissue, "_Correlation.PEER", num_factors))
  dev.off()
  
  # Saving results
  write.table(Alpha, paste0(r_peer_out, "PEER", num_factors, "_Alpha.txt"), sep="\t", quote=F, row.names=T, col.names=T)
  write.table(weights, paste0(r_peer_out, "PEER", num_factors, "_Weights.txt"), sep="\t", quote=F, row.names=F, col.names=T)
  write.table(residuals, paste0(r_peer_out, "PEER", num_factors, "_residuals.scaled.txt"), sep="\t", quote=F, row.names=T, col.names=T)
  write.table(r_cor_mat, paste0(r_peer_out,".PEER", num_factors, "_factors_corr.txt"), sep="\t", quote=F, row.names=T, col.names=T)
  
}

# Main loop for all tissues
tissues <- c("Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Spinal_cord_cervical_c-1", 
             "Brain_Substantia_nigra", "Cells_EBV-transformed_lymphocytes", "Minor_Salivary_Gland", 
             "Small_Intestine_Terminal_Ileum", "Brain_Hippocampus", "Brain_Putamen_basal_ganglia", 
             "Brain_Cerebellar_Hemisphere", "Brain_Hypothalamus", "Brain_Frontal_Cortex_BA9", 
             "Brain_Caudate_basal_ganglia", "Artery_Coronary", "Liver", "Spleen", 
             "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Cortex", "Brain_Cerebellum", 
             "Adrenal_Gland", "Pituitary", "Pancreas", "Stomach", "Colon_Sigmoid", 
             "Esophagus_Gastroesophageal_Junction", "Colon_Transverse", "Heart_Atrial_Appendage", 
             "Artery_Aorta", "Heart_Left_Ventricle", "Breast_Mammary_Tissue", "Esophagus_Muscularis", 
             "Adipose_Visceral_Omentum", "Cells_Cultured_fibroblasts", "Esophagus_Mucosa", 
             "Skin_Not_Sun_Exposed_Suprapubic", "Lung", "Nerve_Tibial", "Artery_Tibial", 
             "Adipose_Subcutaneous", "Thyroid", "Skin_Sun_Exposed_Lower_leg", "Whole_Blood", 
             "Muscle_Skeletal", "Testis", "Uterus", "Vagina", "Ovary", "Prostate")

peer_in <- '/path/to/working/directory/'
ge_outdir <- '/path/to/output/directory/'

for (tissue in tissues) {
  run_peer(tissue, peer_in, ge_outdir)
}
