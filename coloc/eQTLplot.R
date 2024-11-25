library(eQTpLot)

output <- capture.output(eQTpLot(
      Genes.df = r_gene_df_updated, # mapping of nuc-genes
      GWAS.df = r_gw, #nucDNA-trans-eQTL
      eQTL.df = r_qt,  #mtDNA-cis-eQTL
      gene = genes_2_test_list,  # significant cis-eGene
      trait = unique(r_gw$PHE),
      getplot=FALSE,
      gbuild  = "hg38",
      tissue  = unique(r_qt$Tissue),
      range=2000000,
      saveplot=TRUE,
      GeneList = T
    )
    )
