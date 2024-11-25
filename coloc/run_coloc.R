library(coloc)
## read in nucDNA trans-eQTLs for mtDNA genes as r_gwas 
## read in nucDNA cis-eQTLs for nucDNA cis-eGenes as r_re_eGene

## format nucDNA trans-eQTLs for mtDNA genes 
r_gwas_coloc_d1 <- r_gwas[, .(position = Basepair, 
                                beta = Effect, 
                                varbeta = SD^2,    # varbeta is the square of SE
                                snp = Predictor, 
                                MAF = MAF, 
                                n_e_samples = unique(n_e_samples), 
                                p_value = Wald_P, 
                                sdY_gwas = sdY_gwas)]  # Assuming sdY is a constant
## format nucDNA cis-eQTLs for nucDNA cis-eGenes 
r_eqtl_coloc_d1 <- r_re_eGene[, .(position = snp_position, 
                                beta = beta, 
                                varbeta = varbeta,
                                snp = snp_id, 
                                MAF = maf, 
                                n_e_samples = unique(n_e_samples), 
                                p_value = p_value, 
                                sdY_eqtl = sdY_eqtl)]  # Assuming sdY is a constant

## get input for coloc
input <- merge(r_eqtl_coloc_d1, 
                   r_gwas_coloc_d1, by="snp", 
                   all=FALSE, suffixes=c("_eqtl","_gwas"))
    
## run coloc 
r_result <- coloc.abf(
  dataset1=list(snp=input$snp,
                    pvalues=as.numeric(input$p_value_gwas), 
                    type="quant", 
                    N=unique(input$n_e_samples_gwas),
                    beta=input$beta_gwas,
                    varbeta=input$varbeta_gwas,
                    sdY=unique(input$sdY_gwas)), 
  dataset2=list(snp=input$snp,
                    pvalues=as.numeric(input$p_value_eqtl), 
                    type="quant", 
                    N=unique(input$n_e_samples_eqtl),
                    beta=(input$beta_eqtl),
                    varbeta=input$varbeta_eqtl,
                    sdY=unique(input$sdY_eqtl)), 
  MAF=input$MAF_gwas)
  



