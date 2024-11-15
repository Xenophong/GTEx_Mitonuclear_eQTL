# GTEx Mitonuclear eQTL
Code used in analyses in the paper: Giannoulis X., et al. Interplay between mitochondrial and nuclear DNA in gene expression regulation, BioRxiv 2024

# pre-processing and quality control 

Scripts for identifying EUR samples in [GTEx v8](https://gtexportal.org/home/) using genotype and population data from [1000G Phase 3](https://www.internationalgenome.org/category/phase-3/) are in the ```Preprocessing``` directory
Scripts for obtaining per-tissue [PEER factors](https://www.nature.com/articles/nprot.2011.457) that account for batch effects and other unknown confoundings present in  ```PEER``` directory

# Genetic relatedness matrix (GRM)
For all eQTLs analyses the same genetic relatedness matrix (GRM) is used; this GRM is calculated using the following steps from 5,523,421 common SNPs (MAF > 5%, missingness < 0.1, P value for HWE > 10-6) in 684 European samples in GTEx who we use in all analyses we perform in this paper, using [LDAK v5](https://dougspeed.com/).

1. Calculation of SNP weights
```
ldak --bfile $bfile --cut-weights $wdir/sections
sections=$(cat $wdir/sections/section.number)
for section in $(seq 1 $sections)
do echo $section
ldak --bfile $bfile --calc-weights $wdir/sections --section $section
done
ldak --bfile $bfile --join-weights $wdir/sections
```

2. Calculation of GRM
```
ldak --bfile $bfile --cut-kins $wdir/partitions --by-chr YES
partitions=$(cat $wdir/partitions/partition.number)
for partition in $(seq 1 $partitions)
do echo $partition
ldak --bfile $bfile --calc-kins $wdir/partitions --partition $partition --weights $wdir/ldak/sections/weights.all --power -1
done
ldak --bfile $bfile --join-kins $wdir/partitions --kinship-raw YES
```

3. Perform PCA on GRM (for analyses that are run in R)
```
ldak --bfile $bfile --grm $wdir/partitions/kinships.all --pca $wdir/partitions/kinships.all --axes 20
```

# mtDNA cis-eQTLs 

All trans-eQTL analyses in the paper are performed per tissue using [LIMIX eQTL pipeline](https://github.com/single-cell-genetics/limix_qtl)
All scripts for primary and secondary analyses are shown in the ```mtDNA_cis_eQTL``` directory

# mtDNA trans-eQTLs and nucDNA trans-eQTLs 

All trans-eQTL analyses in the paper are performed per tissue using [LDAK v5](https://dougspeed.com/).
For mtDNA trans-eQTLs between mtDNA SNPs and nucDNA gene expression, genotype $bfile contains mtDNA SNPs, and gene expression $featurefile contains nucDNA gene expression values. For nucDNA trans-eQTLs between nucDNA SNPs and mtDNA gene expression, genotype ```$bfile``` contains nucDNA SNPs, and gene expression ```$featurefile``` contains mtDNA gene expression values. The same GRM estimated using nucDNA SNPs (as above) is used in both analysis as ```$kin```. ```$covfile``` contains sex, age and PEER factors specific for each tissue. 

```
ldak --bfile $bfile --pheno $featurefile --covar $covfile --linear $outdir/$z --grm $kin --max-iter 50000
```

# Celltype interaction QTL analysis (ct-iQTLs) 

We performed cell-type interaction QTL analysis (ct-iQTL) on all significant mtDNA cis-eQTL, mtDNA trans-eQTL and nucDNA trans-eQTL. For these analyses we use [xCell scores](https://github.com/dviraran/xCell) obtained per sample per tissue, for 7 cell types in total, as previously described in [Kim-Hellmuth et al Science 2020](https://www.science.org/doi/10.1126/science.aaz8528) for GTEx samples. 

Following recommendation in [Kim-Hellmuth et al Science 2020](https://www.science.org/doi/10.1126/science.aaz8528), we only performed cell-type interaction QTL analysis on cell types in each tissue where the median xCell score across all samples is > 0.1. 

All scripts performing the ct-iQTL analysis are shown in the ```ct_iQTL``` directory. 

# Colocalization analyses 

We performed statistical colocalization between nucDNA trans-eQTLs on mtDNA encoded genes and nucDNA cis-eQTLs on nucDNA-eGenes (which we hypothesize may mediate the nucDNA trans-eQTLs on mtDNA encoded genes), using the [coloc R package](https://cran.r-project.org/web/packages/coloc/index.html). All script for formatting the data as well as performing the analyses are shown in the ```Coloc``` directory.

# Enrichment analyses 

We performed enrichment analysis for nucDNA genes identified to be associated with mtDNA SNPs in the mtDNA trans-eQTL analysis. We performed two different enrichment analysis, one for annotated pathways, and one for previously identified disease associations. 

The former is performed using the [pathfindR R package](https://github.com/egeulgen/pathfindR), which identifies enrichment in specified databases. We used the following databases: 

```
genesets=c("KEGG","Reactome","BioCarta","GO-All")
for (geneset in genesets){
out=run_pathfindR(input,gene_sets=geneset,enrichment_threshold=0.05)
}
```

The latter is performed using the [disgenet2r R package](https://github.com/jinfar/disgenet2r), we show scripts for running and processing of disgenet2r results in the ```post_eQTL``` directory. 

# Mendelian randomization

To perform Mendelian Randomization (MR) between nucDNA trans-eQTLs on mtDNA encoded genes and GWAS on complex traits and diseases, we obtained summary statistics at all nucDNA trans-eQTLs on mtDNA encoded genes at complex traits and diseases submitted to the [GWAS Cataloge](https://www.ebi.ac.uk/gwas/) using the [LDlinkR R package](https://cran.r-project.org/web/packages/LDlinkR/index.html). 

We then performed analyses between nucDNA trans-eQTLs on mtDNA encoded genes and the GWAS on complex traits and diseases identified, using the [Mendelian Randomization R package](https://cran.r-project.org/web/packages/MendelianRandomization/index.html). 

All script for extracting the GWAS on complex traits and diseases data, as well as MR analyses, are shown in the ```post_eQTL``` directory

