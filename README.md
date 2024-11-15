# GTEx_Mitonuclear_eQTL
Code used in analyses in the paper: Giannoulis X., et al. Interplay between mitochondrial and nuclear DNA in gene expression regulation, BioRxiv 2024

# pre-processing and quality control 

Scripts for identifying EUR samples in GTEx using genotype and population data from [1000G Phase 3](https://www.internationalgenome.org/category/phase-3/) are in the ```Preprocessing``` directory
Scripts for obtaining per-tissue [PEER factors](https://www.nature.com/articles/nprot.2011.457) that account for batch effects and other unknown confoundings present in  ```PEER``` directory

# Genetic relatedness matrix (GRM)
For all eQTLs analyses the same genetic relatedness matrix (GRM) is used; this GRM is calculated using the following steps from 5,523,421 common SNPs (MAF > 5%, missingness < 0.1, P value for HWE > 10-6) in 684 European samples in GTEx who we use in all analyses we perform in this paper

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
All scripts for primary and secondary analyses are shown in the ```mtDNA-cis-eQTL``` directory

# mtDNA trans-eQTLs and nucDNA trans-eQTLs 

All trans-eQTL analyses in the paper are performed per tissue using [LDAK v5](https://dougspeed.com/).
For mtDNA trans-eQTLs between mtDNA SNPs and nucDNA gene expression, genotype $bfile contains mtDNA SNPs, and gene expression $featurefile contains nucDNA gene expression values. For nucDNA trans-eQTLs between nucDNA SNPs and mtDNA gene expression, genotype $bfile contains nucDNA SNPs, and gene expression $featurefile contains mtDNA gene expression values. The same GRM estimated using nucDNA SNPs (as above) is used in both analysis as $kin. $covfile contains sex, age and PEER factors specific for each tissue. 

```
ldak --bfile $bfile --pheno $featurefile --covar $covfile --linear $outdir/$z --grm $kin --max-iter 50000
```

# Celltype interaction QTL analysis (ct-iQTLs) 

# Enrichment analyses 

# Colocalization analyses 

# Mendelian randomization
