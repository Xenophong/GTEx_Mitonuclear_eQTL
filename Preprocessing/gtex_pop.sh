## extract only those SNPs from 1000G in GTEX  
wdir="/GTExV8/genotypes"
bfile=$wdir/allchr.maf5hwe6missing10
chromosomes=$(seq 1 22)
for chr in $chromosomes 
do echo $chr 
plink --bfile $bfile --make-bed --extract $wdir/chr$chr.maf5hwe6missing10.kgp.extract --out $wdir/chr$chr.maf5hwe6missing10.kgp; \
mv $wdir/chr$chr.maf5hwe6missing10.kgp.rsids $wdir/chr$chr.maf5hwe6missing10.kgp.bim; \
plink --bfile $wdir/chr$chr.maf5hwe6missing10.kgp --indep-pairwise 1000 10 0.2 --out $wdir/chr$chr.maf5hwe6missing10.kgp.ld2; \
plink --bfile $wdir/chr$chr.maf5hwe6missing10.kgp --extract $wdir/chr$chr.maf5hwe6missing10.kgp.ld2.prune.in --make-bed --out $wdir/chr$chr.maf5hwe6missing10.kgp.ld2
done 

## extract those SNPs from 1000g in 1000G
wdir="GTExV8/genotypes"
rm $wdir/mergeld2-list
chromosomes=$(seq 2 22)
bfile=$wdir/chr1.maf5hwe6missing10.kgp.ld2
for chr in $chromosomes
do echo $chr 
echo -e "$wdir/chr$chr.maf5hwe6missing10.kgp.ld2.bed $wdir/chr$chr.maf5hwe6missing10.kgp.ld2.bim $wdir/chr$chr.maf5hwe6missing10.kgp.ld2.fam">>$wdir/mergeld2-list
done 
plink --bfile $bfile --merge-list $wdir/mergeld2-list --make-bed --out $wdir/allchr.maf5hwe6missing10.kgp.ld2

## remove those SNPs with discordant alleles 
wdir="/GTExV8/genotypes"
bfile="$wdir/allchr.maf5hwe6missing10.kgp.ld2"
plink --bfile $bfile --exclude $wdir/discordant.exclude --make-bed --out $bfile.nodis

## use LDAK for GRM and PCs
bfile="$wdir/allchr.maf5hwe6missing10.kgp.ld2.nodis"
ldak5.linux --bfile $bfile --calc-kins-direct $bfile --ignore-weights YES --power -1 
ldak5.linux --pca $bfile --grm $bfile --axes 20 
ldak5.linux --calc-pca-loads $bfile --grm $bfile --bfile $bfile --pcastem $bfile

## do projection
wdir="/GTExV8/genotypes"
bfile="$wdir/allchr.maf5hwe6missing10.kgp.ld2.nodis"
kgpdir="/1000G/GRCh38"
kgpbfile="$kgpdir/allchr.snps.maf001.hwe.geno9.biallelic.gtex.ld2.nodis"
/nfs/research1/stegle/users/na/bin/ldak5.linux --calc-scores $bfile --scorefile $kgpbfile.load --bfile $bfile --allow-flips YES 
