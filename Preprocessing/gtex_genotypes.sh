## get genotypes at common variants from VCF in GTEx V8
vcfdir="/GTExV8/ncbi/files/phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1"
wdir="/GTExV8"
mkdir $wdir/genotypes
mkdir $wdir/genotypes/src
mkdir $wdir/genotypes/log
chromosomes=$(seq 1 22)
for chr in $chromosomes 
do echo $chr
plink2 \
--vcf $vcfdir/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz \
--maf 0.05 --hwe 10e-6 --geno 0.1 --snps-only --chr $chr --make-bed \
--out $wdir/genotypes/chr$chr.maf5hwe6missing10
done 

## merge data from all chromosomes
wdir="/GTExV8/genotypes"
rm $wdir/merge-list
chromosomes=$(seq 2 22)
bfile=$wdir/chr1.maf5hwe6missing10
for chr in $chromosomes
do echo $chr 
echo -e "$wdir/chr$chr.maf5hwe6missing10.bed $wdir/chr$chr.maf5hwe6missing10.bim $wdir/chr$chr.maf5hwe6missing10.fam">> $wdir/merge-list
done 
plink --bfile $bfile --merge-list $wdir/merge-list --make-bed --out $wdir/allchr.maf5hwe6missing10

## get relatedness using KING 
bfile=$wdir/allchr.maf5hwe6missing10
king -b $bfile.bed --kinship --related --degree 3
