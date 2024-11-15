## check overlap between SNPs in the GTEX WGS data and 1000G
kgdir="/1000G/GRCh38"
wdir="/GTExV8"
gtex=read.table(paste0(wdir,"/genotypes/allchr.maf5hwe6missing10.bim"),as.is=T)
gtex$snp=paste(gtex$V1,gtex$V4,sep=":")
chromosomes=c(1:22)
for (chr in chromosomes){
	print(chr)
	gtexchr=gtex[which(gtex$V1==chr),]
	kgp=read.table(paste0(kgdir,"/chr",chr,".snps.maf001.hwe.geno9.biallelic.bim"),as.is=T)
	kgp$snp=paste(kgp$V1,kgp$V4,sep=":")
	gtexkgpsnps=intersect(gtex$snp,kgp$snp)
	gtexchr=gtexchr[which(gtexchr$snp %in% gtexkgpsnps),]
	dups=which(table(gtexchr$snp)>=2)
	if (length(dups)>=1){
		gtexchr=gtexchr[-which(gtexchr$snp%in%names(dups)),]
		kgp=kgp[-which(kgp$snp%in%names(dups)),]
	}
	kgp=kgp[which(kgp$snp %in% gtexkgpsnps),]
	gtexchr=gtexchr[order(gtexchr$V4,decreasing=F),]
	kgp=kgp[order(kgp$V4,decreasing=F),]
	gtexchr$rsid=kgp$V2
	write.table(gtexchr[,c(1:6)],paste0(wdir,"/genotypes/chr",chr,".maf5hwe6missing10.kgp.extract"),quote=F, row.names=F, col.names=F, sep="\t")
	write.table(gtexchr[,c(1,8,3,4,5,6)],paste0(wdir,"/genotypes/chr",chr,".maf5hwe6missing10.kgp.rsids"),quote=F, row.names=F, col.names=F, sep="\t")
}
