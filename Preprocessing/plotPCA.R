## in R plot the projections
kgp=read.table("/1000G/GRCh38/allchr.snps.maf001.hwe.geno9.biallelic.gtex.ld2.nodis.vect",as.is=T)
pop=read.table("/1000G/GRCh38/samples.populations.txt",header=T,as.is=T)
score=read.table("/GTExV8/genotypes/allchr.maf5hwe6missing10.kgp.ld2.nodis.profile",head=T,as.is=T)
colnames(kgp)=c("Sample", "FID", paste("PC", 1:20,sep=""))
kgp$pop=pop$Population[match(kgp$Sample,pop$Sample)]
kgp$col=NA
pops=levels(as.factor(kgp$pop))
plotcol=rainbow(n=nlevels(as.factor(kgp$pop)))
plotcol=rev(plotcol)
for (p in 1:length(pops)){
	kgp$col[which(kgp$pop==pops[p])]=plotcol[p]
}

kgp$superpop=pop$Superpopulation[match(kgp$Sample,pop$Sample)]
kgp$supercol=NA
kgp$supercol[which(kgp$superpop=="AFR")]="#e41a1c"
kgp$supercol[which(kgp$superpop=="EUR")]="#377eb8"
kgp$supercol[which(kgp$superpop=="ASN")]="#4daf4a"
kgp$supercol[which(kgp$superpop=="SAS")]="#984ea3"
kgp$supercol[which(kgp$superpop=="AMR")]="#ff7f00"
super=c("AFR","EUR","ASN","SAS","AMR")
supercol=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")

pdf("/nfs/research1/stegle/stegle_secure/GTExV8/genotypes/allchr.maf5hwe6missing10.kgp.ld2.nodis.projections.pdf",height=8,width=8)
par(mfrow=c(2,2))
##first plot all
plot(kgp[,3:4],col=kgp$supercol)
points(score[,c(5,7)],col="black")
legend('topleft', legend = super,col = supercol, cex = 0.8, pch = 1, ncol=2)

##then plot all EUR
eurpops=unique(kgp$pop[which(kgp$superpop=="EUR")])
popind=match(eurpops,pops)
eurcols=plotcol[popind]
plot(kgp[,3:4],col=kgp$col,ylim=c(-0.027,-0.015),xlim=c(0.016,0.024))
points(score[,c(5,7)],col="black")
legend('bottom', legend = eurpops,col = eurcols, cex = 0.8, pch = 1, ncol=5)

##then plot all ASN
asnpops=unique(kgp$pop[which(kgp$superpop=="ASN")])
popind=match(asnpops,pops)
asncols=plotcol[popind]
plot(kgp[,3:4],col=kgp$col,ylim=c(0.029,0.037),xlim=c(0.002,0.006))
points(score[,c(5,7)],col="black")
legend('bottom', legend = asnpops,col = asncols, cex = 0.8, pch = 1, ncol=5)

##then plot all AFR
afrpops=unique(kgp$pop[which(kgp$superpop=="AFR")])
popind=match(afrpops,pops)
afrcols=plotcol[popind]
plot(kgp[,3:4],col=kgp$col,ylim=c(-0.02,-0.005),xlim=c(-0.04,-0.005))
points(score[,c(5,7)],col="black")
legend('bottom', legend = afrpops,col = afrcols, cex = 0.8, pch = 1, ncol=4)
dev.off()

eur=score[which(score$Profile2<=-0.019&score$Profile1>=0.017),c("ID1", "ID2")]
write.table(eur, "/GTExV8/genotypes/eursample.keep", sep = "\t", quote=F, row.names=F, col.names=F )
