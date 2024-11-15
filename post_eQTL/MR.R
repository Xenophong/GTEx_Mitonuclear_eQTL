## read in GWAS sumstats at nuc trans-eQTLs 
gwas=read.delim("/final_gtex_nuctraneqtls_gwas_ld8.txt")
## convert OR into BETA
getlci<-function(x){
lci=as.numeric(unlist(strsplit(x,"-",fixed=T,perl=F))[1])
return(lci)
}
gwas$LCI=sapply(gwas$Beta_or_OR,getlci)
or=which(gwas$Effect_Size_95_CI>1)
gwas$Effect_Size_95_CI[or]=log(gwas$Effect_Size_95_CI[or])
gwas$LCI[or]=log(gwas$LCI[or])
gwas$SE=(gwas$Effect_Size_95_CI-gwas$LCI)/1.96
## keep the top LD proxy when there are many 
gwasuniq=NULL
for (snp in unique(gwas$Query)){
        print(snp)
        sgwas=gwas[which(gwas$Query==snp),]
        phenos=unique(sgwas$GWAS_Trait)
        for (pheno in phenos){
                print(pheno)
                keep=sgwas[which(sgwas$GWAS_Trait==pheno),]
                keep=keep[which(keep$R2==max(keep$R2)),]
                keep=keep[which(keep$P_value==min(keep$P_value)),]
                keep=keep[1,]
                gwasuniq=rbind(gwasuniq,keep)
        }
}
write.table(gwasuniq, "/final_gtex_nuctraneqtls_gwas_ld8topsnp.txt",sep="\t",quote=F, row.names=F,col.names=T)

## read in nuc trans-eQTLs 
trans=read.delim("/final_gtex_nuctraneqtls.txt")
withgwas=trans[which(trans$snp_id%in%gwas$Query),]

## MR 
library(MendelianRandomization)
eQTLsnps=c()
tissues=c()
genes=c()
xes=c()
xses=c()
xps=c()
assocsnps=c()
ld=c()
phenos=c()
yes=c()
yses=c()
yps=c()
mr_p=c()
mr_estimate=c()
mr_se=c()
mr_f=c()

for (snp in unique(gwasuniq$Query)){
        eqtls=which(withgwas$snp_id==snp)
        print(snp)
        print(eqtls)
        for (i in eqtls){
                tissue=withgwas$tissue[i]
                gene=withgwas$gene_name[i]
                xe=withgwas$beta[i]
                xse=withgwas$se[i]
                xp=withgwas$p_value[i]
                assocs=which(gwasuniq$Query==snp&!is.na(gwasuniq$SE))
                print(assocs)
                for (j in assocs){
                        print(j)
                        assocsnp=gwasuniq$RS_Number[j]
                        pheno=gwasuniq$GWAS_Trait[j]
                        ye=gwasuniq$Effect_Size_95_CI[j]
                        yse=gwasuniq$SE[j]
                        yp=gwasuniq$P_value[j]
                        ## fill in results
                        eQTLsnps=c(eQTLsnps,snp)
                        tissues=c(tissues,tissue)
                        genes=c(genes,gene)
                        xes=c(xes,xe)
                        xses=c(xses,xse)
                        xps=c(xps, xp)
                        assocsnps=c(assocsnps,assocsnp)
                        ld=c(ld,gwasuniq$R2[j])
                        phenos=c(phenos,pheno)
                        yes=c(yes,ye)
                        yses=c(yses,yse)
                        yps=c(yps,yp)
                        ## fill in MR results
                        mr=mr_ivw(mr_input(bx = xe, bxse = xse, by = ye, byse = yse))
                        mr_p=c(mr_p,mr$Pvalue)
                        mr_estimate=c(mr_estimate,mr$Estimate)
                        mr_se=c(mr_se,mr$StdError)
                        mr_f=c(mr_f,mr$Fstat)
                }
        }
}
results=data.frame(eQTLsnps,tissues,genes,xes,xses,xps,assocsnps,ld,phenos,yes,yses,yps,mr_p,mr_estimate,mr_se,mr_f)
results$annot=annot$annotation[match(results$tissues,annot$tissue)]
results=results[,c("annot","tissues","genes","eQTLsnps","xes","xses","xps","assocsnps","ld","phenos","yes","yses","yps","mr_p","mr_estimate","mr_se","mr_f")]
colnames(results)=c("annot","tissue","eqtl_gene","eqtl_snp_id","eqtl_beta","eqtl_beta_se","eqtl_p","gwas_snp_id","gwas_snp_ld","gwas_pheno","gwas_beta","gwas_se","gwas_p","mr_p","mr_estimate","mr_se","mr_f")
write.table(results, "/final_gtex_nuctraneqtls_gwas_ld8topsnp_MR.txt",sep="\t",quote=F, row.names=F,col.names=T)

## do multiple testing correction across all studies accessible 
library(gwasrapidd)
allphenos=0
allstudies=0
for (esnp in results$eqtl_snp_id){
        assocsnps=gwas$RS_Number[which(gwas$Query==esnp)]
        astudies=NULL
        aphenos=NULL
        ## get all studies that have any of the assocsnps
        for (assocsnp in assocsnps){
                studies=get_studies(variant_id = assocsnp, warnings = FALSE)
                astudy=studies@studies$study_id
                phenos=studies@studies$reported_trait
                astudies=union(astudies,astudy)
                aphenos=union(aphenos,phenos)
        }
        nstudy=length(astudies)
        npheno=length(aphenos)
        allstudies=allstudies+nstudy
        allphenos=allphenos+npheno
}
# > allphenos
# [1] 1492
# > allstudies
# [1] 1976
sigthreshold=0.05/allphenos
results=read.delim("./final/final_gtex_nuctraneqtls_gwas_ld8topsnp_MR.txt",as.is=T)
results$sig="NO"
results$sig[which(results$mr_p<sigthreshold)]="YES"
write.table(results, "/final_gtex_nuctraneqtls_gwas_ld8topsnp_MR.txt",sep="\t",quote=F, row.names=F,col.names=T)
results$eQTL=paste(results$tissue,results$eqtl_snp_id,results$eqtl_gene,sep=".")
sig=results[which(results$sig=="YES"),]

