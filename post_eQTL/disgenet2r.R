## disgenet 
library(disgenet2r)
api_key="yourkey"
Sys.setenv(DISGENET_API_KEY= api_key)

## read in mtDNA trans-eQTLs
trans=read.delim("/final_gtex_mttranseqtls.txt",as.is=T)
genes=unique(trans$ensembl_gene_id)
results=gene2disease(
       genes,
       vocabulary = "ENSEMBL",
       n_pags = 100,
       api_key = api_key,
       chemical = NA,
       database = "CURATED",
       score = c(0.5, 1),
       verbose = FALSE,
       warnings = TRUE
     )
tab=results@qresult
tab=tab[,c("gene_symbol","ensemblid","geneNcbiType","geneDSI","geneDPI","genepLI","protein_class_name","disease_name","diseaseType","diseaseClasses_MSH","score","yearInitial","yearFinal","evidence_index")]

## get diseases classes
getclasses<-function(x){
	y=paste0(unlist(x),collapse=", ")
	return(y)
}
tab$disease_MSH=sapply(tab$diseaseClasses_MSH,getclasses)
tab$protein_class_name=sapply(tab$protein_class_name,getclasses)
tab$ensemblid=sapply(tab$ensemblid,getclasses)

## get disease codes
getcodes<-function(x){
	ys=unlist(x)
	cs=c()
	for (y in ys){
		c=regmatches(y, gregexpr("(?<=\\().*?(?=\\))", y, perl=T))[[1]]
		cs=c(cs,c)
	}
	cs=paste0(unlist(cs),collapse=",")
	return(cs)
}
tab$disease_MSH_codes=sapply(tab5$diseaseClasses_MSH,getcodes)
allcodes = c()
uniqcodes = c()
for (codes in tab5$disease_MSH_codes){
	codes = unlist(strsplit(codes,",",fixed=T,perl=F))
	allcodes = c(allcodes,codes)
	uniqcodes = union(uniqcodes,codes)
}

write.table(tab,"./final_gtex_mttranseqtls_disgenet.txt",sep="\t",quote=F,row.names=F,col.names=T)
