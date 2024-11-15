library(LDlinkR)
trans=read.delim("final_gtex_nuctraneqtls.txt")
allassocs=NULL
n=ceiling(nrow(trans)/50) ## query 50 at a time 
for (i in 2:n-1){
        ii=i*50+1
        ij=(i+1)*50
        if (ij > 260){
        ij = 260
        }
        print(c(i,ii,ij))
        SNP=trans$snp_id[ii:ij]
        a=LDtrait(SNP,
                pop = "CEU",
                r2d = "r2",
                r2d_threshold = 0.8, ## LD r2 for GWAS proxy SNP
                win_size = 1000000,
                token = "tokenid",
                file = FALSE,
                genome_build = "grch38"
        )
        allassocs=rbind(allassocs,a)
}

write.table(allassocs,"final_gtex_nuctraneqtls_gwas_ld8.txt",sep="\t",quote=F,row.names=F,col.names=T)
