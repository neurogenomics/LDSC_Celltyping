#!/usr/bin/env Rscript
# file_name="celltype_data_allKImouse_level1_thresh0_trim0"
args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
sep <- args[2]
if(sep=="tab"){sep="\t"}

library("data.table")
gwasSumStats = fread(file_name,stringsAsFactors=FALSE,sep=sep,fill=TRUE)

#snps = fread(sprintf("/%s/MAGMA/g1000_eur.bim",Sys.getenv("HOME")),stringsAsFactors=FALSE)
#save(snps,file=sprintf("/%s/MAGMA/g1000_eur.bim.rda",Sys.getenv("HOME")))
load(file=sprintf("/%s/MAGMA/g1000_eur.bim.rda",Sys.getenv("HOME")))
colnames(snps)=c("CHR","SNP","GeneticDistance","BP","A1","A2")
snps = snps[,c("SNP","CHR","BP")]
#colnames(adhd)[1]="MarkerName"
snps$MarkerName=sprintf("%s:%s",snps$CHR,snps$BP)
gwasSumStats2 = merge(gwasSumStats,snps,by="MarkerName")
keepCols = setdiff(colnames(gwasSumStats),c("Direction","MarkerName"))
dropCols = setdiff(colnames(gwasSumStats2),c("SNP","CHR","BP",keepCols))
gwasSumStats3 = gwasSumStats2[, c("SNP","CHR","BP",keepCols),with=FALSE]
#gwasSumStats3 = setcolorder(gwasSumStats2, keepCols)

colnames(gwasSumStats3) = gsub("Weight","N",colnames(gwasSumStats3))
colnames(gwasSumStats3) = gsub("ref","A1",colnames(gwasSumStats3))
colnames(gwasSumStats3) = gsub("alt","A2",colnames(gwasSumStats3))
colnames(gwasSumStats3) = gsub("Allele1","A1",colnames(gwasSumStats3))
colnames(gwasSumStats3) = gsub("Allele2","A2",colnames(gwasSumStats3))
colnames(gwasSumStats3) = gsub("P.value","P",colnames(gwasSumStats3))
colnames(gwasSumStats3) = gsub("P-value","P",colnames(gwasSumStats3))

gwasSumStats4 = gwasSumStats3
gwasSumStats4$P = as.numeric(gwasSumStats4$P)
gwasSumStats4 = gwasSumStats4[!is.na(gwasSumStats4$P),]

write.table(gwasSumStats4,file=file_name,quote=FALSE,row.names=FALSE)
