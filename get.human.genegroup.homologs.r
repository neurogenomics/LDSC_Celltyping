get.human.genegroup.homologs <- function(gg,mouse_to_human_homologs){
    # If the gene set already has hgnc_symbol, just drop any without homologs
    if(sum(colnames(gg)=="hgnc_symbol")==1){ gg = gg[gg$hgnc_symbol %in% mouse_to_human_homologs$HGNC.symbol,] }
    if(sum(colnames(gg)=="mgi_symbol")==1){ 
        m2h = unique(mouse_to_human_homologs[,c("MGI.symbol","HGNC.symbol")])
        colnames(m2h)[colnames(m2h)=="MGI.symbol"]="mgi_symbol"
        gg = merge(gg,m2h,by="mgi_symbol")
        gg = gg[gg$HGNC.symbol!="",]
        gg = gg[,c("HGNC.symbol","proportion","percentile")]
        colnames(gg)[1]="hgnc_symbol"
    }
    return(gg)
}