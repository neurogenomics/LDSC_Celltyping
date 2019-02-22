# gg is the output of load_ldsc_gene_groups()
get.celltype.annot <- function(gg){
    library(stringi)
    files=names(gg)
    count=0
    for(ff in files){
        count=count+1
        tmp = data.frame(file=ff,celltypes=names(gg[[ff]]),writtenName="")
        tmp$writtenName = get.written.cellnames(files=tmp$file,cellnames=tmp$celltypes)
        if(count==1){
            print("creating annew")
            celltype_annot = tmp
        }else{
            print("appending")
            celltype_annot = rbind(celltype_annot,tmp)
        }
    }
    write.csv(celltype_annot,file="celltype_annot.csv")
    return(celltype_annot)
}