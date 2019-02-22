plot_slopePval_graphs <- function(slope_res,saveFolder,p_Thresh = 0.001){
    slope_res$log10p = log10(slope_res$p)
    slope_res$Significance = sprintf("P > %s",p_Thresh)
    slope_res$Significance[slope_res$p<p_Thresh] = sprintf("P <  %s",p_Thresh)
    myColors = c("red","deepskyblue3")
    names(myColors) = c(sprintf("P <  %s",p_Thresh),sprintf("P > %s",p_Thresh))
    #Some test data
    #dat <- data.frame(x=runif(10),y=runif(10),
    #                  grp = rep(LETTERS[1:5],each = 2),stringsAsFactors = TRUE)
    #Create a custom color scale
    #library(RColorBrewer)
    #myColors <- brewer.pal(5,"Set1")
    #names(myColors) <- levels(dat$grp)
    #colScale <- scale_colour_manual(name = "grp",values = myColors)
    for(gwas in as.character(unique(slope_res$gwas))){
        for(file in as.character(unique(slope_res$file))){
            subRes = slope_res[slope_res$file==file & slope_res$gwas==gwas,]
            subRes$writtenName = factor(subRes$writtenName,levels=subRes$writtenName[order(subRes$p,decreasing=TRUE)])
            # Determine appropriate width/height
            longest=get.length.of.longest.string(as.character(unique(subRes$writtenName)))
            width=4+((8/3)/45)*longest # Account for width of cell names
            height=(6.5*length(unique(subRes$celltype))/35)+1 # Account for number of cells
            # Generate the graph
            base_graph = ggplot(subRes,aes(x=writtenName,fill=Significance))+xlab("Celltype")+graph_theme+coord_flip()+graph_theme
            pdf(file=sprintf("%s/%s_%s.pdf",saveFolder,gwas,file),width=width,height=height)
            print(base_graph + geom_bar(aes(y=log10p),stat="identity")+coord_flip()+ylab("log10 P-value")+xlab("Celltype")+theme(legend.position = c(0.25, 0.18))+scale_fill_manual(values=myColors))
            dev.off()        
        }
    }
}