# ldsc.zslope.plot
# - Generate the scatter plot of z-scores against enrichment percentile, with line of best fit
#' @param res The LDSC results dataframe
#' @param ct Name of the celltype to plot. If NULL, uses facet_wrap to show all the celltypes
#' @param gwas Filename of the GWAS study used
#' @param file Nametag of the celltype_data file
#' @param Psize Size of the points
#' @param orderBY (optional) Only used if ct=NULL (showing all celltypes as facets). Determines whether to sort the facets by z-score or slope. Set "z" for z-score or "slope" for slope.
#' @param plotType (optional) While originally intended to show Z-Score, it can also be used to show Enrichment (with confidence intervals) or Coefficient (with confidence intervals). Should be either "z", "coef" or "enrichment"
ldsc.zslope.plot <- function(res,ct,gwas,file,Psize,orderBY="z",plotType="z",colourLastDecile=FALSE){
    if(!plotType %in% c( "z", "coef", "enrichment")){
        stop("ERROR: plotType passed to ldsc.zslope.plot must be either: z, coef or enrichment")
    }
    if(!is.null(ct)){
        # IF ONLY SHOWING A SINGLE CELLTYPE
        subRes = res[res$gwas==gwas & res$file==file & res$celltype==ct,]
        subRes = subRes[!is.na(subRes$percentile),]
        if(dim(subRes)[1]==0){    stop(sprintf("ERROR: there are no results entries for %s, %s and %s",gwas,file,ct))     }
        ct_name=unique(res[res$celltype==ct,]$writtenName)
    }else{
        # IF SHOWING ALL CELLTYPES
        subRes = res[res$gwas==gwas & res$file==file,]
        subRes = subRes[!is.na(subRes$percentile),]
        if(dim(subRes)[1]==0){    stop(sprintf("ERROR: there are no results entries for %s and %s",gwas,file))     }
        # Divide celltype writtenNames if too long
        for(wN in unique(subRes$writtenName)){  subRes$writtenName[subRes$writtenName==wN] = paste(strwrap(wN, width=17),collapse="\n")  }
        # Find metric for ordering celltypes
        if(orderBY!="p"){
            subRes = get.celltype.order(subRes,orderBY,decreasing=TRUE)
        }else{
            subRes = get.celltype.order(subRes,orderBY,decreasing=FALSE)
        }
        # Apply that order to the celltypes
        #ct_mapping = unique(subRes[,c("celltype","writtenName")])
        #rownames(ct_mapping)=ct_mapping$celltype
        #subRes$celltype = factor(as.character(subRes$celltype),levels=ct_order)
        #subRes$writtenName = factor(as.character(subRes$writtenName),levels=ct_mapping[ct_order,]$writtenName)
        
    }
    subRes$percentile_colour = subRes$percentile
    if(colourLastDecile==TRUE){subRes$percentile_colour[subRes$percentile<1]=0}
    facLEVELS=unique(sort(a.n(subRes$percentile_colour)))
    colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))
    subRes$coef_ci = subRes$Coefficient_std_error*1.96
    subRes$enrich_ci = subRes$Enrichment_std_error*1.96
    subRes$dotColor = "black"
    subRes$dotColor[subRes$percentile==1] = "red"
    
    library(cowplot)
    if(plotType=="z"){
        p2 = ggplot(subRes,aes(x=percentile,y=Coefficient_z.score))+geom_point(pch=21,aes(fill=factor(percentile_colour,levels=facLEVELS)),size=Psize)+
            xlab("Specificity Decile")+ylab("Z-Score") + geom_smooth(method = "lm", se = TRUE)+
           # graph_theme +scale_fill_brewer(palette="Paired")+ theme(legend.position="none") + scale_x_continuous(breaks=unique(sort(a.n(subRes$percentile))))+
            scale_fill_manual(values=colfunc(50))+ theme(legend.position="none") + scale_x_continuous(breaks=unique(sort(a.n(subRes$percentile))))+
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }
    if(plotType=="enrichment"){
        p2 = ggplot(subRes,aes(x=percentile,y=Enrichment))+geom_point(pch=21,aes(fill=factor(percentile_colour,levels=facLEVELS)),size=Psize)+
            xlab("Specificity Decile")+ylab("Enrichment") + geom_smooth(method = "lm", se = TRUE, alpha=0.5)+
            # graph_theme +scale_fill_brewer(palette="Paired")+ theme(legend.position="none") + scale_x_continuous(breaks=unique(sort(a.n(subRes$percentile))))+
            scale_fill_manual(values=colfunc(50))+ theme(legend.position="none") + scale_x_continuous(breaks=unique(sort(a.n(subRes$percentile))),labels=c("X","N",as.character(1:10)))+
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
        if(colourLastDecile){
            colNames = c("black","red"); names(colNames)=c("0","1")
            p2 = p2 + geom_errorbar(aes(ymin=Enrichment-enrich_ci, ymax=Enrichment+enrich_ci, colour=factor(percentile_colour,levels=facLEVELS)), width=.1, position="dodge")#  + scale_color_manual(values=colNames)
        }else{    
            p2 = p2 + geom_errorbar(aes(ymin=Enrichment-enrich_ci, ymax=Enrichment+enrich_ci), colour="black", width=.1, position="dodge") 
        }
    }
    if(plotType=="coef"){
        p2 = ggplot(subRes,aes(x=percentile,y=Coefficient))+geom_point(pch=21,aes(fill=factor(percentile_colour,levels=facLEVELS)),size=Psize)+
            xlab("Specificity Decile")+ylab("Coefficient") + geom_smooth(method = "lm", se = TRUE)+
            # graph_theme +scale_fill_brewer(palette="Paired")+ theme(legend.position="none") + scale_x_continuous(breaks=unique(sort(a.n(subRes$percentile))))+
            scale_fill_manual(values=colfunc(50))+ theme(legend.position="none") + scale_x_continuous(breaks=unique(sort(a.n(subRes$percentile))))+
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_errorbar(aes(ymin=Coefficient-coef_ci, ymax=Coefficient+coef_ci), colour="black", width=.1, position="dodge") 
    }    
    
    # IF ONLY SHOWING ONE CELLTYPES... 
    if(!is.null(ct)){
        p2 = p2 + ggtitle(ct_name)
    }
    # IF SHOWING ALL CELLTYPES... apply facets
    if(is.null(ct)){
        p2 = p2+facet_wrap(~writtenName)+scale_colour_brewer(palette="Dark2")
    }
    #p2 = p2 + scale_x_discrete(labels=0:10)
    return(list(plot=p2,subRes=subRes))
}