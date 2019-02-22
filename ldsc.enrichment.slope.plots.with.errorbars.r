#' ldsc_celltype_quantile_enrichments_wtControls
#'
#' \code{ldsc_celltype_quantile_enrichments_wtControls} Generates a plot showing:
#' (1) Density plot of the proportion statistic for 'celltype1' in 'sct_file'
#' (2) Z-Score scatter plot with line of best fit for each quantile gene set of celltype1 in gwas1
#' (3) Z-Score scatter plot with line of best fit for each quantile gene set of celltype2 in gwas2
#'
#' @param all_res Output of load_ldsc_res_files()
#' @param gene_groups Output of load_ldsc_gene_groups()
#' @param sct_file Name of a celltype_data file, i.e. "celltype_data_allKImouse_level1_thresh0_trim0"
#' @param celltype1 The celltype whose proportion density will be plotted, and the enrichments for which will be shown in inset1
#' @param celltype2 The celltype whose enrichments for which will be shown in inset2
#' @param plotType (optional) While originally intended to show Z-Score, it can also be used to show Enrichment (with confidence intervals) or Coefficient (with confidence intervals). Should be either "z", "coef" or "enrichment"
#' @param gwas1
#' @param gwas2 
#' the base name for the celltype_data file. For instance, the following name:
#'     "GeneGroups_celltype_data_allKImouse_level1_thresh0_trim0_CutOff0_percentile_StepSize_0.1.csv"
#' Can be reduced to: "celltype_data_allKImouse_level1_thresh0_trim0" using "GeneGroups_|_CutOff0_percentile_StepSize_0.1.csv"
#' @return A list of lists containing the quantile group of each gene in the 'cell_groups' dataframe, and also the maximum 
#' proportion within each quantile group within 'max_prop_inGroup'.
#' @examples
#' ldsc_celltype_quantile_enrichments_wtControls(all_res=all_gwas,gene_groups=gene_groups,sct_file="celltype_data_allKImouse_level1_thresh0_trim0",celltype1="MediumSpinyNeuron",gwas1="schiz.qjecp.clozukpgc2.gz",gwas2="GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt",gwas1_title="Schizophrenia",gwas2_title="Height",folder=folder)
#' @export
ldsc.enrichment.slope.plots.with.errorbars <- function(all_res,gene_groups,sct_file="",celltype1,celltype2,gwas1="",gwas2="",gwas1_title="",gwas2_title="",folder,maxY=NA,plotType="z"){
    if(!plotType %in% c( "z", "coef", "enrichment")){
        stop("ERROR: plotType passed to ldsc.zslope.plot must be either: z, coef or enrichment")
    }
    if(sum(c(celltype1,celltype2) %in% names(gene_groups[[sct_file]]))!=2){stop("ERROR: celltype1 or celltype2 is not found within  names(gene_groups[[sct_file]])")}

    # all_res=all_gwas
    # gene_groups=gene_groups
    # sct_file="celltype_data_allKImouse_level1_thresh0_trim0"
    # celltype1="MediumSpinyNeuron"
    # celltype2="MediumSpinyNeuron"#"endothelial-mural"
    # gwas1="schiz.qjecp.clozukpgc2.gz"
    # gwas2="GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt"
    # gwas1_title="Schizophrenia"
    # gwas2_title="Height"
    
    print(celltype1)
    print(celltype2)
    
    ccIN = gene_groups[[sct_file]][[celltype1]]$cell_groups
    dens <- density(ccIN$proportion) #density(dt$y)
    df <- data.frame(x=dens$x, y=dens$y)
    probs <- seq(from=0.1,to=1,by=0.1)#c(0.1, 0.25, 0.5, 0.75, 0.9)
    quantiles <- quantile(ccIN$proportion, prob=probs)
    df$quant <- factor(findInterval(df$x,quantiles))
    #display.brewer.all() 
    ct_name=unique(all_res[all_res$celltype==celltype1 & all_res$file==sct_file,]$writtenName)
    all_res$percentile[is.na(all_res$percentile)]=-0.1
    df$quant=as.numeric(as.character(df$quant))
    df$quant[1]=-2
    df$quant[2]=-1
    df$quant=as.factor(df$quant)
    p1 = ggplot(df, aes(x,y)) + geom_line() + geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) + scale_x_continuous() + 
        #scale_fill_brewer(guide="none",palette="Paired") + xlab("Proportion of Expression in Celltype") +
        scale_fill_brewer(guide="none",palette="Paired") + xlab("Specificity of Expression") + 
        ylab("Density")+graph_theme+ggtitle(ct_name)
    #p2 = ldsc.zslope.plot(all_gwas,cc,gwas="schiz.qjecp.clozukpgc2.gz",fName,Psize=6)
    #p3 = ldsc.zslope.plot(all_gwas,"astrocytes_ependymal",gwas="schiz.qjecp.clozukpgc2.gz",fName,Psize=6)
    relData = all_res[all_res$celltype %in% c(celltype1,celltype2) & all_res$gwas %in% c(gwas1,gwas2),]
    if(plotType=="z"){
        minY = min(relData$Coefficient_z.score)-0.2
        maxY = max(relData$Coefficient_z.score)+0.2
    }
    if(plotType=="coef"){
        border = max(relData$Coefficient_std_error)*2
        y1dat = relData[relData$gwas==gwas1 & relData$celltype==celltype1 & relData$file==sct_file,]
        border = max(y1dat$Coefficient_std_error)*2
        minY1 = min(y1dat$Coefficient)-border
        maxY1 = max(y1dat$Coefficient)+border
        y2dat = relData[relData$gwas==gwas2 & relData$celltype==celltype2 & relData$file==sct_file,]
        border = max(y2dat$Coefficient_std_error)*2
        minY2 = min(y2dat$Coefficient)-border
        maxY2 = max(y2dat$Coefficient)+border        
    }    
    if(plotType=="enrichment"){
        border = max(relData$Enrichment_std_error)*2
        minY = min(relData$Enrichment)-border
        maxY = max(relData$Enrichment)+border
    }        
    # res,ct,gwas,file,Psize,orderBY="z"
    p2 = ldsc.zslope.plot(res=all_res,ct=celltype1,gwas=gwas1,file=sct_file,Psize=6,plotType=plotType)
    if(plotType=="coef"){minY=minY1;maxY=maxY1}
    p2$plot = p2$plot + ggtitle(gwas1_title) + scale_fill_brewer(guide="none",palette="Paired") + ylim(c(minY,maxY))
    p3 = ldsc.zslope.plot(res=all_res,ct=celltype2,gwas=gwas2,file=sct_file,Psize=6,plotType=plotType)
    if(plotType=="coef"){minY=minY2;maxY=maxY2}
    p3$plot = p3$plot + ggtitle(gwas2_title) + scale_fill_brewer(guide="none",palette="Paired") + ylim(c(minY,maxY))
    #yRange=range(c(p2$subRes_ct$Coefficient_z.score,p3$subRes_ct$Coefficient_z.score))
    #yRange[1]=yRange[1]-0.2
    #yRange[2]=yRange[2]+0.2
    #if(yRange[2]<4){yRange[2]=4}
    #if(yRange[1]>-2){yRange[1]=-2}
    #p2$plot=p2$plot+ylim(yRange)
    #p3$plot=p3$plot+ylim(yRange)
    #p1 ; print(p2,vp=viewport(.8, .75, .4, 0.4) )
    # X2;Y2;width;height
    
    pd_folder = sprintf("%s/Figures/PropDensity/",folder)
    if(!file.exists(pd_folder)){  dir.create(pd_folder)  }
    pd_folderF = sprintf("%s/Figures/PropDensity/%s",folder,sct_file)
    if(!file.exists(pd_folderF)){  dir.create(pd_folderF)  }
    
    newFNAME = sprintf("%s/Fig_PropDensity_%s_%s_%s_%s_%s.pdf",pd_folderF,plotType,celltype1,gwas1_title,celltype2,gwas2_title)
    #pdf(file=newFNAME,width=8,height=3.5)
    pdf(file=newFNAME,width=6,height=3.5)
    #print(p1) ; print(p2$plot,vp=viewport(.4, .65, .3, .6) ) ; print(p3$plot,vp=viewport(.8, .65, .3, .6) )
    print(p1) ; print(p2$plot,vp=viewport(.4, .55, .4, .6) ) ; print(p3$plot,vp=viewport(.8, .55, .4, .6) )
    dev.off()
    
    library(cowplot)
    newFNAME = sprintf("%s/Fig_PropDensity_%s_%s_%s_%s_%s_3Plot.pdf",pd_folderF,plotType,celltype1,gwas1_title,celltype2,gwas2_title)
    pdf(file=newFNAME,width=9,height=3)
    print(plot_grid(p1,p2$plot,p3$plot,labels = c("A", "B", "C"),ncol=3))
    dev.off()
    
    newFNAME = sprintf("%s/Fig_LogPropDensity_%s_%s_%s_%s_%s_3Plot.pdf",pd_folderF,plotType,celltype1,gwas1_title,celltype2,gwas2_title)
    pdf(file=newFNAME,width=9,height=3)
    print(plot_grid(p1+scale_x_log10(),p2$plot,p3$plot,labels = c("A", "B", "C"),ncol=3))
    dev.off()
}