#' plot.gwas.diff
#'
#' \code{plot.gwas.diff} Generates scatter plots showing the difference in celltype enrichments between two GWAS studies. These plots
#' are saved to file and returned from the function. One of the plots is just a scatter plot with a line of best fit. The second plot
#' names all celltypes with z-scores over 2 that have residuals greater than 1 from the line of best fit.
#'
#' @param all_res Output of load_ldsc_res_files()
#' @param sct_file Name of a celltype_data file, i.e. "celltype_data_allKImouse_level1_thresh0_trim0"
#' @param gwas1 The GWAS study that will form the x-axis (as denoted in all_res)
#' @param gwas2 The GWAS study that will form the y-axis (as denoted in all_res)
#' @param gwas1_title A publishable name for the GWAS study to use on the graph axes
#' @param gwas2_title A publishable name for the GWAS study to use on the graph axes
#' @param folder The root directory for the results analysis. Should contain the 'Figures' folder.
#' @return A list of scatter plots. Also saves two figures as pdf.
#' @examples
#' plot.gwas.diff(all_res=all_gwas,file=ff,gwas1="schiz.qjecp.clozukpgc2.gz",gwas2="scz2.snp.results.txt",
#'      gwas1_title="Schizophrenia (CLOZUK)",gwas2_title="Schizophrenia (108 loci)",folder=folder)
#' @export
plot.gwas.diff <- function(all_res,sct_file,gwas1,gwas2,gwas1_title="",gwas2_title="",folder,plotType="z"){
    library(ggrepel)
    
    # Check arguments
    if(sum(c(gwas1,gwas2) %in% unique(all_res$gwas))<2){stop("ERROR: either gwas1 or gwas2 not found in all_res")}
    if(gwas1_title==""){gwas1_title=gwas1}
    if(gwas2_title==""){gwas2_title=gwas2}
    if(gwas1==gwas2){stop("ERROR: cannot compare a GWAS study to itself")}
    #if(!file.exists(sprintf("%s/Figures/EnrichmentComparisons",folder))){stop("ERROR: 'folder' must contain Figures/EnrichmentComparisons subfolder")}
    
    # Create new directory for 'file'
    newDir=sprintf("%s/Figures/EnrichmentComparisons/%s",folder,sct_file)
    #newDir=sprintf("%s/%s",folder)
    if(!file.exists(newDir)){dir.create(newDir)}
    
    a.c<-function(x){as.character(x)}
    aR = all_res[!is.na(all_res$percentile),]
    dat1 = aR[aR$file==sct_file & aR$gwas==gwas1 & aR$percentile==1,]
    dat2 = aR[aR$file==sct_file & aR$gwas==gwas2 & aR$percentile==1,]
    rownames(dat1) = dat1$celltype
    rownames(dat2) = dat2$celltype
    cts = unique(c(a.c(dat1$celltype),a.c(dat1$celltype)))
    dat1=dat1[cts,]
    dat2=dat2[cts,]
    #datBoth = data.frame(celltype=dat1$celltype,gwas1_z=dat1$Coefficient_z.score,gwas2_z=dat2$Coefficient_z.score)
    datBoth = data.frame(celltype=dat1$writtenName,gwas1_z=dat1$Coefficient_z.score,gwas2_z=dat2$Coefficient_z.score,gwas1_p=dat1$Enrichment_p,gwas2_p=dat2$Enrichment_p)
    if(plotType=="z"){    mod=lm(gwas2_z~gwas1_z,data=datBoth) }
    if(plotType=="p"){    mod=lm(gwas2_p~gwas1_p,data=datBoth) }
    residuals = residuals(mod)
    names(residuals) = datBoth$celltype
    datBoth$residuals = residuals
    datBoth$label = ""
    if(plotType=="z"){    showLabel = abs(datBoth$residuals)>1 & (datBoth$gwas1_z>2 | datBoth$gwas2_z>2)    }
    #if(plotType=="p"){    showLabel = abs(datBoth$residuals)>1 & (datBoth$gwas1_p<0.05 | datBoth$gwas2_p<0.05)    }
    if(plotType=="p"){    showLabel = (datBoth$gwas1_p<0.05 | datBoth$gwas2_p<0.05)    }
    datBoth$label[showLabel] = as.character(datBoth$celltype)[showLabel]
    datBoth$sL = "Not-signif"
    datBoth$sL[showLabel] = "Signif"

    xL = sprintf("%s\nz-score",gwas1_title)
    yL = sprintf("%s\nz-score",gwas2_title)
    Rlabel = sprintf("R = %.2f",cor(datBoth$gwas1_z,datBoth$gwas2_z))
    
    if(plotType=="z"){    
        thePlot = ggplot(datBoth,aes(x=gwas1_z,y=gwas2_z,label=label))+xlab(gwas1)+ylab(gwas2) + annotate("rect", xmin = -8, xmax = 2, ymin = -8, ymax = 2, alpha = .2)+
            coord_cartesian(xlim=c(min(datBoth$gwas1_z)-0.2,max(datBoth$gwas1_z)+0.2),ylim=c(min(datBoth$gwas2_z)-0.2,max(datBoth$gwas2_z)+0.2))+
            geom_smooth(method = "lm", se = FALSE)+
            geom_point(aes(colour=sL))    +    xlab(xL) + ylab(yL)+        graph_theme+
            geom_abline(intercept=mod$coefficients[1]+1,slope=mod$coefficients[2])+
            geom_abline(intercept=mod$coefficients[1]-1,slope=mod$coefficients[2])+ 
            scale_colour_manual(name="",values = c("Not-signif"="black", "Signif"="red"))+ theme(legend.position="none")+
            geom_text(aes(label=Rlabel, x = -1.3, y = 2.1))
    }else{
        minPlim = min(datBoth$gwas1_p,datBoth$gwas2_p)/10
        thePlot = ggplot(datBoth,aes(x=gwas1_p,y=gwas2_p,label=label))+xlab(gwas1)+ylab(gwas2) + scale_y_log10() + scale_x_log10() +
            geom_point() + geom_smooth(method = "lm", se = TRUE) + coord_cartesian(xlim=c(minPlim,1),ylim=c(minPlim,1)) + 
            annotate("rect", xmin = minPlim, xmax = 0.05, ymin = minPlim, ymax = 0.05, alpha = .2)
            
        
        yMin = min(datBoth$gwas1_p,datBoth$gwas2_p)/10
        thePlot = thePlot
    }
    
    thePlot = thePlot 
    
    fName=sprintf("%s/UNLABELED_%s_%s.pdf",newDir,gwas1,gwas2)
    pdf(file=fName,width=4,height=4);    print(thePlot);    dev.off()
    
    fName=sprintf("%s/LABELED_%s_%s.pdf",newDir,gwas1,gwas2)
    thePlotLabelled = thePlot
    if(sum(showLabel)>0){   thePlotLabelled = thePlot+geom_label_repel(fill="yellow",box.padding = unit(1, "lines"),point.padding = unit(1.5, "lines"))    }
    pdf(file=fName,width=4,height=4);    print(thePlotLabelled);    dev.off()
    return(list(unlabelled=thePlot,labelled=thePlotLabelled))
}