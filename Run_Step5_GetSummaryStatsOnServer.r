#!/bin/sh

# Set parameters:
#folder_set = "Results_Phase6"
#folder_set = "Results_Phase 7-May2017"
#folder_set = "Results-Phase12-NewDataSet"
folder_set = "Results-Phase29-AllanMTG-lvl1"
thresh=0
stepSize=0.1

# Create base folder for Summarised Results (if it doesn't exist)
sumResFolder = file.path(Sys.getenv("LDSCSRC"), "SummaryResults")
dir.create(sumResFolder, showWarnings = FALSE)
folder=file.path(sumResFolder, folder_set)
dir.create(folder, showWarnings = FALSE)
resfolder=file.path(folder, "ResFiles")
dir.create(folder, showWarnings = FALSE)

#cd $LDSCSRC
#source $LDSCCODE/InitVars.sh



library(ggplot2)
library(gridExtra)
library(grid)
split_mode = "percentile"
a.n<-function(x){as.numeric(as.character(x))}

# mv *allKImouse_MergedStriatal* ../Results_Phase6/ResFiles/
# mv *TASIC_1to1only_level1* ../Results_Phase6/ResFiles/ 


#folder=sprintf("/lustre/scratch117/unknown/team141/ns9/LDSC/%s",folder_set)
setwd(folder)
geneGroupFolder=sprintf("%s/GeneGroups",Sys.getenv("LDSCSRC"))
gwas_annot = read.csv(sprintf("%s/gwas_set.csv",Sys.getenv("LDSCCODE")),stringsAsFactors = FALSE)
rownames(gwas_annot)=gwas_annot$Filename
gwas_set = trimws(gwas_annot$Filename)

ewce_folder = Sys.getenv("LDSCCODE")
source(sprintf("%s/graph_theme.r",ewce_folder))
source(sprintf("%s/ldsc.zslope.plot.r",ewce_folder))
source(sprintf("%s/load_ldsc_res_files.r",ewce_folder))
source(sprintf("%s/plot.gwas.diff.r",ewce_folder))
source(sprintf("%s/get.celltype.annot.r",ewce_folder))
source(sprintf("%s/get.comparison.hgnc.genesets.r",ewce_folder))
source(sprintf("%s/get.human.genegroup.homologs.r",ewce_folder))
source(sprintf("%s/load_ldsc_gene_groups.r",ewce_folder))
source(sprintf("%s/get.written.cellnames.r",ewce_folder))

#source(sprintf("%s/ldsc_celltype_quantile_enrichments_wtControls.r",ewce_folder))
source(sprintf("%s/get.length.of.longest.string.r",ewce_folder))
source(sprintf("%s/get.celltype.order.r",ewce_folder))
source(sprintf("%s/plot_slopePval_graphs.r",ewce_folder))
source(sprintf("%s/ldsc.enrichment.slope.plots.with.errorbars.r",ewce_folder))


# CREATE DIRECTORIES FOR FIGURES
dir.create(sprintf("%s/Figures",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures/SlopePval/",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures/EnrichmentZscores_BySlope/",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures/EnrichmentPVals_BySlope/",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures/EnrichmentZscores_ByZSCORE/",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures/EnrichmentPVals_ByZSCORE/",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures/Enrichment_ByP/",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures/Celltype_ZSCORE/",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures/Celltype_PVALS/",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures/Celltype_SLOPE/",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures/EnrichmentComparisons/",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures/EnrichmentHeatmaps/",folder),showWarnings=FALSE)
dir.create(sprintf("%s/Figures/CoefficientPlots/",folder),showWarnings=FALSE)

# Aggregate the results files
#allGwasFiles = list.files(sprintf("%s/Results",Sys.getenv("LDSCSRC")))
#for(gwas in allGwasFiles){
#}
system(sprintf("cp -r %s %s",sprintf("%s/Results",Sys.getenv("LDSCSRC")),resfolder))

#for(stepSize in c(0.1,0.025)){
for(stepSize in c(0.1)){
    # LOAD THE LDSC RESULTS, THEN SAVE TO FILE
    all_gwas = load_ldsc_res_files(resPath=folder,gwas_set=gwas_set,validFileGrep=sprintf("Thresh0.Stepsize%s.results",stepSize),stepSize=stepSize)
    all_gwas$stepSize=stepSize
    all_gwas = all_gwas[order(all_gwas$Enrichment_p),]
    write.csv(all_gwas,file=sprintf("%s/LDSC_GWAS_Results_stepSize_%s.csv",folder,stepSize))
    save(all_gwas,file=sprintf("%s/LDSC_GWAS_Results_stepSize_%s.Rda",folder,stepSize))
    
    # PROCESS THE GENEGROUP FILES, DROPPING THE UNNECCESARY GENES
    gene_groups = load_ldsc_gene_groups(ggPath=sprintf("%s",geneGroupFolder))
    save(gene_groups,file=sprintf("%s/gene_groups_stepSize_%s.rda",folder,stepSize,stepSize),compress="xz",compression_level=9)
}

originalCT=all_gwas$celltype
# Add 'writtenNames' to all_gwas
for(StepSize in c(0.1)){
    load(file=sprintf("%s/LDSC_GWAS_Results_stepSize_%s.Rda",folder,stepSize))
    load(file=sprintf("%s/gene_groups_stepSize_%s.rda",folder,stepSize))
    celltype_annot = get.celltype.annot(gene_groups)
    all_gwas$str = sprintf("%s_%s",all_gwas$file,all_gwas$celltype)
    celltype_annot$str = sprintf("%s_%s",celltype_annot$file,celltype_annot$celltypes) 
    cta = celltype_annot[,c("writtenName","str")]
    all_gwas = merge(all_gwas,cta,by="str")
    all_gwas$gwas=as.character(all_gwas$gwas)
}


  
# library(cowplot)
# for(gwas in unique(all_gwas$gwas)){
#     sub_gwas = all_gwas[all_gwas$gwas==gwas,]
#     all_files = unique(sub_gwas$file)
#     dir.create(sprintf("%s/Figures/CoefficientPlots/%s",folder,gwas),showWarnings=FALSE)
#     for(ff in all_files){
#         sub_folder = sprintf("%s/Figures/CoefficientPlots/%s/%s",folder,gwas,ff)
#         dir.create(sub_folder,showWarnings=FALSE)
#         FigSX_SlopeEnrichmentsAllCells = ldsc.zslope.plot(res=all_gwas,ct=NULL,gwas=gwas,file=ff,Psize=1,orderBY="p",plotType="enrichment",colourLastDecile=TRUE)
#         pdf(file=sprintf("%s/SlopeEnrichmentsAllCells_%s_%s.pdf",sub_folder,gwas,ff),width=10,height=8)
#         print(FigSX_SlopeEnrichmentsAllCells$plot)
#         dev.off()
#     }
# }

# unique(all_gwas[all_gwas$file=="celltype_data_DroncSeqOldFormat_level1",]$celltype)
# for(pT in c("z","coef","enrichment")){
#     ldsc.enrichment.slope.plots.with.errorbars(all_res=all_gwas,gene_groups=gene_groups,sct_file="celltype_data_allKImouse_wtHypo_MergedStriatal_1to1only_level1_thresh0_trim0",celltype1="MediumSpinyNeuron",celltype2="MediumSpinyNeuron",gwas1="Schiz_Clozuk.txt",gwas2="GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",gwas1_title="Schizophrenia",gwas2_title="Height",folder=folder,plotType=pT)
#     ldsc.enrichment.slope.plots.with.errorbars(all_res=all_gwas,gene_groups=gene_groups,sct_file="celltype_data_allKImouse_wtHypo_MergedStriatal_1to1only_level1_thresh0_trim0",celltype1="endothelial-mural",celltype2="endothelial-mural",gwas1="Schiz_Clozuk.txt",gwas2="GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",gwas1_title="Schizophrenia",gwas2_title="Height",folder=folder,plotType=pT)
#     ldsc.enrichment.slope.plots.with.errorbars(all_res=all_gwas,gene_groups=gene_groups,sct_file="celltype_data_allKImouse_wtHypo_MergedStriatal_1to1only_level1_thresh0_trim0",celltype1="astrocytes_ependymal",celltype2="astrocytes_ependymal",gwas1="Schiz_Clozuk.txt",gwas2="GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",gwas1_title="Schizophrenia",gwas2_title="Height",folder=folder,plotType=pT)
#     ldsc.enrichment.slope.plots.with.errorbars(all_res=all_gwas,gene_groups=gene_groups,sct_file="celltype_data_allKImouse_wtHypo_MergedStriatal_1to1only_level1_thresh0_trim0",celltype1="pyramidalCA1",celltype2="pyramidalCA1",gwas1="Schiz_Clozuk.txt",gwas2="GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",gwas1_title="Schizophrenia",gwas2_title="Height",folder=folder,plotType=pT)
#     ldsc.enrichment.slope.plots.with.errorbars(all_res=all_gwas,gene_groups=gene_groups,sct_file="celltype_data_allKImouse_wtHypo_MergedStriatal_1to1only_level1_thresh0_trim0",celltype1="pyramidalSS",celltype2="pyramidalSS",gwas1="Schiz_Clozuk.txt",gwas2="GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",gwas1_title="Schizophrenia",gwas2_title="Height",folder=folder,plotType=pT)
#     ldsc.enrichment.slope.plots.with.errorbars(all_res=all_gwas,gene_groups=gene_groups,sct_file="celltype_data_DroncSeqOldFormat_level1",celltype1="GABA2",celltype2="GABA2",gwas1="Schiz_Clozuk.txt",gwas2="GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",gwas1_title="Schizophrenia",gwas2_title="Height",folder=folder,plotType=pT)
#     ldsc.enrichment.slope.plots.with.errorbars(all_res=all_gwas,gene_groups=gene_groups,sct_file="celltype_data_DroncSeqOldFormat_level1",celltype1="exPFC1",celltype2="exPFC1",gwas1="Schiz_Clozuk.txt",gwas2="GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",gwas1_title="Schizophrenia",gwas2_title="Height",folder=folder,plotType=pT)
#     ldsc.enrichment.slope.plots.with.errorbars(all_res=all_gwas,gene_groups=gene_groups,sct_file="celltype_data_DroncSeqOldFormat_level1",celltype1="OPC",celltype2="OPC",gwas1="Schiz_Clozuk.txt",gwas2="GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",gwas1_title="Schizophrenia",gwas2_title="Height",folder=folder,plotType=pT)
#     ldsc.enrichment.slope.plots.with.errorbars(all_res=all_gwas,gene_groups=gene_groups,sct_file="celltype_data_DroncSeqOldFormat_level1",celltype1="END",celltype2="END",gwas1="Schiz_Clozuk.txt",gwas2="GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",gwas1_title="Schizophrenia",gwas2_title="Height",folder=folder,plotType=pT)
# }



all_gwas$Enrichment_p[all_gwas$Enrichment>1] <- all_gwas$Enrichment_p[all_gwas$Enrichment>1]/2
all_gwas$Enrichment_p[all_gwas$Enrichment<1] <- 1-(all_gwas$Enrichment_p[all_gwas$Enrichment<1]/2)
minP = min(all_gwas$Enrichment_p)
minQ = min(p.adjust(all_gwas$Enrichment_p,method="BH"))

library(cowplot)
files = unique(all_gwas$file)
for(gwas in gwas_set[gwas_set %in% as.character(unique(all_gwas$gwas))]){
    for(file in files[files %in% all_gwas$file]){
        set_name = sprintf("%s_%s",gwas,file)
        subRes = all_gwas[all_gwas$gwas==gwas & all_gwas$file==file & !is.na(all_gwas$percentile),]

        subRes$Enrichment_q = p.adjust(subRes$Enrichment_p,method="BH")

        subMinQ = minQ
        if(min(subRes$Enrichment_q)<subMinQ){subMinQ=min(subRes$Enrichment_q)}	

        if(dim(subRes)[1]==0){
            print(sprintf("Warning: No data available for %s and %s",file,gwas))
        }else{
            #q_Thresh = 0.05
            subRes$Significance = sprintf("Q > %s",0.05)
            subRes$Significance[subRes$Enrichment_q<=0.05] = sprintf("Q <= %s",0.05)
            subRes$Significance[subRes$Enrichment_q<=0.005] = sprintf("Q <= %s",0.005)
            subRes$Significance[subRes$Enrichment_q<=0.0005] = sprintf("Q <= %s",0.0005)

            manualColours = c("blue4","blue","blueviolet","red")	
            names(manualColours) = c("Q > 0.05","Q <= 0.05","Q <= 0.005","Q <= 0.0005")

            #possibleBreaks = sort(5/10^(1:50),decreasing=TRUE)
            

            # Find the length of longest celltype name, so width of graphs can be set
            longest=get.length.of.longest.string(as.character(unique(subRes$writtenName)))

            if(dim(subRes)[1]>1){
                width=4+((8/3)/45)*longest # Account for width of cell names
                height=(6.5*length(unique(subRes$celltype))/35)+1

                # JUST PLOT Z-SCORES
                subRes = get.celltype.order(subRes,"z")
                #subRes$celltype = factor(as.character(subRes$celltype),levels=get.celltype.order(subRes,"z"))
                base_graph = ggplot(subRes[subRes$percentile==1,],aes(x=writtenName,fill=Significance))+xlab("Celltype")+graph_theme+coord_flip() + scale_fill_manual(values=manualColours)
                pdf(file=sprintf("%s/Figures/Celltype_ZSCORE/%s_%s.pdf",folder,gwas,file),width=width,height=height)
                print(base_graph + geom_bar(aes(y=Coefficient_z.score),stat="identity")+coord_flip()+ylab("Z-Score")+xlab("Celltype")+theme(legend.position = c(0.75, 0.18)))
                dev.off()

                # JUST PLOT PVALS
                subRes = get.celltype.order(res=subRes,orderBY="p",decreasing=TRUE)
                #subRes$celltype = factor(as.character(subRes$celltype),levels=get.celltype.order(subRes,"p",decreasing=TRUE))
                base_graph = ggplot(subRes[subRes$percentile==1,],aes(x=writtenName,fill=Significance))+xlab("Celltype")+graph_theme+coord_flip() + scale_fill_manual(values=manualColours) #+ylim(limits = c(subMinQ, 1))#+
                	#scale_y_continuous()
                pdf(file=sprintf("%s/Figures/Celltype_PVALS/%s_%s.pdf",folder,gwas,file),width=width,height=height)
                print(base_graph + geom_bar(aes(y=Enrichment_q),stat="identity")+ylab(expression(log[10](p)))+scale_y_log10(limits=c(subMinQ,1))+theme(legend.position = c(0.25, 0.18)))
                dev.off()

                # JUST PLOT SLOPE
                subRes = get.celltype.order(res=subRes,orderBY="slope",decreasing=FALSE)
                #subRes$celltype = factor(as.character(subRes$celltype),levels=get.celltype.order(subRes,"p",decreasing=TRUE))
                base_graph = ggplot(subRes[subRes$percentile==1,],aes(x=writtenName,fill=Significance))+xlab("Celltype")+graph_theme+coord_flip() + scale_fill_manual(values=manualColours)
                pdf(file=sprintf("%s/Figures/Celltype_SLOPE/%s_%s.pdf",folder,gwas,file),width=width,height=height)
                print(base_graph + geom_bar(aes(y=slope),stat="identity")+ylab("Slope of Enrichment")+theme(legend.position = c(0.25, 0.18)))
                dev.off()
            }
        }
    }
}

# FIT LINEAR MODEL
count=0
for(gwas in gwas_set[gwas_set %in% as.character(unique(all_gwas$gwas))]){
    for(file in files[files %in% all_gwas$file]){
        subRes1 = all_gwas[all_gwas$file==file & all_gwas$gwas==gwas,]
        if(dim(subRes1)[1]==0){
            print(sprintf("Warning: No data available for %s and %s",file,gwas))
        }else{
            for(ct in as.character(unique(subRes1$celltype))){
                subRes = subRes1[subRes1$celltype==ct,]
                subRes[is.na(subRes$percentile),]$percentile=-0.1
                subRes$half = "Top"
                subRes$half[subRes$percentile<=0.5]="Bottom"
                wilcox_p = wilcox.test(subRes[subRes$half=="Top",]$Coefficient_z.score,subRes[subRes$half=="Bottom",]$Coefficient_z.score)$p.value
                count=count+1
                mod = lm(Coefficient_z.score~percentile,data=subRes)
                p = anova(mod)$"Pr(>F)"[1]
                slope = mod$coefficients[2]
                if(slope>0){p=p/2}
                if(slope<0){p=1-p/2}
                slope = lm(Coefficient_z.score~percentile,data=subRes)$coefficients[2]
                tmp = data.frame(gwas=gwas,file=file,celltype=ct,stepSize=stepSize,p=p,slope=slope,writtenName=unique(subRes$writtenName),wilcox_p=wilcox_p)
                if(count==1){
                    slope_res = tmp
                }else{slope_res = rbind(slope_res,tmp)}
            }
        }
    }
}
save(slope_res,file=sprintf("%s/slope_res.rda",folder))
slope_res = slope_res[order(slope_res$p),]
slope_res$q=p.adjust(slope_res$p,method="BH")
#slope_res$p[slope_res$slope<0]=1
write.csv(slope_res,file=sprintf("%s/slope_res_stepSize_%s.csv",folder,stepSize))
slope_res_q = slope_res
slope_res_q$p=slope_res_q$q
plot_slopePval_graphs(slope_res,sprintf("%s/Figures/SlopePval/",folder,stepSize),p_Thresh = 0.05)
