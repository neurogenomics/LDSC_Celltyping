detect.species <- function(genenames,mouse_to_human_homologs,opt){
	numHuman = sum(genenames %in% mouse_to_human_homologs$HGNC.symbol)
	numMouse = sum(genenames %in% mouse_to_human_homologs$MGI.symbol)
	species = "mouse"
	if(numHuman>numMouse){species="human"}
	opt$species = species
	return(opt)
}		


merge.genePercentiles <- function(folder_path){
	input_files=c("celltype_data_allKImouse_level1_thresh0_trim0","celltype_data_allKImouse_level2_thresh0_trim0","celltype_data_allKImouse_noLevels_thresh0_trim0","celltype_data_cortex_and_hippocampus_KI_mouse_level1_thresh0_trim0","celltype_data_cortex_and_hippocampus_KI_mouse_level2_thresh0_trim0","celltype_data_striatal_and_midbrain_KI_mouse_level1_thresh0_trim0","celltype_data_striatal_and_midbrain_KI_mouse_level2_thresh0_trim0","celltype_data_TASIC_level1_thresh1_trim0","celltype_data_TASIC_level2_thresh1_trim0","celltype_data_SCAPT_noLevels_NOthresh")	
	step_sizes=c(0.1,0.01)
	split_mode="percentile"
	
	if(length(grep("KI",input_file))>0){cutoff_thresh=0.5}
	if(length(grep("SCAPT",input_file))>0){cutoff_thresh=30}
	if(length(grep("TASIC",input_file))>0){cutoff_thresh=15}
	
	for(input_file in input_files){
		for(step_size in step_sizes){
			iFiles = list.dirs(path=folder_path,full.names=FALSE)
			dirPattern=sprintf("%s.*_CutOff%s_%s_StepSize_%s.*",input_file,cutoff_thresh,split_mode,step_size)
			iFiles = iFiles[grep(dirPattern,iFiles)]
			iFiles = gsub(sprintf("%s/",folder_path),"",iFiles)
			
			celltype_count=0
			for(ff in iFiles){
				celltype_count=celltype_count+1
				celltype = gsub(sprintf("%s_AnnotLevel_1_",input_file),"",ff)
				celltype = gsub(sprintf("_CutOff%s_%s_StepSize_%s.*",cutoff_thresh,split_mode,step_size),"",celltype)
				ffName = sprintf("%s_AnnotLevel_1_%s_CutOff%s_%s_StepSize_%s/genePercentiles.csv",input_file,celltype,cutoff_thresh,split_mode,step_size)
				geneGroupTable = cbind(celltype=celltype,read.csv(ffName,stringsAsFactors=FALSE))
				
				if(celltype_count==1){
		   			allGeneGroups=geneGroupTable
		   		}else{
		   			allGeneGroups=rbind(allGeneGroups,geneGroupTable)
		   		}
			}
			write.csv(allGeneGroups,file=sprintf("%s/GeneGroups_%s_CutOff%s_%s_StepSize_%s.csv",folder_path,input_file,cutoff_thresh,split_mode,step_size))
		}	
	}
}

check.ct.names <- function(folder_path,opt){
	load(sprintf("%s%s.rda",opt$ctd_folder,opt$input_file))
	#colnames(celltype_data[[1]]$all_scts) = gsub(";|\\(|\\)| |\\/","",colnames(celltype_data[[1]]$all_scts))
	colnames(celltype_data[[1]]$all_scts) = gsub(";|\\(|\\)| |\\/|,|+|-","",colnames(celltype_data[[1]]$all_scts))
    #colnames(celltype_data[[1]]$cell_dists) = gsub(";|\\(|\\)| |\\/","",colnames(celltype_data[[1]]$cell_dists))
    colnames(celltype_data[[1]]$cell_dists) = gsub(";|\\(|\\)| |\\/|,|+|-","",colnames(celltype_data[[1]]$cell_dists))
    outfile=sprintf("%s%s.rda",opt$ctd_folder,opt$input_file)
    print(outfile)
    print("Corrected celltype names:")
    print(colnames(celltype_data[[1]]$cell_dists))
	save(celltype_data,file=outfile)
	print("Cell type names corrected succesfully")
}

generate.go.enrichment.files <- function(folder_path,opt){
	allGeneGroups = read.csv(file=sprintf("%s/GeneGroups_%s_CutOff%s_%s_StepSize_%s.csv",folder_path,opt$input_file,opt$cutoff_thresh,opt$split_mode,opt$step_size),stringsAsFactors=FALSE)
    # Test GO enrichments for specific gene sets
    celltype_count=0
    for(target_celltype in as.character(unique(allGeneGroups$celltype))){
    	celltype_count=celltype_count+1
    	subDat = allGeneGroups[allGeneGroups$celltype==target_celltype,]
    	GOres = cbind(celltype=target_celltype,run.go.analysis(subDat))
   		if(celltype_count==1){
   			allGOres=GOres
   		}else{
   			allGOres=rbind(allGOres,GOres)
   		}    	
    }
	write.csv(allGOres,file=sprintf("%s/GOres_%s_CutOff%s_%s_StepSize_%s.csv",folder_path,opt$input_file,opt$cutoff_thresh,opt$split_mode,opt$step_size))
}


get.mouse_to_human_homologs <- function(opt){
	load(sprintf("%s/mouse_to_human_homologs.rda",opt$base_folder))
	#library("biomaRt")
	#human = useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
	#mouse = useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
	#mgi_symbols = getBM(attributes=c("mgi_symbol"), mart=mouse)
	#mouse_to_human_homologs = getLDS(attributes = c("mgi_symbol","entrezgene","ensembl_gene_id"),
	#     filters = "mgi_symbol", values = mgi_symbols,
	#     mart = mouse,
	#     attributesL = c("hgnc_symbol","ensembl_gene_id","entrezgene"), martL = human)
	
	# If the single cell dataset is actually human... then trick the code into treating it as mouse
	#if(length(grep("SCAPT",opt$input_file))>0){
	#	uniq = unique(rownames(celltype_data[[annot_level]]$cell_dists))
	#	mouse_to_human_homologs = data.frame(MGI.symbol=uniq,HGNC.symbol=uniq)
	#}	
	
	return(mouse_to_human_homologs)
}

prepare.new.template.matchingGenes.files <- function(path,base,chr){
	print(sprintf("Preparing chromosome %s",chr))
	
	# Load the massive SNP database
	# NOTE:
	# - Loading the actual SNP database is unrealistic, so preprocess it first:
	# - First download it:
	# -   wget ftp://ftp.ncbi.nih.gov/snp/database/organism_data/human_9606/b147_SNPContigLocusId_107.bcp.gz
	# - Then remove the SNP and Gene ID columns:
	# -   gzip -cd b147_SNPContigLocusId_107.bcp.gz | cut -f 1,7 > SNPPartialLocus.csv
	# - For running tests, if the code is modified, grab the first 50000 lines of SNPPartialLocus.csv
	#     head -n 50000 SNPPartialLocus.csv > SNPPartialLocus_head.csv
	print("Reading the SNP database")
	all.snps <- read.csv("/lustre/scratch117/unknown/team141/ns9/LDSC/SNPPartialLocus.csv", as.is = TRUE,sep="\t", header=FALSE, stringsAsFactors=FALSE)
	print("Finished reading the SNP database")
	
	
	template_files = list.files(sprintf("%s",path))
	template_files = template_files[grep(".annot.gz",template_files)]
	numSteps=10000
	lenSNPS = length(all.snps$V2)
	#quantile(1:lenSNPS)
	percentiles = c(1,lenSNPS * (1:numSteps/numSteps))

			
    # Load the template file
    fffName=sprintf("%s%s",path,sprintf(base,chr))
    print(sprintf("Loading %s",fffName))
    template.file <- read.csv(gzfile(fffName), as.is = TRUE,sep="\t")
    template.snps <- gsub("rs","",template.file$SNP)
    count=0
    countJ=0
    keptSNPS = all.snps[1:100000,]
    #tmpSTORE=list()
    for(i in 1:(numSteps-1)){
    #for(i in 1:50){
        print(sprintf("chr: %s, i: %s",chr,i))
        #print(unique(as.character(all.snps$V2[(percentiles[i]-200):percentiles[i]])))
        #all.snps$V2[(percentiles[i]-20):percentiles[i]]
        keepSNPS = all.snps$V1[percentiles[i]:percentiles[i+1]]  %in% as.numeric(template.snps)
        numKept = sum(keepSNPS)
        #print(numKept)
    #}
        if(numKept>0){
            # To avoid having to create a new dataframe each time, sequentially fill a 1mil long matrix
            if(count+numKept>100000){
                if(countJ==0){
                    allKeptSNPS = keptSNPS[1:count,]
                }else{
                    allKeptSNPS = rbind(allKeptSNPS,keptSNPS[1:count,])
                }
                countJ=countJ+1
                count=0
            }
            keptSNPS[(count+1):(count+numKept),] = all.snps[percentiles[i]:percentiles[i+1],][keepSNPS,] 
            #print(as.character(unique(all.snps[percentiles[i]:percentiles[i+1],][keepSNPS,2] )))
            #if(count==0){
            #    keptSNPS = all.snps[percentiles[i]:percentiles[i+1],][keepSNPS,]
            #}else{
            #    keptSNPS = rbind(keptSNPS,all.snps[percentiles[i]:percentiles[i+1],][keepSNPS,])
            #}
            count=count+numKept
        }
    }
    allKeptSNPS = rbind(allKeptSNPS,keptSNPS[1:count,])
    
    templateName = sprintf("%s%s",path,sprintf(base,chr))
    newFName = gsub("annot.gz","matchingGenes.rda",templateName)
    save(allKeptSNPS,file=newFName)
}


add.col.to.template <- function(template.file,inputSNPS,colName){
	inAnnot = rep(0,dim(template.file)[1])
	inAnnot[template.snps %in% inputSNPS] = 1
    template.file = cbind(template.file,inAnnot)
    colnames(template.file)[dim(template.file)[2]] = colName
    return(template.file)
}

# tempFileOrig=template.file

##############################################################################################################################
# find.snps.without.genes
	# Find all SNPs which don't match to either:
	#  (A) gene regions [that is, the SNP is in an intergenic region]
	#  (B) genes with mouse-->human homologs
# - Inputs:
#  * allKeptSNPS: data.frame with two columns, V1 storing SNP IDs, V2 storing HGNC symbols
#  * mouse_to_human_homologs: dataframe storing all mouse to human homologs, NOT just those found in the dataset
find.snps.without.genes <- function(allKeptSNPS,mouse_to_human_homologs){
	SNPS_with_GENES = unique(allKeptSNPS$V1[allKeptSNPS$V2 %in% mouse_to_human_homologs$HGNC.symbol])
	SNPS_without_GENES = unique(allKeptSNPS$V1[!allKeptSNPS$V2 %in% mouse_to_human_homologs$HGNC.symbol])
	SNPS_without_GENES = setdiff(SNPS_without_GENES,SNPS_with_GENES)
	percentWithoutGenes = length(SNPS_without_GENES)/length((unique(SNPS_with_GENES,SNPS_without_GENES)))*100
	print(sprintf("%s SNPS map to genes with human<-->homologs... %s (%.0f%%) SNPS do not",length(SNPS_with_GENES),length(SNPS_without_GENES),percentWithoutGenes))	
	return(list(SNPS_with_GENES=SNPS_with_GENES,SNPS_without_GENES=SNPS_without_GENES,percentWithoutGenes=percentWithoutGenes))
}
##############################################################################################################################


##############################################################################################################################
# find.unexpressed.genes

find.unexpressed.genes <- function(cell_prop,allKeptSNPS,mouse_to_human_homologs,cutOFF=0.001,opt){
	if(opt$species=="mouse"){
		# Find all SNPS associated with genes that are... in homologous genes but:
		# (A) Not expressed in SCT_dataset_ dataset
		#Genes_NotInSCT_dataset   = as.character(unique(mouse_to_human_homologs[!as.character(mouse_to_human_homologs$MGI.symbol) %in% names(cell_prop),"HGNC.symbol"]))
		Genes_NotInSCT_dataset   = as.character(unique(mouse_to_human_homologs[!as.character(mouse_to_human_homologs$MGI.symbol) %in% names(cell_prop),"MGI.symbol"]))
		SCT_dataset_Homologs   = as.character(unique(mouse_to_human_homologs[mouse_to_human_homologs$MGI.symbol %in% names(cell_prop),"MGI.symbol"]))
		human.SCT_dataset_Homologs   = as.character(unique(mouse_to_human_homologs[mouse_to_human_homologs$MGI.symbol %in% names(cell_prop),"HGNC.symbol"]))

		SNPS_notExp_SCT_dataset_ = unique(allKeptSNPS$V1[allKeptSNPS$V2 %in% mouse_to_human_homologs$HGNC.symbol & (!allKeptSNPS$V2 %in% human.SCT_dataset_Homologs)])
		percent_SNPS_notExp_SCT_dataset_ = length(SNPS_notExp_SCT_dataset_) / length(unique(allKeptSNPS$V1))*100
		print(sprintf("%s genes have known human homologs, but are not expressed in any cell in the single cell dataset",length(Genes_NotInSCT_dataset)))
		print(sprintf("%s SNPS (%.1f %%) map to genes with known human homologs, but are not expressed in any cell in the single cell dataset",length(SNPS_notExp_SCT_dataset_),percent_SNPS_notExp_SCT_dataset_))
	
		# (B) Not expressed in the cell type
		notExpressedGenes = names(cell_prop[cell_prop<cutOFF])
		human.notExpressedGenes = unique(mouse_to_human_homologs[mouse_to_human_homologs$MGI.symbol %in% notExpressedGenes,]$HGNC.symbol)
		expressedGenes = names(cell_prop[cell_prop>cutOFF])	
		human.ExpressedGenes = unique(mouse_to_human_homologs[mouse_to_human_homologs$MGI.symbol %in% expressedGenes,]$HGNC.symbol)	
	}else if(opt$species=="human"){
		# (A) Not expressed in SCT_dataset_ dataset
		# hgnc_symbols = getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), mart=human)
		# save(hgnc_symbols,file="hgnc_symbols.rda")
		load(sprintf("%s/hgnc_symbols.rda",opt$base_folder))
		AllHumanGenes = unique(hgnc_symbols$hgnc_symbol)
		Genes_NotInSCT_dataset   = as.character(unique(AllHumanGenes[!AllHumanGenes %in% names(cell_prop)]))
		print(sprintf("Of %s human genes, %s are not expressed in the single cell dataset",length(AllHumanGenes),length(Genes_NotInSCT_dataset)))
		
		SNPS_notExp_SCT_dataset_ = unique(allKeptSNPS$V1[allKeptSNPS$V2 %in% Genes_NotInSCT_dataset])
		SCT_dataset_Homologs=NA
		
		# (B) Not expressed in the cell type
		notExpressedGenes = names(cell_prop[cell_prop<cutOFF])
		expressedGenes = names(cell_prop[cell_prop>cutOFF])	
		human.notExpressedGenes = notExpressedGenes
		human.ExpressedGenes = expressedGenes
	}
	SNPS_Exp_inCT = allKeptSNPS$V1[allKeptSNPS$V2 %in% human.ExpressedGenes]
	SNPS_not_Exp_inCT = allKeptSNPS$V1[allKeptSNPS$V2 %in% human.notExpressedGenes]
	SNPS_not_Exp_inCT = setdiff(SNPS_not_Exp_inCT,SNPS_Exp_inCT)
	percent_SNPS_notExpInCT_SCT_dataset_ = length(unique(SNPS_not_Exp_inCT)) / length(unique(allKeptSNPS$V1))*100
	if(opt$species=="mouse"){
		print(sprintf("%s genes (with human homologs) are expressed in the single cell dataset BUT not in this celltype",length(human.notExpressedGenes)))
		print(sprintf("%s SNPS (%.1f %%) map to genes with known human homologs, are expressed in brain but not this celltype",length(SNPS_not_Exp_inCT),percent_SNPS_notExpInCT_SCT_dataset_))		
	}else if(opt$species=="human"){
		print(sprintf("%s genes are expressed in the single cell dataset BUT not in this celltype",length(human.notExpressedGenes)))
		print(sprintf("%s SNPS (%.1f %%) map to genes, are expressed in brain but not this celltype",length(SNPS_not_Exp_inCT),percent_SNPS_notExpInCT_SCT_dataset_))		
	}
	
	# Amend cell_prop to not include unexpressed genes
	cell_prop2 = cell_prop[cell_prop>cutOFF]

	# Return data
	return(list(cell_prop_reduced=cell_prop2,SNPS_not_Exp_inCT=SNPS_not_Exp_inCT,SNPS_Exp_inCT=SNPS_Exp_inCT,SNPS_notExp_SCT_dataset_=SNPS_notExp_SCT_dataset_,SCT_dataset_Homologs=SCT_dataset_Homologs,Genes_NotInSCT_dataset=Genes_NotInSCT_dataset,human.notExpressedGenes=human.notExpressedGenes,mouse.notExpressedGenes=notExpressedGenes))
}
##############################################################################################################################

##############################################################################################################################
get.mouse.id.lookup <- function(input_mgi_symbols,opt){
#	library("biomaRt")
#	mouse = useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
#	mgi_symbols = getBM(attributes=c("mgi_symbol"), mart=mouse)
#	attrib_mus = listAttributes(mouse)
#	idlookup = getBM(attributes=c("mgi_symbol","ensembl_gene_id"), filters="mgi_symbol", values=mgi_symbols, mart=mouse)		
#	save(idlookup,file="idlookup.rda")
	load(sprintf("%s/idlookup.rda",opt$base_folder))
	input_mgi_symbols = input_mgi_symbols[input_mgi_symbols!=""]	
	idlookup = idlookup[idlookup$mgi_symbol %in% input_mgi_symbols,]
	return(idlookup)
}
#############################################################################################################################

##############################################################################################################################
get.human.id.lookup <- function(input_hgnc_symbols,opt){
#	library("biomaRt")
#	mouse = useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
#	mgi_symbols = getBM(attributes=c("mgi_symbol"), mart=mouse)
#	attrib_mus = listAttributes(mouse)
#	idlookup = getBM(attributes=c("mgi_symbol","ensembl_gene_id"), filters="mgi_symbol", values=mgi_symbols, mart=mouse)		
#	save(idlookup,file="idlookup.rda")
	load(sprintf("%s/hgnc_symbols.rda",opt$base_folder))
	input_hgnc_symbols = input_hgnc_symbols[input_hgnc_symbols!=""]
	hgnc_symbols = hgnc_symbols[hgnc_symbols$hgnc_symbol %in% input_hgnc_symbols,]
	return(hgnc_symbols)
}
#############################################################################################################################


#############################################################################################################################
write.annot.desc.files <- function(cell_prop_reduced,opt,annotDescFile){
	# Get ensembl IDs
	if(opt$species=="mouse"){
		idlookup = get.mouse.id.lookup(names(cell_prop_reduced),opt)		
		cell_prop_table = data.frame(mgi_symbol=names(cell_prop_reduced),proportion=cell_prop_reduced,percentile=rep(0,length(cell_prop_reduced)))
	}else if(opt$species=="human"){
		idlookup = get.human.id.lookup(names(cell_prop_reduced),opt)		
		cell_prop_table = data.frame(hgnc_symbol=names(cell_prop_reduced),proportion=cell_prop_reduced,percentile=rep(0,length(cell_prop_reduced)))		
	}
	
	countKeptGenes = 0
	for(min_percentile in seq(from=0,to=1-opt$step_size,by=opt$step_size)){
		max_percentile = min_percentile+opt$step_size
		# Find the target genes & map target genes to human ensembl IDs	
    	prop_lims = quantile(cell_prop_reduced,probs=c(min_percentile,max_percentile))
    	target_genes = names(cell_prop_reduced[cell_prop_reduced>=prop_lims[1] & cell_prop_reduced<=prop_lims[2]])
    	# Setup values in the table
    	if(opt$species=="mouse"){
	    	cell_prop_table[cell_prop_table$mgi_symbol %in% target_genes,"percentile"]=min_percentile
    	}else if(opt$species=="human"){
    		cell_prop_table[cell_prop_table$hgnc_symbol %in% target_genes,"percentile"]=min_percentile	
    	}
	}
	
	if(opt$species=="mouse"){cell_prop_table_ensembl = merge(cell_prop_table,idlookup,by="mgi_symbol")}
	if(opt$species=="human"){cell_prop_table_ensembl = merge(cell_prop_table,idlookup,by="hgnc_symbol")}
	write.csv(cell_prop_table_ensembl,file=annotDescFile)	
	return(cell_prop_table_ensembl)
}
#############################################################################################################################


#############################################################################################################################
run.go.analysis <- function(data){
	# CHECK IF THE TOPGO PACKAGE IS ALREADY INSTALLED, IF NOT INSTALL IT
	if(!"topGO" %in% rownames(installed.packages())){
		source("https://bioconductor.org/biocLite.R")
		biocLite("topGO")
	}
	library(topGO)
	all_genes <- data$percentile
	
	names(all_genes) <- data$ensembl_gene_id
	
	GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(p) p == 0.9, description = "Test", annot = annFUN.org, mapping="org.Mm.eg.db", ID="Ensembl")
	GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = all_genes, geneSel = function(p) p == 0.9, description = "Test", annot = annFUN.org, mapping="org.Mm.eg.db", ID="Ensembl")
	GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = all_genes, geneSel = function(p) p == 0.9, description = "Test", annot = annFUN.org, mapping="org.Mm.eg.db", ID="Ensembl")
	
	resultFisher_BP <- runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")
	resultFisher_CC <- runTest(GOdata_CC, algorithm = "classic", statistic = "fisher")
	resultFisher_MF <- runTest(GOdata_MF, algorithm = "classic", statistic = "fisher")
#	resultFisher_BP <- runTest(GOdata_BP, algorithm = "weight01", statistic = "fisher")
#	resultFisher_CC <- runTest(GOdata_CC, algorithm = "weight01", statistic = "fisher")
#	resultFisher_MF <- runTest(GOdata_MF, algorithm = "weight01", statistic = "fisher")

	
	a=cbind(GenTable(GOdata_BP, classicFisher = resultFisher_BP, topNodes = 10),ontology="BP")
	b=cbind(GenTable(GOdata_CC, classicFisher = resultFisher_CC, topNodes = 10),ontology="CC")
	c=cbind(GenTable(GOdata_MF, classicFisher = resultFisher_MF, topNodes = 10),ontology="MF")
	d=rbind(a,b,c)
	d=d[order(as.numeric(d$classicFisher)),]
	return(d)
}
#############################################################################################################################


#############################################################################################################################
add.cell.to.template <- function(target_celltype,cell_prop,mouse_to_human_homologs,allKeptSNPS,template.file,template.snps,opt,annotDescFile){
	# First, find SNPs that map to outside genes, or which map to genes without mouse->human homologs
	unusedSNPS = find.snps.without.genes(allKeptSNPS,mouse_to_human_homologs)
	colName = sprintf("%s_NoHomologOrNotGene",gsub(" ","",target_celltype))	
	template.file = add.col.to.template(template.file,unusedSNPS$SNPS_without_GENES,colName)

	# Find all SNPS associated with genes that are: in homologous genes but not expressed in the celltype
	# cell_prop,allKeptSNPS,mouse_to_human_homologs,cutOFF=0.001,opt
	unExp = find.unexpressed.genes(cell_prop,allKeptSNPS,mouse_to_human_homologs,opt=opt)
	colName = sprintf("%s_NotExpressedInCelltype",gsub(" ","",target_celltype))	
	template.file = add.col.to.template(template.file,c(unExp$SNPS_notExp_SCT_dataset_,unExp$SNPS_not_Exp_inCT),colName)
	cell_prop2 = unExp$cell_prop_reduced

	# Generate a table listing which genes are in which percentile (for GO analysis etc)
	if(opt$split_mode=="percentile"){ 
		
		if(opt$species=="mouse"){
			GenesNotExpInSCT=cbind(get.mouse.id.lookup(unExp$Genes_NotInSCT_dataset,opt),proportion="GenesNotExpInSCT",percentile=-1)
			GenesNotExpInCell=cbind(get.mouse.id.lookup(unExp$mouse.notExpressedGenes,opt),proportion="GenesNotExpInCell",percentile=-1)
			#human.GenesWithoutHomologs = as.character(unique(allKeptSNPS$V2[!allKeptSNPS$V2 %in% mouse_to_human_homologs$HGNC.symbol]))
			#mouse.GenesWithoutHomologs = cbind(get.mouse.id.lookup(as.character(unique(allKeptSNPS$V2[!allKeptSNPS$V2 %in% mouse_to_human_homologs$HGNC.symbol]))),proportion="GenesWithoutHomologs",percentile=-1)
		}else if(opt$species=="human"){
			GenesNotExpInSCT=cbind(get.human.id.lookup(unExp$Genes_NotInSCT_dataset,opt),proportion="GenesNotExpInSCT",percentile=-1)
			GenesNotExpInCell=cbind(get.human.id.lookup(unExp$human.notExpressedGenes,opt),proportion="GenesNotExpInCell",percentile=-1)			
			#GenesWithoutHomologs=cbind(get.human.id.lookup(as.character(unique(allKeptSNPS$V2[!allKeptSNPS$V2 %in% mouse_to_human_homologs$HGNC.symbol]))),proportion="GenesWithoutHomologs",percentile=-1)			
		}
		#notGenes = rbind(GenesNotExpInSCT,GenesNotExpInCell,GenesWithoutHomologs)
		notGenes = rbind(GenesNotExpInSCT,GenesNotExpInCell)
		notGenes = notGenes[,c(1,3,4,2)]
		geneGroupTable = write.annot.desc.files(cell_prop_reduced = cell_prop2,opt,annotDescFile=annotDescFile) 
		geneGroupTable = rbind(geneGroupTable,notGenes)
		geneGroupTable = cbind(geneGroupTable,celltype=target_celltype)
	}

	countKeptGenes = 0
	for(min_percentile in seq(from=0,to=1-opt$step_size,by=opt$step_size)){
	    max_percentile = min_percentile+opt$step_size
	    annot_name = sprintf("%s_%s_%s",gsub(" ","",target_celltype),min_percentile,max_percentile)
	    
	    # Find the target genes & map target genes to human ensembl IDs
	    if(opt$split_mode=="percentile"){
	    	prop_lims = quantile(cell_prop2,probs=c(min_percentile,max_percentile))
	    	target_genes = names(cell_prop2[cell_prop2>=prop_lims[1] & cell_prop2<=prop_lims[2]])
	    }else if(opt$split_mode=="proportion"){
		    target_genes = names(cell_prop2[cell_prop2>=min_percentile & cell_prop2<=max_percentile])
	    }
	    
	    # Get human homologs of the target genes, then restrict to those on CHR and with SNPs
	    if(opt$species=="human"){
	    	human.targets = unique(target_genes)	    	
	    	keptGenes = as.character(unique(allKeptSNPS$V2)[unique(allKeptSNPS$V2) %in% human.targets])	    	
	    }else if(opt$species=="mouse"){
	    	human.targets = unique(mouse_to_human_homologs[mouse_to_human_homologs$MGI.symbol %in% target_genes,])
	    	keptGenes = as.character(unique(allKeptSNPS$V2)[unique(allKeptSNPS$V2) %in% human.targets$HGNC.symbol])
	    }
	    
	    # Find SNPs associated with genes in the target set
	    numKeptGenes = length(keptGenes)
	    countKeptGenes = countKeptGenes + numKeptGenes
	    target.out = unique(allKeptSNPS[allKeptSNPS$V2 %in% keptGenes,])
	    
	    # Print stats about genes and SNPs
	    print(sprintf("%s - genes with proportion %s -- %s ",target_celltype,min_percentile,max_percentile))
	    print("NUMBER OF SNPs PER GENE:")
	    print(sort(table(as.character(target.out$V2))))
	    print(sprintf("Mean number of SNPs per gene: %s",mean(table(as.character(target.out$V2)))))
	    print(sprintf("Min number of SNPs per gene: %s",min(table(as.character(target.out$V2)))))
	    
	    target.snps = unique(allKeptSNPS$V1[allKeptSNPS$V2 %in% keptGenes])
	    inAnnot = rep(0,dim(template.file)[1])
	    inAnnot[template.snps %in% target.snps] = 1
	    # If the column is empty, set 50 random SNPs to 1 to avoid numerical errors
	    if(sum(inAnnot)>0){inAnnot[sample(1:length(inAnnot),50)] = 1} 
	    
	    template.file = cbind(template.file,inAnnot)
	    colnames(template.file)[dim(template.file)[2]] = annot_name
	 }       
	 
	 # 
	 print(sprintf("Number of genes matched to SNPs in template file: %s",length(unique(allKeptSNPS$V2))))
	 print(sprintf("Number of genes used from chromosome: %s",countKeptGenes))
	 
	 return(list(template.file=template.file,geneGroupTable=geneGroupTable))
}

# Setup folder paths & variable names
print(.libPaths())
base_folder   = Sys.getenv("LDSCCODE") #"/lustre/scratch117/unknown/team141/ns9/LDSC/"
folder_path   = Sys.getenv("SNPANNOT") # sprintf("%s/snp_annot_files",base_folder)
ctd_path      = Sys.getenv("CTDFOLDER") 
#template_path = "1000g phase 1 template"
template_path = "1000G_Phase3_cell_type_groups"
template_files = list.files(sprintf("%s/%s",folder_path,template_path))
templateBaseRoot = "cell_type_group.1."
templateBaseName = paste(templateBaseRoot,"%s.annot.gz",sep="")
#template_files = template_files[grep(".annot.gz",template_files)]
template_files = template_files[grep("cell_type_group\\.1\\..*\\.annot.gz",template_files)]

# Get command line arguments
if(!"optparse" %in% rownames(installed.packages())){
	install.packages("optparse")
}
library("optparse")

option_list = list(
  make_option(c("-T", "--pre_tI"), type="integer", default=1,
              help="Index value of template annotation file to be loaded (normally 1--22)"),
  make_option(c("-C", "--cutoff_thresh"), default=0,
              help="Genes must have a total summed value of all_scts greater than the cutoff"),
  make_option(c("-M", "--split_mode"), type="character", default="percentile", 
              help="Determines whether annotations are based on percentiles of proportions, or the actual proportions [default=percentile]"),
  make_option(c("-I", "--input_file"), type="character", default="celltype_data_allKImouse_wtHypo_MergedStriatal_1to1only_level1_thresh0_trim0", 
              help="File name (without extension) of the single cell annotation file. File should be an rda containing a list with all_scts and cell_dists"),
              # celltype_data_allKImouse_level1_thresh0_trim0 celltype_data_allKImouse_level2_thresh0_trim0 celltype_data_allKImouse_noLevels_thresh0_trim0 celltype_data_cortex_and_hippocampus_KI_mouse_level1_thresh0_trim0 celltype_data_cortex_and_hippocampus_KI_mouse_level2_thresh0_trim0 celltype_data_striatal_and_midbrain_KI_mouse_level1_thresh0_trim0 celltype_data_striatal_and_midbrain_KI_mouse_level2_thresh0_trim0 celltype_data_TASIC_level1_thresh1_trim0 celltype_data_TASIC_level2_thresh1_trim0 celltype_data_SCAPT_noLevels_NOthresh
  make_option(c("-S", "--step_size"), default=0.1,
              help="Determines how many percentile groups to divide into"),
  make_option(c("-P", "--prep_template"), action="store_true", default=FALSE,
              help="Run in 'Prepare Template' mode? [default %default]"),
  make_option(c("-N", "--check_names"), action="store_true", default=FALSE,
              help="Run in 'Check Celltype Names' mode? [default %default]"),              
  make_option(c("-G", "--gen_go"), action="store_true", default=FALSE,
              help="Run in 'Generate GO enrichments' mode? [default %default]")              
)
opt = parse_args(OptionParser(option_list=option_list))
opt$base_folder = base_folder
opt$ctd_folder = sprintf("%s/",ctd_path)
opt$cell_annot_folder = sprintf("%s/CellAnnots",Sys.getenv("LDSCSRC"))
opt$tmp_folder = sprintf("%s/TMP",Sys.getenv("LDSCSRC"))
opt$gene_groups_folder = sprintf("%s/GeneGroups",Sys.getenv("LDSCSRC"))
print(opt)

# Create Tmp folder if it doesn't already exist
if(!file.exists(opt$tmp_folder)){
    dir.create(opt$tmp_folder)
}

#opt$input_file="celltype_data_allKImouse_level1_thresh0_trim0"

# - Convert chromosome number to index of template_files
# :: so pre_tI is the actual chromosome number
# :: tI is the index
tI = which(template_files==sprintf(templateBaseName,opt$pre_tI))

if(opt$prep_template){
	print("Running Prepare Template mode!")	
   	foldPath=sprintf("%s/%s/",folder_path,template_path)
   	prepare.new.template.matchingGenes.files(path=foldPath,base=templateBaseName,chr=tI)
}else if(opt$check_names){
	print("Running 'check celltype names' mode!")	
	check.ct.names(folder_path,opt)
}else if(opt$gen_go){
	print("Running 'Generate GO enrichments' mode!")	
	generate.go.enrichment.files(folder_path,opt)
}else{
	print("Running Generate Annotations mode!")	
	
	# Load EWCE proportions data
	proportion_file = opt$input_file#"celltype_data_MBB_(wtLEVELS)_thresh(0)_trim(0)"
	load(sprintf("%s%s.rda",opt$ctd_folder,proportion_file))
	annot_level = 1
	
	# Ensure no celltype has a wierd character in it's name
	a = colnames(celltype_data[[annot_level]]$all_scts)
	b = gsub(" ","",a)
	c = gsub("/","and",b)
	colnames(celltype_data[[annot_level]]$all_scts) = c
	colnames(celltype_data[[annot_level]]$cell_dists) = c
		
	# Drop genes with low numbers of total reads
	sumTotal = apply(celltype_data[[annot_level]]$all_scts,1,sum)
	celltype_data[[annot_level]]$all_scts = celltype_data[[annot_level]]$all_scts[sumTotal>=opt$cutoff_thresh,]
	celltype_data[[annot_level]]$cell_dists = celltype_data[[annot_level]]$cell_dists[sumTotal>=opt$cutoff_thresh,]
	
	# Load mouse-->human homologs
	mouse_to_human_homologs = get.mouse_to_human_homologs(opt)
		
	# Determine the species
	genenames=rownames(celltype_data[[annot_level]]$cell_dists)
	opt = detect.species(genenames=genenames,mouse_to_human_homologs,opt)
		
	print(tI)
	print(template_files[tI])

    # Load the snp-->gene file
    templateName = sprintf("%s/%s/%s",folder_path,template_path,template_files[tI])
    newFName = gsub("annot.gz","matchingGenes.rda",templateName)
    print(newFName)
    # - First, check that it actually exists 
    #    * if you've switched to a new set of SNP annot template files, it will need regenerating
    if(!file.exists(newFName)){	stop("Prepare Template mode needs to be run first. Use -PR argument.")   }
    load(newFName)	
	
	# Load the template file
    template.file <- read.csv(gzfile(sprintf("%s/%s/%s",folder_path,template_path,template_files[tI])), as.is = TRUE,sep="\t")
    #template.file <- template.file[grep("rs",template.file$SNP),]
    #template.snps <- as.numeric(gsub("rs","",template.file$SNP))
    template.snps <- as.numeric(gsub("rs|ss","",template.file$SNP))
    template.file = template.file[,-5]
    template.file.O = template.file 

	# Which SCT_dataset_ genes are not listed as having human homologs? Many miRNA, Riken genes, Fam genes etc
	#cPn = names(cell_prop)
	#cPn_noHomolog = sort(cPn[!cPn %in% mouse_to_human_homologs$MGI.symbol])
	#write.csv(cPn_noHomolog,file="cPn_noHomolog.csv")

	# Which genes in cell_prop (SCT_dataset_ dataset) do not have any SNPs (based on CNS.* template)?
	#chrNO = as.numeric(gsub(".annot.gz","",gsub("CNS.","",template_files[tI])))
	chrNO = opt$pre_tI#tI
	chrDATA = read.csv(sprintf("%s/mart_export_genesNchr.csv",opt$base_folder))
	allGenesOnCHR = as.character(unique(chrDATA[chrDATA$Chromosome.Name==chrNO,3]))
	if(opt$species=="mouse"){
		SCT_dataset_Homologs   = unique(mouse_to_human_homologs[mouse_to_human_homologs$MGI.symbol %in% genenames,"HGNC.symbol"])
		SCT_dataset_GenesOnCHR = SCT_dataset_Homologs[SCT_dataset_Homologs %in% allGenesOnCHR]
	}else if(opt$species=="human"){
		SCT_dataset_GenesOnCHR = unique(genenames[genenames %in% allGenesOnCHR])			
	}
	print(sprintf("Biomart lists %s genes as being located on this chromosome",length(allGenesOnCHR)))	
	genesWithSNPSonCHR = unique(allKeptSNPS$V2)
	SCT_dataset_GenesOnCHRwithSNPs = SCT_dataset_GenesOnCHR[SCT_dataset_GenesOnCHR %in% genesWithSNPSonCHR]
	SCT_dataset_GenesOnCHRwithoutSNPs = SCT_dataset_GenesOnCHR[!SCT_dataset_GenesOnCHR %in% genesWithSNPSonCHR]
	print(sprintf("Of %s genes from single cell dataset on chromosome, %s have SNPs",length(SCT_dataset_GenesOnCHR),length(SCT_dataset_GenesOnCHRwithSNPs)))
	print(sprintf("So %s genes from single cell dataset on chromosome do NOT have SNPs",length(SCT_dataset_GenesOnCHRwithoutSNPs)))	
		 print("Most the genes without SNPS are wierd rubbish, take a look at the first 20:")   
		 print(SCT_dataset_GenesOnCHRwithoutSNPs[1:20])   	
		       	
	# Ensure that cell_dists / all_scts are data.matrix formatted (otherwise the gene names get lost)
	if(class(celltype_data[[annot_level]]$cell_dists)=="data.frame"){
		celltype_data[[annot_level]]$cell_dists = data.matrix(celltype_data[[annot_level]]$cell_dists)
	}	 
	if(class(celltype_data[[annot_level]]$all_scts)=="data.frame"){
		celltype_data[[annot_level]]$all_scts = data.matrix(celltype_data[[annot_level]]$all_scts)
	}	 

    # For each celltype
    # Get the genes in the set
    # target_celltype = colnames(celltype_data[[annot_level]]$cell_dists)[1]
    celltype_count=0
    print("All celltypes: ")
    print(colnames(celltype_data[[annot_level]]$cell_dists))
    for(target_celltype in colnames(celltype_data[[annot_level]]$cell_dists)){
    	print(target_celltype)
   	    celltype_count=celltype_count+1
    	template.file = template.file.O
		for(jjj in 1:5){print("")}
		for(jjj in 1:5){print(target_celltype)}
        cell_prop = celltype_data[[annot_level]]$cell_dists[,target_celltype]

        # CATCH POTENTIAL ERROR: PROGRAM CRASHES IF ANY OF CELL_PROP ARE NA... SO JUST DROP NA GENES... THEY RESULT FROM DIVIDE BY ZERO
        numNA = sum(is.na(cell_prop))
        if(numNA>0){
        	print(sprintf("REMOVING %s GENES WHICH HAVE SPECIFICITY METRICS: NA",numNA))
        	cell_prop = cell_prop[!is.na(cell_prop)]
        }
        
        # Check for potential error
        if(is.null(names(cell_prop))){stop("ERROR: names(cell_prop) is null. It should be a list of the gene names.")}

        ctd_file = gsub("\\(|\\)","",proportion_file)
        new_gzName = sprintf("%s_AnnotLevel_%s_%s_CutOff%s_%s_StepSize_%s",ctd_file,annot_level,gsub(" ","",target_celltype),opt$cutoff_thresh,opt$split_mode,opt$step_size)
        new_path = sprintf("%s/%s",ctd_file,new_gzName)
        

        # Create new results folder
        if(!file.exists(sprintf("%s/%s",opt$cell_annot_folder,ctd_file))){
            dir.create(sprintf("%s/%s",opt$cell_annot_folder,ctd_file))
        }
        if(!file.exists(sprintf("%s/%s",opt$cell_annot_folder,new_path))){
            dir.create(sprintf("%s/%s",opt$cell_annot_folder,new_path))
        }
        
        # Print quantiles
        #print("Quantiles:")
        #print(quantile(cell_prop2,probs=seq(from=0,to=0.9,by=0.1)))
        #sort(cell_prop[cell_prop>0.317])
        
        # Get the 
        # function(target_celltype,cell_prop,mouse_to_human_homologs,allKeptSNPS,template.file,template.snps,opt$split_mode)
        annotDescFile = sprintf("%s/%s/%s",opt$cell_annot_folder,new_path,"genePercentiles.csv")
        
        annotF = sprintf("%s/%s/%s",opt$cell_annot_folder,new_path,gsub("cell_type_group\\.1",new_gzName,template_files[tI]))

        # Generate the new annotation file ONLY if it doesn not already exist
        #if(!file.exists(annotF)){
	        #target_celltype,cell_prop,mouse_to_human_homologs,allKeptSNPS,template.file,template.snps,opt,annotDescFile

	        tmpFNAME = sprintf("%s/%.0f.rda",opt$tmp_folder,floor(runif(1)*100000000000000000))
	        print(tmpFNAME)
	        save.image(file=tmpFNAME)
	        cellAdd = add.cell.to.template(target_celltype,cell_prop,mouse_to_human_homologs,allKeptSNPS,template.file,template.snps,opt,annotDescFile)
	        template.file = cellAdd$template.file
	   		if(celltype_count==1){
	   			allGeneGroups=cellAdd$geneGroupTable
	   		}else{
	   			allGeneGroups=rbind(allGeneGroups,cellAdd$geneGroupTable)
	   		}
	   
	   		if(is.null(dim(template.file))){
	   			stop("For some reason the template file generation has failed")	
	   		}
	   
	        # Save the annot
	        #annotF = sprintf("%s/%s/%s",folder_path,annot_name,gsub("CNS",annot_name,gsub(".gz","",template_files[tI])))
	        
	        #annotF = sprintf("%s/%s/%s",folder_path,new_path,gsub("CNS",new_path,template_files[tI]))
	        
	        gz1 <- gzfile(annotF, "w")
	        write.table(template.file,sep="\t",quote=FALSE,row.names=FALSE,gz1)
	        close(gz1)
    	#}else{
    	#	warning(sprintf("Not generating annotation for %s as the file already exists",target_celltype))
    	#}
    }
    if(!file.exists(sprintf("%s/%s",opt$gene_groups_folder,opt$input_file))){
        dir.create(sprintf("%s/%s",opt$gene_groups_folder,opt$input_file))
    }    
    write.csv(allGeneGroups,file=sprintf("%s/%s/GeneGroups_%s_CutOff%s_%s_StepSize_%s.csv",opt$gene_groups_folder,opt$input_file,opt$input_file,opt$cutoff_thresh,opt$split_mode,opt$step_size))
}    

warnings()
