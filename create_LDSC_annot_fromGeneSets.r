get.mouse_to_human_homologs <- function(opt){
	load("mouse_to_human_homologs.rda")
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


# Get command line arguments
library("optparse")

option_list = list(
  make_option(c("-c", "--chr"), type="integer", default=1,
              help="Index value of template annotation file to be loaded (normally 1--22)"),
  make_option(c("-N", "--num_partitions"), type="integer", default=20,
              help="Number of partitions to divide the gene sets between"),
  make_option(c("-P", "--partition"), type="integer", default="1", 
              help="Which partition should this instance work on")
)
opt = parse_args(OptionParser(option_list=option_list))

print(sprintf("Chr: %s",opt$chr))
print(sprintf("Partition: %s / %s",opt$partition,opt$num_partitions))

## SETUP FOLDER VARIABLES
folder_path   = "/lustre/scratch117/unknown/team141/ns9/LDSC/snp_annot_files"
template_path = "1000G_Phase3_cell_type_groups"
template_files = list.files(sprintf("%s/%s",folder_path,template_path))
templateBaseRoot = "cell_type_group.1."
templateBaseName = paste(templateBaseRoot,"%s.annot.gz",sep="")
template_files = template_files[grep("cell_type_group\\.1\\..*\\.annot.gz",template_files)]

## LOAD MOUSE->HUMAN HOMOLOGS
mouse_to_human_homologs = get.mouse_to_human_homologs(opt)

## LOAD GENE LISTS
print("Reading in gene lists...")
allGeneLists = list()
mouse.gene.lists = list.files(path="/lustre/scratch117/unknown/team141/ns9/LDSC/GeneSets/Mouse")
for(ms in mouse.gene.lists){
	print(ms)
	genelist = unlist(read.csv(sprintf("/lustre/scratch117/unknown/team141/ns9/LDSC/GeneSets/Mouse/%s",ms),stringsAsFactors=FALSE,sep="!")[,1][-(1:2)])
	human.genelist = unique(mouse_to_human_homologs[mouse_to_human_homologs$MGI.symbol %in% genelist,"HGNC.symbol"])
	allGeneLists[[length(allGeneLists)+1]] = list(name=ms,genelist=human.genelist)
}
human.gene.lists = list.files(path="/lustre/scratch117/unknown/team141/ns9/LDSC/GeneSets/Human")	
for(hs in human.gene.lists){
	print(hs)
	human.genelist = unlist(read.csv(sprintf("/lustre/scratch117/unknown/team141/ns9/LDSC/GeneSets/Human/%s",hs),stringsAsFactors=FALSE)[,1][-(1:2)])
	allGeneLists[[length(allGeneLists)+1]] = list(name=hs,genelist=human.genelist)
}

## LOOP OVER THE TEMPLATE FILES
#for(chr in 1:length(template_files)){
print(sprintf("Chr: %s",opt$chr))
# MAP CHR to TEMPLATE FILE
tI = which(template_files==sprintf(templateBaseName,opt$chr))

## LOAD THE SNP->GENE FILE
print("Loading the SNP->gene file")
templateName = sprintf("%s/%s/%s",folder_path,template_path,template_files[tI])
newFName = gsub("annot.gz","matchingGenes.rda",templateName)
print(newFName)
# - First, check that it actually exists 
#    * if you've switched to a new set of SNP annot template files, it will need regenerating
if(!file.exists(newFName)){	stop("Prepare Template mode needs to be run first. Use -PR argument.")   }
load(newFName)	

## LOAD THE TEMPLATE FILE
print("Loading the template file")
template.file <- read.csv(gzfile(sprintf("%s/%s/%s",folder_path,template_path,template_files[tI])), as.is = TRUE,sep="\t")
template.snps <- as.numeric(gsub("rs|ss","",template.file$SNP))
template.file = template.file[,-5]
template.file.O = template.file 

## PARTITION THE GENE LISTS
partition_size    = ceiling(length(allGeneLists)/opt$num_partitions)
partition_start   = (opt$partition-1)*partition_size+1
partition_end     = opt$partition*partition_size
if(partition_end>length(allGeneLists)){partition_end=length(allGeneLists)}
partition_indices = partition_start:partition_end

## PAD THE TEMPLATE FILE
template.file = cbind(template.file,matrix(0,nrow=dim(template.file)[1],ncol=length(partition_indices)))

## ADD COLUMN TO TEMPLATE FILE
count=0
for(i in partition_indices){
	count=count+1
	print(sprintf("Adding genelist[%s] to the template",i))
	print(allGeneLists[[i]]$name)
	# KEEP ONLY GENES ON THE RIGHT CHROMOSOME
	keptGenes = as.character(unique(allKeptSNPS$V2)[unique(allKeptSNPS$V2) %in% allGeneLists[[i]]$genelist])	    	
	
	# Find SNPs associated with genes in the target set
    #numKeptGenes = length(keptGenes)
    #countKeptGenes = countKeptGenes + numKeptGenes
    target.out = unique(allKeptSNPS[allKeptSNPS$V2 %in% keptGenes,])
    
    # Print stats about genes and SNPs
    target.snps = unique(allKeptSNPS$V1[allKeptSNPS$V2 %in% keptGenes])
    inAnnot = rep(0,dim(template.file)[1])
    inAnnot[template.snps %in% target.snps] = 1
    # If the column is empty, set 50 random SNPs to 1 to avoid numerical errors
    if(sum(inAnnot)==0){inAnnot[sample(1:length(inAnnot),50)] = 1} 
    
    template.file[,4+count] = inAnnot
    # Drop spaces from the colname
    cName = gsub(" ","",allGeneLists[[i]]$name)
    colnames(template.file)[4+count] = gsub("\\.txt","",cName)
	print(sprintf("Finished adding genelist[%s] to the template",i))	    
}

tag=sprintf("GeneSetAnnot.p%sn%s",opt$partition,opt$num_partitions)
# Create new results folder
if(!file.exists(sprintf("%s/%s",folder_path,tag))){
    dir.create(sprintf("%s/%s",folder_path,tag))
}
annotF = sprintf("%s/%s/%s",folder_path,tag,gsub("cell_type_group\\.1","GeneSetAnnot",template_files[tI]))
gz1 <- gzfile(annotF, "w")
write.table(template.file,sep="\t",quote=FALSE,row.names=FALSE,gz1)
close(gz1)	

warnings()
