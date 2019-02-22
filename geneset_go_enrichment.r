library(topGO)
library("biomaRt")

# Get command line arguments
library("optparse")

option_list = list(
    make_option(c("-i", "--listi"), type="integer", default=1,
                help="Index value of template annotation file to be loaded (normally 1--22)")
)
opt = parse_args(OptionParser(option_list=option_list))


#human = useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
#mouse = useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
#mgi_symbols = getBM(attributes=c("mgi_symbol"), mart=mouse)
#hgnc_symbols = getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), mart=human)
#save(hgnc_symbols,file="hgnc_symbols.rda")
#mouse_to_human_homologs = getLDS(attributes = c("mgi_symbol","entrezgene","ensembl_gene_id"),
#                                 filters = "mgi_symbol", values = mgi_symbols[,1],
#                                 mart = mouse,
#                                 attributesL = c("hgnc_symbol","ensembl_gene_id","entrezgene"), martL = human)
load("/lustre/scratch117/unknown/team141/ns9/LDSC/mouse_to_human_homologs.rda")

all_ensembl = unique(mouse_to_human_homologs$Ensembl.Gene.ID.1)
all_ensembl = all_ensembl[all_ensembl!=""]

all_genes <- rep(0,length(all_ensembl))
names(all_genes) <- unique(all_ensembl)

#compare_folder = "/Users/ns9/Desktop/LDSC Results/Phase3bRes/CellTypeComparison_text_files"
compare_folder = "/lustre/scratch117/unknown/team141/ns9/LDSC/GeneSets/Human"
all_files = list.files(compare_folder,pattern="*.txt")
#ff = "EmbryonicDopaminergicNeuron0_allKImouse_level2_Specific_VS_EmbryonicDopaminergicNeuron0_HumanMidbrainKI_level2_Specific to Both.txt"
#ff = "EmbryonicDopaminergicNeuron0_HumanMidbrainKI_level2_Specific_VS_EmbryonicDopaminergicNeuron0_allKImouse_level2_Low Expressed.txt"
ff  = all_files[opt$listi]
ff_dat = read.csv(sprintf("%s/%s",compare_folder,ff),stringsAsFactors = FALSE)[-c(1),1]
ff_dat_ensembl = unique(mouse_to_human_homologs[mouse_to_human_homologs$HGNC.symbol %in% ff_dat,]$Ensembl.Gene.ID.1)

all_genes[ff_dat_ensembl] = 1

GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = all_genes, geneSel = function(p) p > 0.9, description = "Test", annot = annFUN.org, mapping="org.Hs.eg.db", ID="Ensembl")
GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = all_genes, geneSel = function(p) p > 0.9, description = "Test", annot = annFUN.org, mapping="org.Hs.eg.db", ID="Ensembl")
GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = all_genes, geneSel = function(p) p > 0.9, description = "Test", annot = annFUN.org, mapping="org.Hs.eg.db", ID="Ensembl")

#resultFisher_BP <- runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")
#resultFisher_CC <- runTest(GOdata_CC, algorithm = "classic", statistic = "fisher")
#resultFisher_MF <- runTest(GOdata_MF, algorithm = "classic", statistic = "fisher")
resultFisher_BP <- runTest(GOdata_BP, algorithm = "weight01", statistic = "fisher")
resultFisher_CC <- runTest(GOdata_CC, algorithm = "weight01", statistic = "fisher")
resultFisher_MF <- runTest(GOdata_MF, algorithm = "weight01", statistic = "fisher")


a=cbind(GenTable(GOdata_BP, classicFisher = resultFisher_BP, topNodes = 30),ontology="BP")
b=cbind(GenTable(GOdata_CC, classicFisher = resultFisher_CC, topNodes = 30),ontology="CC")
c=cbind(GenTable(GOdata_MF, classicFisher = resultFisher_MF, topNodes = 30),ontology="MF")
d=rbind(a,b,c)
d=d[order(as.numeric(d$classicFisher)),]
d$genes=""
#d = d[as.numeric(d$classicFisher)<0.05 & d$Significant>2,]
print(dim(d))
d = d[as.numeric(d$Significant>2),]

# Get genes associated with each term
for(i in 1:dim(d)[1]){
    if(as.character(d[i,"ontology"])=="BP"){  ann.genes <- genesInTerm(GOdata_BP, d[i,"GO.ID"]) }
    if(as.character(d[i,"ontology"])=="CC"){  ann.genes <- genesInTerm(GOdata_CC, d[i,"GO.ID"]) }
    if(as.character(d[i,"ontology"])=="MF"){  ann.genes <- genesInTerm(GOdata_MF, d[i,"GO.ID"]) }
    #ff_dat_ensembl = unique(mouse_to_human_homologs[mouse_to_human_homologs$HGNC.symbol %in% ff_dat,]$Ensembl.Gene.ID.1)
    ens_in_go = ff_dat_ensembl[ff_dat_ensembl %in% ann.genes[[1]]]
    go_term_genes = mouse_to_human_homologs[mouse_to_human_homologs$Ensembl.Gene.ID.1 %in% ff_dat_ensembl[ff_dat_ensembl %in% ann.genes[[1]]],"HGNC.symbol"]
    d[i,]$genes = paste(go_term_genes,collapse=",")
}

if(dim(d)[1]>2){
	isDup = rep(FALSE,dim(d)[1])
	for(i in 2:dim(d)[1]){
	    for(j in 1:(i-1)){
	        list1=strsplit(d[i,]$genes,split=",")[[1]]
	        list2=strsplit(d[j,]$genes,split=",")[[1]]
	        if(sum(list1 %in% list2)/length(list1)>0.5){
	            isDup[i] = TRUE
	            print(sprintf("%s %s",i,j))
	        }
	    }    
	}
	e=d[!isDup,]
}else{
    e=d	
}

res_path = "/lustre/scratch117/unknown/team141/ns9/LDSC/GeneSets/Human_GOterms"
ff2 = gsub(".txt",".GOterms.csv",ff)
write.csv(e,file=sprintf("%s/%s",res_path,ff2))