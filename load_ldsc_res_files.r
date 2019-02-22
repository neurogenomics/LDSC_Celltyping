#' load_ldsc_res_files
#'
#' \code{load_ldsc_res_files} Loads in the output of the LDSC software used to partition heritability.
#' All the .results files created by the LDSC software should be placed into a 'ResFiles' directory within the 'resPath' folder.
#' This function expects that the filenames of the .results files indicates the GWAS they pertain to ($GWAS), and the celltype_data file
#' used to generate them ($FILE), as well as the celltype they are associated with ($CELLTYPE). The .results files should be named as follows:
#'
#'  '$GWAS_$FILE_annot1.$CELLTYPE.percentile.Thresh0.results'
#'
#' @param resPath Root folder for the LDSC results. Assumes that this folder contains a subfolder called 'ResFiles'. 
#' @param gwas_set Array of names of the GWAS files (which then form part of the .results filenames)
#' @return res Expression matrix with gene names as rownames.
#' @examples
#' exp2 = drop.uninformative.genes(exp=rpkm,level2annot=all_annot$Tasic_et_al_2016_label)
#' @export
load_ldsc_res_files <- function(resPath,gwas_set,resFolder="ResFiles",validFileGrep="Thresh0.results",stepSize=0.1){
    library("stringr")
    # Check that 'ResFiles' folder exists
    if(!file.exists(sprintf("%s/%s",resPath,resFolder))){
        stop(sprintf("ERROR: the .results files should be stored in a folder named %s within the folder referenced within the resPath argument",resFolder))
    }
    
    #all_res_files = list.files(path=sprintf("%s/%s",resPath,resFolder))
    all_res_files = list.files(path=sprintf("%s/%s",resPath,resFolder),recursive=TRUE)
    all_res_files = all_res_files[grep(".log",all_res_files,invert=TRUE)]
    all_res_files = all_res_files[grep(validFileGrep,all_res_files,invert=FALSE)]
    
    # Get all the GWAS for which there are results
    #all_gwas_present = list.files(sprintf("%s/Results",Sys.getenv("LDSCSRC")))
    all_gwas_present = list.files(sprintf("%s/%s",resPath,resFolder))
    # Check if all present GWAS are in GWAS_SET
    gwas_set=trimws(gwas_set)
    allPresent=TRUE
    if(sum(!all_gwas_present %in% gwas_set)>0){allPresent=FALSE}
    missing_gwas = all_gwas_present[!all_gwas_present %in% gwas_set]
    if(allPresent==FALSE){
        stop(sprintf("ERROR: the following GWAS studies are missing from the gwas_set table: %s",missing_gwas))
    }
    #gwas_set[gwas_set %in% all_gwas_present,]

    # Check that there are results files for each GWAS
    real_gwas_present = rep(0,length(all_gwas_present))
    count=0
    for(gwas in all_gwas_present){
        count=count+1
        if(length(grep(gwas,all_res_files))>0){ # If there is at least one results file for this GWAS 
            real_gwas_present[count]=1
        }
    }    
    all_gwas_present = all_gwas_present[which(real_gwas_present>0)]

    for(gwas in all_gwas_present){
        print(gwas)
        if(length(grep(gwas,all_res_files))>0){
            all_gwas_res_files = all_res_files[grep(gwas,all_res_files)]
            if(resFolder=="ResFiles"){
                ### IF LOADING RESULTS FOR CELL TYPES                
                gwas_res_files = gsub(sprintf("%s_",gwas),"",all_gwas_res_files)
                #files = unique(gsub("_annot.*","",gwas_res_files))
                files=list.files(sprintf("%s/Results/%s",Sys.getenv("LDSCSRC"),gwas))
                count=0
                for(file in files){
                    if(length(grep(file,gwas_res_files))>0){ # If there is at least one results file for this GWAS & CTD FILE
                        print(file)

                        #all_file_res_files = gwas_res_files[grep(file,gwas_res_files)]
                        #file_res_files = gsub(sprintf("%s_annot1.",file),"",all_file_res_files)
                        #file_res_files2 = gsub(sprintf("\\.percentile\\.Thresh%s\\.results",thresh),"",file_res_files)
                        #file_res_files2 = gsub(sprintf("\\.percentile\\.%s",validFileGrep),"",file_res_files)
                        #celltypes = unique(file_res_files2)

                        # Extract names of all celltypes
                        relFiles = list.files(sprintf("%s/Results/%s/%s",Sys.getenv("LDSCSRC"),gwas,file),pattern=".results")
                        celltypes = unlist(str_match(relFiles, "annot\\d.(.*).percentile")[,2])

                        for(ct in celltypes){
                            print(ct)
                            filename = sprintf("%s/ResFiles/%s/%s/%s_%s_annot1.%s.%s.%s",folder,gwas,file,gwas,file,ct,split_mode,validFileGrep)
                            data = read.csv(filename,stringsAsFactors=FALSE,header=TRUE,sep="\t")
                            percentile=c(NA,seq(from=0,to=1,by=stepSize))
                            temp = data.frame(gwas=gwas,file=file,thresh=thresh,celltype=ct,split=split_mode,percentile=percentile,data[1:length(percentile),])#data[1:12,])
                            if(count==0){
                                allRes = temp
                            }else{
                                allRes = rbind(allRes,temp)
                            }
                            count=count+1
                        }
                    }
                }
            }else{
                ### IF LOADING RESULTS FOR GENE SETS
                count=0
                for(file in all_gwas_res_files){
                    print(file)
                    filename = sprintf("%s/%s/%s",resPath,resFolder,file)
                    data = read.csv(filename,stringsAsFactors=FALSE,header=TRUE,sep="\t")
                    data = cbind(data,gwas=gwas)
                    if(count==0){
                        allRes = data
                    }else{
                        allRes = rbind(allRes,data)
                    }
                    count=count+1                    
                }
            }
            if(gwas==all_gwas_present[1]){
                all_gwas = allRes
            }else{
                all_gwas = rbind(all_gwas,allRes)
            }
        }else{
            print(sprintf("No results files for: %s",gwas))
            gwas_set = gwas_set[!gwas_set %in% gwas]
        }
    }
    
    if("files" %in% colnames(all_gwas) ){
        # CORRECT FOR MULTIPLE TESTING
        all_gwas$Enrichment_q = 0
        for(file in all_gwas$files){
            p = all_gwas[all_gwas$file==file,]$Enrichment_p
            all_gwas[all_gwas$file==file,]$Enrichment_q = p.adjust(p,method="BH")
        }
    }
    #all_gwas$Enrichment_q = p.adjust(all_gwas$Enrichment_p,method="BH")
    return(all_gwas)
}