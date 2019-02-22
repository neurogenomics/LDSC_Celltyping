#!/usr/bin/env Rscript
# file_name="celltype_data_allKImouse_level1_thresh0_trim0"
args <- commandArgs(trailingOnly = TRUE)
file_name <- args[1]
path=Sys.getenv("CTDFOLDER") #"/lustre/scratch117/unknown/team141/ns9/LDSC/"
load(sprintf("%s/%s.rda",path,file_name))
cat(paste(gsub(" ","",colnames(celltype_data[[1]]$cell_dists)),collapse=" "))