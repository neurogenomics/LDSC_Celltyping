#!/bin/bash
# This file sets the paths to all folders
# It is intended that any changes required to swap the code onto a new SLURM based server can be set here
# It should be loaded at the start of all other scripts

alias pythonenv="/nas/longleaf/home/nskene/env/bin/python"

export LDSCSRC=$SCRATCH/LDSC
if [ ! -d $LDSCSRC ]; then
  mkdir $LDSCSRC
fi

export SNPANNOT=$HOME/LDSC/snp_annot_files
#echo "LDSC_Celltyping variables set"

# FOLDER WHERE CODE (i.e. this file) IS LOCATED
export LDSCCODE=$HOME/LDSC_Celltyping
export LDSCPROG=$HOME/LDSC

# FOLDER WHERE CELLTYPE_DATA/CTD FILES ARE LOCATED
export CTDFOLDER=$SCRATCH/Celltype_data

# FOLDER WHERE EDITED SUMSTATS ARE LOCATED
export SUMSTATS_FOLDER=$SCRATCH/GWAS_SumStats/Edited
export MUNGED_FOLDER=$SCRATCH/GWAS_SumStats/Munged


# CHECK THAT SNP_MatchingGenes files have been prepped
#echo "Testing if matchingGenes files have been prepped..."
if [ ! -f "$SNPANNOT/1000G_Phase3_cell_type_groups/cell_type_group.1.1.matchingGenes.rda" ]
then
	#echo "Files need prepping"
	unzip SNP_MatchingGenes.zip
	mv SNP_MatchingGenes/* $SNPANNOT/1000G_Phase3_cell_type_groups/
	rm -r SNP_MatchingGenes
fi

#module load r/3.4.1
module load python/2.7.12