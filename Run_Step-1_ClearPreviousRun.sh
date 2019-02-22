#!/bin/bash
source $LDSCCODE/InitVars.sh

cd $LDSCSRC
for LOGFOLDER in Logs_Step1_MakeAnnots Logs_Step1b_GenerateGORes Logs_Step2_PrepAnnots Logs_Step4_PartitionHeritability Results
do
	if [ -d $LOGFOLDER ]; then
		rm -r $LOGFOLDER
	fi
done
mkdir $LDSCSRC/Results
mkdir $LDSCSRC/CellAnnots
#rm $SNPANNOT/GeneGroups_cell*
#rm $SNPANNOT/genegroups.txt