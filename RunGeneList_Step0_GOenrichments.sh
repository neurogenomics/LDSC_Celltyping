#!/bin/sh
cd /nfs/team156/ns9/LDSC
export TEMP=Logs_GeneList_Step0_GOenrichments
mkdir $TEMP

export QUEUE=normal
export NUMLISTS=1539
#export NUMLISTS=10

for (( LIST=1; LIST<=$NUMLISTS; LIST=LIST+1 ))
do
	export LOGNAME='log_list_'$LIST'_numLists_'$NUMLISTS
	#echo $LOGNAME
	bsub -q $QUEUE -R 'select[type==X86_64 && mem>16000] rusage[mem=16000]' -M 16000 -o $TEMP'/'$LOGNAME'.output' -e $TEMP'/'$LOGNAME'.error' 'R-3.0.0 --file=geneset_go_enrichment.r --args --listi '$LIST
done
		

