#!/bin/sh
cd /nfs/team156/ns9/LDSC
export TEMP=Logs_GeneList_Step1_MakeAnnots
mkdir $TEMP

export QUEUE=normal
export NUMPARTITIONS=20

for (( PARTITION=1; PARTITION<=$NUMPARTITIONS; PARTITION=PARTITION+1 ))
do
	for (( TI=1; TI<=22; TI=TI+1 ))
	do
		export LOGNAME='log_p_'$PARTITION'_n_'$NUMPARTITIONS'_chr_'$TI
		#echo $LOGNAME
		bsub -q $QUEUE -R 'select[type==X86_64 && mem>32000] rusage[mem=32000]' -M 32000 -o $TEMP'/'$LOGNAME'.output' -e $TEMP'/'$LOGNAME'.error' 'R-3.0.0 --file=create_LDSC_annot_fromGeneSets.r --args --chr '$TI' --num_partitions '$NUMPARTITIONS' --partition '$PARTITION
	done
done
		

