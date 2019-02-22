#!/bin/sh
cd /nfs/team156/ns9/LDSC

export TEMP=Logs_GeneList_Step2_PrepAnnots
mkdir $TEMP

export NUMPARTITIONS=20

for (( PARTITION=1; PARTITION<=$NUMPARTITIONS; PARTITION=PARTITION+1 ))
do
	export ANNOTPATH='/nfs/team156/ns9/LDSC/snp_annot_files/GeneSetAnnot.p'$PARTITION'n'$NUMPARTITIONS
	for (( i=1; i<=22; i=i+1 ))
	do		
		export QUEUE=normal
	
		export CMDSTR="/nfs/team156/ns9/LDSC/env/bin/python ldsc.py\
		    --l2\
		    --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.$i\
		    --ld-wind-cm 1\
		    --annot ${ANNOTPATH}/GeneSetAnnot.$i.annot.gz\
		    --out ${ANNOTPATH}/GeneSetAnnot.$i\
		    --print-snps hapmap3_snps/hm.$i.snp"
		
		export LOGOUT=$TEMP'/test_chr'$i'_p'$PARTITION'_n'$NUMPARTITIONS'.output'
		export LOGERR=$TEMP'/test_chr'$i'_p'$PARTITION'_n'$NUMPARTITIONS'.error'					
		
		if [ -f "${ANNOTPATH}/GeneSetAnnot.$i.annot.gz" ];
		then
		    bsub -n2 -q $QUEUE -R 'select[type==X86_64 && mem>16000] rusage[mem=16000] span[hosts=1]' -M 16000 -o $LOGOUT -e $LOGERR $CMDSTR
		    #gunzip ${ANNOTPATH}/GeneSetAnnot.$i.annot.gz 
		    #sed -e '1s/Specific to Both/SpecificToBoth/' ${ANNOTPATH}/GeneSetAnnot.$i.annot
		    #sed -e '1s/High Expressed/HighExpressed/' ${ANNOTPATH}/GeneSetAnnot.$i.annot		    
		    #sed -e '1s/Low Expressed/LowExpressed/' ${ANNOTPATH}/GeneSetAnnot.$i.annot		    		    
		    #gzip ${ANNOTPATH}/GeneSetAnnot.$i.annot
		    #echo ${ANNOTPATH}/GeneSetAnnot.$i.annot
		else
		    echo ${ANNOTPATH}/GeneSetAnnot.$i.annot.gz;
				    echo "File does not exist" 
		fi
	done
done