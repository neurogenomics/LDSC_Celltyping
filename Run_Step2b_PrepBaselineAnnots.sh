#!/bin/sh
cd $LDSCSRC
source $LDSCCODE/InitVars.sh

export TEMP=Logs_Step2b_PrepBaselineAnnots
mkdir $TEMP

#export ANNOTPATH="/lustre/scratch117/unknown/team141/ns9/LDSC/baseline_1000G_Phase3/"
export ANNOTPATH="$HOME/LDSC/1000G_EUR_Phase3_baseline/"

for (( i=1; i<=22; i=i+1 ))
do		
	export QUEUE=normal

	export CMDSTR="python $HOME/LDSC/ldsc.py\
	    --l2\
	    --bfile $HOME/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i\
	    --ld-wind-cm 1\
	    --annot ${ANNOTPATH}/baseline.$i.annot.gz\
	    --out ${ANNOTPATH}/baseline.$i\
	    --print-snps $HOME/LDSC/hapmap3_snps/hm.$i.snp"
	
	export LOGOUT=$TEMP'/test_'$i'.output'
	export LOGERR=$TEMP'/test_'$i'.error'					
	
	#bsub -n2 -q $QUEUE -R 'select[type==X86_64 && mem>8000] rusage[mem=8000] span[hosts=1]' -M 8000 -o $LOGOUT -e $LOGERR $CMDSTR
	sbatch -t 120 -n 2 --mem=8000 -o $LOGOUT -e $LOGERR --wrap="${CMDSTR}"	
done