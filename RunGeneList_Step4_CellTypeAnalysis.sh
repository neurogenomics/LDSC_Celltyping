#!/bin/sh
cd $LDSCSRC

export TEMP=Logs_GeneList_Step4_PartitionHeritability
mkdir $TEMP

#export ANNOTPATH="/nfs/team156/ns9/LDSC/snp_annot_files/GeneSets/"

export NUMPARTITIONS=20

for (( PARTITION=1; PARTITION<=$NUMPARTITIONS; PARTITION=PARTITION+1 ))
do
	export ANNOTPATH='/nfs/team156/ns9/LDSC/snp_annot_files/GeneSetAnnot.p'$PARTITION'n'$NUMPARTITIONS'/'
	for GWASFILE in IQ.Sniekers.2017.txt EUR.IBD.gwas_info03_filtered.assoc_edited GWAS_EA_excl23andMe.txt GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq_WITHSNPLOC.txt schiz.qjecp.clozukpgc2.gz Alzh_IGAP_stage_1 BIPOLAR_daner_PGC_BIP32b_mds7a.BIP25_noboma.txt scz2.snp.results.txt EUR.CD.gwas_info03_filtered.assoc EduYears_Main.txt GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt ADHD.pgc.adhd.full.2012-10.txt MDD.pgc.mdd.full.2012-04.txt SubjectiveWellBeing.SWB_Full.txt.gz DepressiveSymptoms.DS_Full.txt.gz SSGAC_EduYears_Rietveld2013_publicrelease.txt EUR.IBD.gwas_info03_filtered.assoc EUR.UC.gwas_info03_filtered.assoc Eating.Disorders.pgc.an.snp.all.13May2016.txt 
	do			
		export QUEUE=normal
	
		export TAGNAME=GeneLists_${GWASFILE}_${PARTITION}_${NUMPARTITIONS}
		
		export CMDSTR="/nfs/team156/ns9/LDSC/env/bin/python ldsc.py --h2 /nfs/team156/ns9/LDSC/GWAS_Processed_Stats/${GWASFILE}.sumstats.gz\
	    --w-ld-chr weights_hm3_no_hla/weights.\
	    --ref-ld-chr ${ANNOTPATH}GeneSetAnnot.,baseline_1000G_Phase3/baseline.\
	    --overlap-annot\
	    --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC.\
	    --out /nfs/team156/ns9/LDSC/Results/${TAGNAME}\
	    --print-coefficients"
	
		bsub -q $QUEUE -R 'select[type==X86_64 && mem>16000] rusage[mem=16000]' -M 16000 -o $TEMP'/'${TAGNAME}'.output' -e $TEMP'/'${TAGNAME}'.error'  $CMDSTR
	done
done