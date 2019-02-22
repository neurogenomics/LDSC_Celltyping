#!/bin/sh
cd $LDSCSRC
source $LDSCCODE/InitVars.sh

export TEMP=Logs_Step3a_MungeSumStats
mkdir $TEMP

export QUEUE=normal

unset PYTHONPATH

#for INPUTFILE in ParkinsonsSumStats2018_23andMe.tab # Education_trait_1.txt HI.txt # BIPOLAR_daner_PGC_BIP32b_mds7a.BIP25_noboma.txt EduYears_Main.txt.gz Alzh_IGAP_stage_1.txt # BroadABC_METALoutput__Males.tbl BroadABC_METALoutput_Females.tbl BroadABC_METALoutput_Combined.tbl scz2.snp.results.txt Schiz_Clozuk.txt # BMI.SNPadjSMK.CombinedSexes.EuropeanOnly.txt BMI.SNPadjSMK.Men.EuropeanOnly.txt BMI.SNPadjSMK.Women.EuropeanOnly.txt # IQ.Sneikers2016.sumstats.txt GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz  #EAGLE_20151029_n10k_GC1_Psort_WITHSNPLOC.tbl als.2016.txt iq_MOD.txt  # Alzh_IGAP_stage_1.txt scz2.snp.results.txt Eating.Disorders.pgc.an.snp.all.13May2016.txt GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt ADHD.pgc.adhd.full.2012-10.txt MDD.pgc.mdd.full.2012-04.txt EduYears_Main.txt SubjectiveWellBeing.SWB_Full.txt.gz DepressiveSymptoms.DS_Full.txt.gz SSGAC_EduYears_Rietveld2013_publicrelease.txt EUR.CD.gwas_info03_filtered.assoc EUR.IBD.gwas_info03_filtered.assoc EUR.UC.gwas_info03_filtered.assoc schiz.qjecp.clozukpgc2.gz parkinsons.phs000089.pha002868.clean.txt
for INPUTFILE in IQ.Sniekers.2017.txt EUR.IBD.gwas_info03_filtered.assoc_edited GWAS_EA_excl23andMe.txt GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq_WITHSNPLOC.txt ParkinsonsSumStats2018_23andMe.tab 2040_RiskTaking_UKBB.assoc.tsv 2030_GuiltyFeelings_UKBB.assoc.tsv iPSYCH-PGC_ASD_Nov2017 CLOZUK_PGC2noclo.METAL.assoc.dosage.fix Alzh_IGAP_stage_1.txt
do
	echo $INPUTFILE
	export IGNORE=""
	case "$INPUTFILE" in
		"ParkinsonsSumStats2018_23andMe.tab")
			export NNN=-1
		;;
		"2040_RiskTaking_UKBB.assoc.tsv")
			export NNN=-1
		;;
		"2030_GuiltyFeelings_UKBB.assoc.tsv")
			export NNN=-1
		;;
		"iPSYCH-PGC_ASD_Nov2017")
			export NNN=46351
		;;
		"CLOZUK_PGC2noclo.METAL.assoc.dosage.fix")
			export NNN=35802
		;;
		"Education_trait_1.txt")
			export NNN=-1
		;;
		"HI.txt")
			export NNN=-1
		;;
		"Alzh_IGAP_stage_1.txt")
			export NNN=54162
		;;
		"Alzh_IGAP_stage_1_2_combined.txt")
			export NNN=74046
		;;
		"Eating.Disorders.pgc.an.snp.all.13May2016.txt")
			export NNN=14477
		;; 
		"GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz") 
			export NNN=253288
		;;
		"GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq_WITHSNPLOC.txt") 
			export NNN=253288
		;;
		"ADHD.pgc.adhd.full.2012-10.txt MDD.pgc.mdd.full.2012-04.txt") 
			export NNN=5415
		;;
		"EduYears_Main.txt.gz")
			export NNN=300000
		;; 
		"SubjectiveWellBeing.SWB_Full.txt.gz") 
			export NNN=298420
		;;
		"DepressiveSymptoms.DS_Full.txt.gz") 
			export NNN=161460
		;;
		"SSGAC_EduYears_Rietveld2013_publicrelease.txt") 
			export NNN=125000
		;;
		"EUR.CD.gwas_info03_filtered.assoc") 
			export NNN=86640
		;;
		"EUR.IBD.gwas_info03_filtered.assoc") 
			export NNN=34652		
		;;
		"EUR.IBD.gwas_info03_filtered.assoc_edited") 
			export NNN=34652		
		;;
		"EUR.UC.gwas_info03_filtered.assoc") 
			export NNN=86640		
		;;
		"Schiz_Clozuk.txt") 
			export NNN=35802
		;;
		"schiz.qjecp.clozukpgc2.gz") 
			export NNN=35802
		;;
		"parkinsons.phs000089.pha002868.clean.txt")
			export NNN=5691
		;;
		"scz2.snp.results.txt")
			export NNN=28918
		;;
		"BIPOLAR_daner_PGC_BIP32b_mds7a.BIP25_noboma.txt")
			export NNN=44572
		;;
		"EAGLE_20151029_n10k_GC1_Psort_WITHSNPLOC.tbl")
			export NNN=17666
		;;
		"als.2016.txt")
			export NNN=36052
		;;
		"iq_MOD.txt")
			export NNN=-1
		;;			
		"BMI.SNPadjSMK.CombinedSexes.EuropeanOnly.txt")
			export NNN=-1
			export IGNORE="markername"
		;;	
		"BMI.SNPadjSMK.Men.EuropeanOnly.txt")
			export NNN=-1
			export IGNORE="markername"
		;;	
		"BMI.SNPadjSMK.Women.EuropeanOnly.txt")
			export NNN=-1
			export IGNORE="markername"
		;;	
		"BroadABC_METALoutput_Combined.tbl")
			export NNN=16400
			export IGNORE="MarkerName"
		;;	
		"BroadABC_METALoutput_Females.tbl")
			export NNN=8535
			export IGNORE="MarkerName"
		;;	
		"BroadABC_METALoutput__Males.tbl")
			export NNN=7772
			export IGNORE="MarkerName"
		;;	
		"IQ.Sneikers2016.sumstats.txt")
			export NNN=78308
			export IGNORE="Beta"
		;;
		"IQ.Sniekers.2017.txt")
			export NNN=78308
			export IGNORE="Beta"
		;;
		"GWAS_EA_excl23andMe.txt")
			export NNN=766345
		;;
		
	esac	
	
	for GWASPATH in "$SCRATCH/GWAS_SumStats/Raw" "$SCRATCH/GWAS_SumStats/Edited"
	do
		if [ -f "$GWASPATH/${INPUTFILE}" ]
		then
			if [ $NNN == -1 ]
			then
				export CMDSTR="/nas/longleaf/home/nskene/anaconda_ete/envs/envLDSC/bin/python $HOME/LDSC/munge_sumstats.py \
					--out $SCRATCH/GWAS_SumStats/Munged/${INPUTFILE} \
					--merge-alleles $LDSCCODE/w_hm3.snplist \
					--sumstats $GWASPATH/${INPUTFILE}" 
			else
				export CMDSTR="/nas/longleaf/home/nskene/anaconda_ete/envs/envLDSC/bin/python $HOME/LDSC/munge_sumstats.py \
					--out $SCRATCH/GWAS_SumStats/Munged/${INPUTFILE} \
					--merge-alleles $LDSCCODE/w_hm3.snplist \
					--N ${NNN} \
					--sumstats $GWASPATH/${INPUTFILE}" 					
			fi
			echo $IGNORE
			if [ "$IGNORE" != "" ]
			then
				export CMDSTR="$CMDSTR --ignore $IGNORE"
			fi
			echo $CMDSTR
			export LOGOUT=$TEMP'/test_'$INPUTFILE'.output' 
			export LOGERR=$TEMP'/test_'$INPUTFILE'.error'
			#bsub -q $QUEUE -R 'select[type==X86_64 && mem>8000] rusage[mem=8000]' -M 8000 -o $TEMP'/test_'$INPUTFILE'.output' -e $TEMP'/test_'$INPUTFILE'.error'  $CMDSTR
			sbatch -t 120 -n 1 --mem=8000 -o $LOGOUT -e $LOGERR --wrap="${CMDSTR}"	
		fi
	done
done