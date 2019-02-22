#!/bin/sh
cd $LDSCSRC
source $LDSCCODE/InitVars.sh

export TEMP=Logs_Step4_PartitionHeritability
mkdir $TEMP

#export ANNOTPATH="/lustre/scratch117/unknown/team141/ns9/LDSC/snp_annot_files/"
export ANNOTPATH="$LDSCSRC/CellAnnots/"

export THRESH=0
export STEPSIZE=0.1
export SPLITMODE=percentile

unset PYTHONPATH

count=0
countAdded=0
needsPrep=0
alreadyExists=0
#for GWASFILE in #EAGLE_20151029_n10k_GC1_Psort_WITHSNPLOC.tbl als.2016.txt iq_MOD.txt scz2.snp.results.txt # schiz.qjecp.clozukpgc2.gz Alzh_IGAP_stage_1.txt BIPOLAR_daner_PGC_BIP32b_mds7a.BIP25_noboma.txt scz2.snp.results.txt EUR.CD.gwas_info03_filtered.assoc EduYears_Main.txt GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt ADHD.pgc.adhd.full.2012-10.txt MDD.pgc.mdd.full.2012-04.txt SubjectiveWellBeing.SWB_Full.txt.gz DepressiveSymptoms.DS_Full.txt.gz SSGAC_EduYears_Rietveld2013_publicrelease.txt EUR.IBD.gwas_info03_filtered.assoc EUR.UC.gwas_info03_filtered.assoc Eating.Disorders.pgc.an.snp.all.13May2016.txt
#for GWASFILE in Education_trait_1.txt HI.txt  # IQ.Sneikers2016.sumstats.txt # BIPOLAR_daner_PGC_BIP32b_mds7a.BIP25_noboma.txt Schiz_Clozuk.txt # Alzh_IGAP_stage_1.txt #EduYears_Main.txt.gz #  # BroadABC_METALoutput__Males.tbl BroadABC_METALoutput_Females.tbl BroadABC_METALoutput_Combined.tbl # scz2.snp.results.txt BMI.SNPadjSMK.CombinedSexes.EuropeanOnly.txt BMI.SNPadjSMK.Men.EuropeanOnly.txt BMI.SNPadjSMK.Women.EuropeanOnly.txt GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz  #  EAGLE_20151029_n10k_GC1_Psort_WITHSNPLOC.tbl als.2016.txt iq_MOD.txt  # Alzh_IGAP_stage_1.txt BIPOLAR_daner_PGC_BIP32b_mds7a.BIP25_noboma.txt scz2.snp.results.txt Eating.Disorders.pgc.an.snp.all.13May2016.txt GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt ADHD.pgc.adhd.full.2012-10.txt MDD.pgc.mdd.full.2012-04.txt EduYears_Main.txt SubjectiveWellBeing.SWB_Full.txt.gz DepressiveSymptoms.DS_Full.txt.gz SSGAC_EduYears_Rietveld2013_publicrelease.txt EUR.CD.gwas_info03_filtered.assoc EUR.IBD.gwas_info03_filtered.assoc EUR.UC.gwas_info03_filtered.assoc schiz.qjecp.clozukpgc2.gz parkinsons.phs000089.pha002868.clean.txt
for GWASFILE in IQ.Sniekers.2017.txt EUR.IBD.gwas_info03_filtered.assoc_edited GWAS_EA_excl23andMe.txt GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq_WITHSNPLOC.txt ParkinsonsSumStats2018_23andMe.tab 2040_RiskTaking_UKBB.assoc.tsv 2030_GuiltyFeelings_UKBB.assoc.tsv iPSYCH-PGC_ASD_Nov2017 CLOZUK_PGC2noclo.METAL.assoc.dosage.fix Alzh_IGAP_stage_1.txt
do
	if [ ! -d "${LDSCSRC}/Results/${GWASFILE}" ]; then
	  mkdir "${LDSCSRC}/Results/${GWASFILE}"
	fi
	FILES=($(ls $CTDFOLDER | grep "celltype_data"))
	for inputfile in "${FILES[@]}"
	do
		INPUTFILE=$(echo ${inputfile} | sed -e 's/\(.*\)\.rda/\1/')
		if [ ! -d "${LDSCSRC}/Results/${GWASFILE}/${INPUTFILE}" ]; then
		  mkdir "${LDSCSRC}/Results/${GWASFILE}/${INPUTFILE}"
		fi	
		export RCALL="Rscript $LDSCCODE/get_celltype_name.r "$INPUTFILE
		export CELLTYPES=($($RCALL))			
		
		for STEPSIZE in 0.1 # 0.025
		do
			export ROOTANNOT=$INPUTFILE"_AnnotLevel_1"
			for CELLTYPE in "${CELLTYPES[@]}"
			do
				export QUEUE=normal

				export LOGNAME=${ROOTANNOT}_${GWASFILE}_${CELLTYPE}_CutOff${THRESH}_${SPLITMODE}_StepSize_${STEPSIZE}
				export TAGNAME=${ROOTANNOT}_${CELLTYPE}_CutOff${THRESH}_${SPLITMODE}_StepSize_${STEPSIZE}
				export FOLDERNAME=${ANNOTPATH}${INPUTFILE}/${TAGNAME}
				export ANNOTNAME=${ROOTANNOT}_${CELLTYPE}_CutOff${THRESH}_${SPLITMODE}_StepSize_${STEPSIZE}
				export RESFOLDER=$LDSCSRC/Results/${GWASFILE}/${INPUTFILE}

				# export CMDSTR="/nas/longleaf/home/nskene/anaconda_ete/envs/envLDSC/bin/python $HOME/LDSC/ldsc.py \
				# --h2 $MUNGED_FOLDER/${GWASFILE}.sumstats.gz\
			 #    --w-ld-chr $LDSCPROG/weights_hm3_no_hla/weights.\
			 #    --ref-ld-chr ${FOLDERNAME}/${ANNOTNAME}.,$LDSCPROG/1000G_EUR_Phase3_baseline/baseline.\
			 #    --overlap-annot\
			 #    --frqfile-chr $LDSCPROG/1000G_Phase3_frq/1000G.EUR.QC.\
			 #    --out $RESFOLDER/${GWASFILE}_${INPUTFILE}_annot1.$CELLTYPE.${SPLITMODE}.Thresh${THRESH}.Stepsize${STEPSIZE}\
			 #    --print-coefficients"
				export CMDSTR="/nas/longleaf/home/nskene/anaconda_ete/envs/envLDSC/bin/python $HOME/LDSC/ldsc.py \
				--h2 '$MUNGED_FOLDER/${GWASFILE}.sumstats.gz'\
			    --w-ld-chr $LDSCPROG/weights_hm3_no_hla/weights.\
			    --ref-ld-chr ${FOLDERNAME}/${ANNOTNAME}.,$LDSCPROG/1000G_EUR_Phase3_baseline/baseline.\
			    --overlap-annot\
			    --frqfile-chr $LDSCPROG/1000G_Phase3_frq/1000G.EUR.QC.\
			    --out '$RESFOLDER/${GWASFILE}_${INPUTFILE}_annot1.$CELLTYPE.${SPLITMODE}.Thresh${THRESH}.Stepsize${STEPSIZE}'\
			    --print-coefficients"			    
				
			    # Check if the results file already exists
			    export RESFILE="$RESFOLDER/${GWASFILE}_${INPUTFILE}_annot1.$CELLTYPE.${SPLITMODE}.Thresh${THRESH}.Stepsize${STEPSIZE}.results"
			    if [ -f "$RESFILE" ]
		    	then
		    		alreadyExists=$(($alreadyExists + 1))
		    		#echo "RESULTS ALREADY EXISTS"
		    		#echo $RESFILE
		    	else
		    		echo ${FOLDERNAME}
					# First, check if the annotation files have been succesfully prepped for each chromosome for this celltype
					ALLGOOD=1
					for (( i=1; i<=22; i=i+1 ))
					do
						if [ -f "${FOLDERNAME}/${ANNOTNAME}.$i.l2.ldscore.gz" ];
						then				    
							export TEST=1					
						else
							ALLGOOD=0
						fi				
					done
					if [ $ALLGOOD == 1 ]
					then
						#echo "READY"
						export TEST=1	
						# CHECK IF IT IS ALREADY QUEUED
						# ... Note, in the below commands the unix IFS (Internal Field Separator) variable is manipulated so that it becomes end of line
						# ... This enables bjobs to be captured as an array, but it must be reset afterwards
						#oldifs="$IFS"
						#IFS=$'\n'
						#myJobs=($(bjobs -w | grep ${GWASFILE}_${INPUTFILE}_annot1.$CELLTYPE.${SPLITMODE}.Thresh${THRESH}.Stepsize${STEPSIZE}))
						#IFS="$oldifs"
						
						#if [ ${#myJobs[@]} -gt 0 ]
						#then
							#echo "JOB ALREADY QUEUED"
						#	count=$(($count + 1))
						#else
							#echo "JOB ADDED TO QUEUE---------------------"
					    	#bsub -q $QUEUE -R 'select[type==X86_64 && mem>16000] rusage[mem=16000]' -M 16000 -o $TEMP'/test_'${LOGNAME}'.output' -e $TEMP'/test_'${LOGNAME}'.error'  $CMDSTR
					    	#echo $CMDSTR
					    	export LOGOUT=$TEMP'/test_'${LOGNAME}'.output' 
					    	export LOGERR=$TEMP'/test_'${LOGNAME}'.error'
					    	sbatch -t 120 -n 1 --mem=16000 -o $LOGOUT -e $LOGERR --wrap="${CMDSTR}"
						 #   export TEST=1						
							countAdded=$(($countAdded + 1))						    		
						#fi						
					else
						#echo "MISSING-----"
						needsPrep=$(($needsPrep + 1))
					fi
				fi
				
				#bsub -q $QUEUE -R 'select[type==X86_64 && mem>16000] rusage[mem=16000]' -M 16000 -o $TEMP'/test_'${LOGNAME}'.output' -e $TEMP'/test_'${LOGNAME}'.error'  $CMDSTR
			done
		done
	done
done
echo "${alreadyExists} results files already exist"
#echo "${count} jobs already in the queue"
echo "${countAdded} jobs added to the queue"
echo "${needsPrep} still need prep work done"