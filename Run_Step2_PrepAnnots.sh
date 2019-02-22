#!/bin/sh
cd $LDSCSRC
unset PYTHONPATH

export TEMP=Logs_Step2_PrepAnnots
mkdir $TEMP
if [ ! -d Results ]; then
  mkdir Results
fi

export ANNOTPATH=$LDSCSRC/CellAnnots/  #"/lustre/scratch117/unknown/team141/ns9/LDSC/snp_annot_files/"

export THRESH=0
export SPLITMODE=percentile

FILES=($(ls $CTDFOLDER | grep "celltype_data"))
count=0
countAdded=0
countNeedMaking=0
countAlreadyPrepped=0
for inputfile in "${FILES[@]}"
do
	INPUTFILE=$(echo ${inputfile} | sed -e 's/\(.*\)\.rda/\1/')

	# Create log folder just for this input file
	export TEMPLOG="${TEMP}/${INPUTFILE}"
	if [ ! -d "${TEMPLOG}" ]; then
	  mkdir "${TEMPLOG}"
	fi

	export RCALL="Rscript $LDSCCODE/get_celltype_name.r "$INPUTFILE
	export CELLTYPES=($($RCALL))
	
	export ROOTANNOT=$INPUTFILE"_AnnotLevel_1"
	for CELLTYPE in "${CELLTYPES[@]}"
	do
		for STEPSIZE in 0.1 
		do
			for (( i=1; i<=22; i=i+1 ))
			do		
				#export QUEUE=normal
				export FOLDERNAME=${ANNOTPATH}${INPUTFILE}/${ROOTANNOT}_${CELLTYPE}_CutOff${THRESH}_${SPLITMODE}_StepSize_${STEPSIZE}
				export ANNOTNAME=${ROOTANNOT}_${CELLTYPE}_CutOff${THRESH}_${SPLITMODE}_StepSize_${STEPSIZE}			

				#export CMDSTR="$HOME/env/bin/python $HOME/LDSC/ldsc.py\
				export CMDSTR="/nas/longleaf/home/nskene/anaconda_ete/envs/envLDSC/bin/python $HOME/LDSC/ldsc.py\
				    --l2\
				    --bfile $HOME/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i\
				    --ld-wind-cm 1\
				    --annot ${FOLDERNAME}/${ANNOTNAME}.$i.annot.gz\
				    --out ${FOLDERNAME}/${ANNOTNAME}.$i\
				    --print-snps $HOME/LDSC/hapmap3_snps/hm.$i.snp"
				
				export LOGOUT=$TEMPLOG'/test_'$i'_'$CELLTYPE'_'$INPUTFILE'_'$SPLITMODE'_Thresh'$THRESH'_StepSize'$STEPSIZE'.output'
				export LOGERR=$TEMPLOG'/test_'$i'_'$CELLTYPE'_'$INPUTFILE'_'$SPLITMODE'_Thresh'$THRESH'_StepSize'$STEPSIZE'.error'					
				
				# First, check that Step1 (MakeAnnots) has completed for this celltype
				if [ -f "${FOLDERNAME}/${ANNOTNAME}.$i.annot.gz" ];
				then				    
					# Then, check whether the celltype has already been 'prepped'
				    #echo ${ANNOTNAME}.$i.l2.ldscore.gz;
					if [ -f "${FOLDERNAME}/${ANNOTNAME}.$i.l2.ldscore.gz" ];
					then				    
						#echo "Annotation already prepped"
						export TEST=1	
						countAlreadyPrepped=$(($countAlreadyPrepped + 1))								
					else
						# Then check if it is already queued
						# ... Note, in the below commands the unix IFS (Internal Field Separator) variable is manipulated so that it becomes end of line
						# ... This enables bjobs to be captured as an array, but it must be reset afterwards
						#oldifs="$IFS"
						#IFS=$'\n'
						#myJobs=($(bjobs -w | grep ${FOLDERNAME}/${ANNOTNAME}.$i.))
						#IFS="$oldifs"
						
						#if [ ${#myJobs[@]} -gt 0 ]
						#then
							#echo "JOB ALREADY QUEUED"
							#count=$(($count + 1))
						#else
							#echo "JOB ADDED TO QUEUE---------------------"
					    	#bsub -n2 -q $QUEUE -R 'select[type==X86_64 && mem>8000] rusage[mem=8000] span[hosts=1]' -M 8000 -o $LOGOUT -e $LOGERR $CMDSTR
					    	sbatch -t 120 -n 2 --mem=8000 -o $LOGOUT -e $LOGERR --wrap="${CMDSTR}"
						    echo "Annotation needs to be prepped"
						    export TEST=1						
							countAdded=$(($countAdded + 1))						    		
						#fi	

					fi
				else
				    #echo ${FOLDERNAME}/${ANNOTNAME}.$i.annot.gz;
  				    echo "Annotation not generated -- repeat Step1" 
  				    export TEST=1			
  				    countNeedMaking=$(($countNeedMaking + 1))	
  				    #echo "${FOLDERNAME}/${ANNOTNAME}.$i.annot.gz"
  				    #break;		
				fi
			done
		done
	done
done
#echo "${count} jobs already in the queue"
echo "${countAlreadyPrepped} jobs already prepared"
echo "${countAdded} jobs added to the queue"
echo "${countNeedMaking} jobs still need their annotations generated"
