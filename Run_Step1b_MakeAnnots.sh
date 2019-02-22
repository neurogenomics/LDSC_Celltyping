#!/bin/sh
unset PYTHONPATH
cd $LDSCSRC
source $LDSCCODE/InitVars.sh
export TEMP=Logs_Step1_MakeAnnots
mkdir $TEMP
export SPLITMODE=percentile

rcode_file=$LDSCCODE/create_LDSC_annot.r

# GET ALL 'celltype_data' FILES LOCATED IN THE ROOT FOLDER... THEN LOOP OVER THEM
FILES=($(ls $CTDFOLDER | grep "celltype_data"))
for inputfile in "${FILES[@]}"
do
	# STRIP '.rda' FROM THE NAME OF THE celltype_data FILE
	echo $INPUTFILE
	INPUTFILE=$(echo ${inputfile} | sed -e 's/\(.*\)\.rda/\1/')
	export RCALL="Rscript $LDSCCODE/get_celltype_name.r "$INPUTFILE
	export CELLTYPES=($($RCALL))
	export ROOTANNOT=$INPUTFILE"_AnnotLevel_1"

	#echo $RCALL
	#echo $INPUTFILE
	#for CELLTYPE in "${CELLTYPES[@]}"
	#do
	#	echo $CELLTYPE
	#done

	for STEPSIZE in 0.1 #0.025
	do
		for (( TI=1; TI<=22; TI=TI+1 ))
		do
			export THRESH=0
			

			# CHECK ALL FILES WERE GENERATED, FOR EACH CELLTYPE
			export ALLFILES=1
			for CELLTYPE in "${CELLTYPES[@]}"
			do
				export FOLDERNAME=${LDSCSRC}/CellAnnots/${INPUTFILE}/${ROOTANNOT}_${CELLTYPE}_CutOff${THRESH}_${SPLITMODE}_StepSize_${STEPSIZE}	
				export FILENAME=${ROOTANNOT}_${CELLTYPE}_CutOff${THRESH}_${SPLITMODE}_StepSize_${STEPSIZE}.$TI.annot.gz
				if [ ! -d ${FOLDERNAME} ];
				then
					export ALLFILES=0
					#echo "Cannot find "${FOLDERNAME}
				else
					if [ ! -f "${FOLDERNAME}/${FILENAME}" ];
					then
						export ALLFILES=0
						#echo "Cannot find ${FOLDERNAME}/${FILENAME}"
						#echo ${FOLDERNAME}/${FILENAME}

					fi
				fi
			done
			
			if [ "$ALLFILES" -eq 1 ]
			then
			   # ASSUME THAT IF THE FILE ALREADY EXISTS, DON'T REGENERATE	
			   export FECK=1
			   echo "FILES EXIST FOR EACH CELLTYPE"
			else
  			    export QUEUE=normal
  			    # DELETE THE LOG FILES IF THEY EXIST
  			    export LOGOUT="${TEMP}/test_TI${TI}_${INPUTFILE}_${SPLITMODE}_Thresh${THRESH}_StepSize${STEPSIZE}.output"
  			    export LOGERR="${TEMP}/test_TI${TI}_${INPUTFILE}_${SPLITMODE}_Thresh${THRESH}_StepSize${STEPSIZE}.error"	  			    
   				if [ -f "$LOGOUT" ];
				then
					rm $LOGOUT
					rm $LOGERR
				fi
				# THEN GENERATE THE ANNOTATION
				#bsub -q $QUEUE -R 'select[type==X86_64 && mem>8000] rusage[mem=8000]' -M 8000 -o $LOGOUT -e $LOGERR 'R-3.0.0 --file=create_LDSC_annot.r --args --pre_tI '$TI' --cutoff_thresh '$THRESH' --split_mode '$SPLITMODE' --input_file '$INPUTFILE' --step_size '$STEPSIZE
				sbatch -t 1400 -n 1 --mem=8000 -o $LOGOUT -e $LOGERR --wrap="Rscript "$rcode_file" --pre_tI "$TI" --cutoff_thresh "$THRESH" --split_mode "$SPLITMODE" --input_file "$INPUTFILE" --step_size "$STEPSIZE
				echo "DOESNT EXIST!!!!!!!!!!!!!!!!!"
				#echo $FOLDERNAME
			fi
			 
		done
	done
done
