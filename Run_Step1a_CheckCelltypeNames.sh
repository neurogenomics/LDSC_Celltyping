#!/bin/bash
source InitVars.sh

cd $LDSCSRC
export TEMP=Logs_Step1_CheckCelltypeNames
#mkdir $TEMP
export SPLITMODE=percentile

rcode_file=$LDSCCODE/create_LDSC_annot.r

echo $SCRATCH/Celltype_data
echo $CTDFOLDER

# GET ALL 'celltype_data' FILES LOCATED IN THE ROOT FOLDER... THEN LOOP OVER THEM
FILES=($(ls $CTDFOLDER | grep "celltype_data"))
for inputfile in "${FILES[@]}"
do
	echo $inputfile

	# STRIP '.rda' FROM THE NAME OF THE celltype_data FILE
	INPUTFILE=$(echo ${inputfile} | sed -e 's/\(.*\)\.rda/\1/')

	# THEN GENERATE THE ANNOTATION
	#bsub -q normal -R 'select[type==X86_64 && mem>8000] rusage[mem=8000]' -M 8000 -o ${TEMP}/check_${INPUTFILE}.o -e ${TEMP}/check_${INPUTFILE}.e 'R-3.0.0 --file=create_LDSC_annot.r --args -N --input_file '$INPUTFILE
	sbatch -t 120 -n 1 --mem=8000 -o ${TEMP}/check_${INPUTFILE}.o -e ${TEMP}/check_${INPUTFILE}.e --wrap="Rscript "$rcode_file" -N --input_file '"$INPUTFILE"'"
	#sbatch -t 120 -n 1 --mem=8000 -o ${TEMP}/check_${INPUTFILE}.o -e ${TEMP}/check_${INPUTFILE}.e --wrap="R CMD BATCH "$rcode_file" -N --input_file '"$INPUTFILE"'"
	#echo $INPUTFILE
done