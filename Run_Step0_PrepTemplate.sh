#!/bin/sh
source InitVars.sh

# This file should be a "run once and forget"
# It runs prepare.new.template.matchingGenes.files() from within create_LDSC_annot.r
# It assumes that you have already done the following:
	# Load the massive SNP database
	# NOTE:
	# - Loading the actual SNP database is unrealistic, so preprocess it first:
	# - First download it:
	# -   wget ftp://ftp.ncbi.nih.gov/snp/database/organism_data/human_9606/b147_SNPContigLocusId_107.bcp.gz
	# - Then remove the SNP and Gene ID columns:
	# -   gzip -cd b147_SNPContigLocusId_107.bcp.gz | cut -f 1,7 > SNPPartialLocus.csv
	# - For running tests, if the code is modified, grab the first 50000 lines of SNPPartialLocus.csv
	#     head -n 50000 SNPPartialLocus.csv > SNPPartialLocus_head.csv
# Note, the main file no longer exists though! So I need to totally replace this section of code.

cd $LDSCSRC
export TEMP=Logs_Step0_PrepTemplate
mkdir $TEMP

rcode_file=$LDSCCODE/create_LDSC_annot.r

for (( TI=1; TI<=22; TI=TI+1 ))
do				
    #export QUEUE=normal
	#bsub -q $QUEUE -R 'select[type==X86_64 && mem>32000] rusage[mem=32000]' -M 32000 -o $TEMP'/test_TI'$TI'.output' -e $TEMP'/test_TI'$TI'.error' 'R-3.0.0 --file=create_LDSC_annot.r --args -P -T '$TI
	sbatch -t 120 -n 1 --mem=32000 -o $TEMP'/test_TI'$TI'.output' -e $TEMP'/test_TI'$TI'.error' --wrap="Rscript "$rcode_file" --args -P -T '"$TI"'"
done
