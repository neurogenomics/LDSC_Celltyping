#!/bin/bash
source InitVars.sh

cd $SCRATCH/GWAS_SumStats/Raw

##==============================================================================================================================
##==============================================================================================================================

# EduYears
wget http://ssgac.org/documents/EduYears_Main.txt.gz

##==============================================================================================================================
##==============================================================================================================================

# Intelligence (Sneikers et al)
# Link from here: http://ctg.cncr.nl/software/summary_statistics
wget http://ctg.cncr.nl/documents/p1651/sumstats.txt.gz
mv $SCRATCH/GWAS_SumStats/Raw/sumstats.txt.gz $SCRATCH/GWAS_SumStats/Raw/IQ.Sneikers2016.sumstats.txt.gz
zcat $SCRATCH/GWAS_SumStats/Raw/IQ.Sneikers2016.sumstats.txt.gz | head -n 3
# The raw file results in an error with LDSC: "Too many signed sumstat columns. Specify which to ignore with the --ignore flag."
# Readme: http://ctg.cncr.nl/documents/p1651/readme.txt
# Columns:
# Chromosome	position	rsid	ref	alt	MAF	Beta	SE	Zscore	p_value	direction
# 1	100000012	rs10875231	T	G	0.234588	0.000298163293384453	0.00596326586768906	0.05	0.9599	+-++--+-
# 1	100000135	rs114947036	A	T	0.00180123	0.0410026320136857	0.0716829231008491	0.572	0.5676	?-+?????
# I tell it to ignore the Beta column
# IT then gives an error: "ValueError: Could not find A1/A2 columns."
# The LDSC tutorial specifies:
# allele_2:	Allele 2, interpreted as non-ref allele for signed sumstat.
# allele_1:	Allele 1, interpreted as ref allele for signed sumstat.
# So presumably 'alt' should be 'A1' and 'ref' should be 'A2'
gzip -d IQ.Sneikers2016.sumstats.txt.gz
sed -i '1s/ref/A1/' IQ.Sneikers2016.sumstats.txt
sed -i '1s/alt/A2/' IQ.Sneikers2016.sumstats.txt
head -n 3 IQ.Sneikers2016.sumstats.txt


##==============================================================================================================================
##==============================================================================================================================

# Anti-social behaviour (Tielbeek, 2017)
# Link from here: http://ctg.cncr.nl/software/summary_statistics
# Paper: https://jamanetwork.com/journals/jamapsychiatry/article-abstract/2656184
# Split into genders, and combined
# The discovery samples comprised 16400 individuals, whereas the target samples consisted of 9381 individuals (all individuals were of European descent), 
#  including child and adult samples (mean age range, 6.7-56.1 years). Three promising loci with sex-discordant associations were found:
#  (8535 female individuals, chromosome 1: rs2764450, chromosome 11: rs11215217; 
#   7772 male individuals, chromosome X, rs41456347). 
#   Polygenic risk score analyses showed prognostication of antisocial phenotypes in an independent Finnish Crime Study (2536 male individuals and 3684 female individuals) 
#   and shared genetic origin with conduct problems in a population-based sample (394 male individuals and 431 female individuals) but not with conduct disorder in a 
#   substance-dependent sample (950 male individuals and 1386 female individuals) (R2 = 0.0017 in the most optimal model, P = 0.03). Significant inverse genetic correlation 
#   of ASB with educational attainment (r = –0.52, P = .005) was detected.
# From text: "In total, BroadABC has genotypic and phe- notypic data from 25 781 individuals across 8 unique samples and, to our knowledge, is the largest collective sample avail- able to estimate the effects of genome-wide genetic variants for ASB and testing for genetic overlap with other traits"
# The size of cohorts is given here:http://broadabc.ctglab.nl/documents/p12/readme_tielbeek_jamapsychiatry2017_antisocial_behavior_sumstats.txt
#  - Combined N=16400
#  - Male: 7772
#  - Female: 8535
wget http://broadabc.ctglab.nl/documents/p12/BroadABC_METALoutput_Combined.tbl
wget http://broadabc.ctglab.nl/documents/p12/BroadABC_METALoutput_Females.tbl
wget http://broadabc.ctglab.nl/documents/p12/BroadABC_METALoutput__Males.tbl
# It gives the "No objects to concatenate" error during MUNGING
head -n 3 $SCRATCH/GWAS_SumStats/Raw/BroadABC_METALoutput__Males.tbl
# MarkerName Allele1 Allele2 Weight Zscore P-value Direction
# X:121854298 a t 6515.00 -5.616 1.954e-08 ---?-
# 6:82075402 a t 7772.00 -5.117 3.111e-07 -----
Rscript $LDSCCODE/get_snpids_for_markername.r $SCRATCH/GWAS_SumStats/Raw/BroadABC_METALoutput__Males.tbl " "
Rscript $LDSCCODE/get_snpids_for_markername.r $SCRATCH/GWAS_SumStats/Raw/BroadABC_METALoutput_Females.tbl " "
# Combined throws an "embedded nul in string" error
#sed 's/\\0//g' $SCRATCH/GWAS_SumStats/Raw/BroadABC_METALoutput_Combined.tbl > $SCRATCH/GWAS_SumStats/Raw/BroadABC_METALoutput_Combined.tbl
sed -i 's/\x0//g' /pine/scr/n/s/nskene/GWAS_SumStats/Raw/BroadABC_METALoutput_Combined.tbl 
Rscript $LDSCCODE/get_snpids_for_markername.r $SCRATCH/GWAS_SumStats/Raw/BroadABC_METALoutput_Combined.tbl "tab"

##==============================================================================================================================
##==============================================================================================================================

# Height (Wood, 2014)
wget https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz

##==============================================================================================================================
##==============================================================================================================================

# BMI (adjusted for smoking, GIANT consortium)
# FROM: https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files
wget https://portals.broadinstitute.org/collaboration/giant/images/3/3a/BMI.SNPadjSMK.zip
unzip BMI.SNPadjSMK.zip
# Contains BMI.SNPadjSMK.CombinedSexes.AllAncestry.txt
# Contains BMI.SNPadjSMK.CombinedSexes.EuropeanOnly.txt
# Contains BMI.SNPadjSMK.Men.AllAncestry.txt
# Contains BMI.SNPadjSMK.Men.EuropeanOnly.txt
# Contains BMI.SNPadjSMK.Women.AllAncestry.txt
# Contains BMI.SNPadjSMK.Women.EuropeanOnly.txt
# N: has an 'N' column in each file
# ERROR! 
# Munging gives the following error: "ValueError: No objects to concatenate"
# This pagea implies it is due to an empty column: https://github.com/bulik/ldsc/issues/66
# There doesn't appear to be any empty columns though.
# Perhaps the error is that it has both 'rs_id' and 'markername' as columns where markername has values like "chr1:161003087"
# Added "export IGNORE="markername"
# It then says it cannot find rsids
sed -i '1s/rs_id/snpid/' BMI.SNPadjSMK.CombinedSexes.EuropeanOnly.txt
sed -i '1s/rs_id/snpid/' BMI.SNPadjSMK.Men.EuropeanOnly.txt
sed -i '1s/rs_id/snpid/' BMI.SNPadjSMK.Women.EuropeanOnly.txt
head -n 3 $SCRATCH/GWAS_SumStats/Raw/BMI.SNPadjSMK.CombinedSexes.EuropeanOnly.txt

##==============================================================================================================================
##==============================================================================================================================

## UK Biobank Risk Taking
wget https://www.dropbox.com/s/27yrfcsngenzab1/2040.assoc.tsv.gz?dl=0%20-O%202040.assoc.tsv.gz
mv 2040.assoc.tsv.gz UKBIOBANK.RiskTaking.2040.assoc.tsv.gz

zcat UKBIOBANK.RiskTaking.2040.assoc.tsv.gz|head -n 3

#[nskene@longleaf-login1 Raw]$ zcat UKBIOBANK.RiskTaking.2040.assoc.tsv.gz|head -n 3
#variant	rsid	nCompleteSamples	AC	ytx	beta	se	tstat	pval
#5:43888254:C:T	rs13184706	325821	1.19508e+04	3.01844e+03	-2.60013e-03	3.98586e-03	-6.52338e-01	5.14184e-01
#5:43888493:C:T	rs58824264	325821	2.34249e+03	5.74553e+02	-1.00442e-02	8.88924e-03	-1.12993e+00	2.58508e-01


##==============================================================================================================================
##==============================================================================================================================

## UK Biobank: Collaboration with David Hill and Ian Deary
## https://www.amazon.co.uk/clouddrive/share/xRX8AAvbA5WZkd8bKv3rCBFQWqV1IKEQC1HUkH2K0Xs
## 9va
 gunzip HI.txt.gz
 gunzip Education_trait_1.txt.gz
 sed -i '1s/rsid/SNPID/' HI.txt
 sed -i '1s/snpid/SNPID/' Education_trait_1.txt
 sed -i '1s/pos/BP/' HI.txt
 sed -i '1s/pos/BP/' Education_trait_1.txt
 sed -i '1s/chr/CHR/' HI.txt
 sed -i '1s/chr/CHR/' Education_trait_1.txt
 sed -i '1s/z/zscore/' Education_trait_1.txt
 sed -i '1s/mtag_z/Z/' Education_trait_1.txt
 sed -i '1s/mtag_pval/P/' Education_trait_1.txt
 sed -i '1s/HI_rsd_beta/BETA/' HI.txt
 sed -i '1s/bBP/BP/' Education_trait_1.txt
 sed -i '1s/zscore/PREzscore/' Education_trait_1.txt
 sed -i '1s/a_1/EFFECT_ALLELE/' HI.txt
 sed -i '1s/a_0/NON_EFFECT_ALLELE/' HI.txt