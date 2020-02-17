What is contained in the repo:
The primary function performed is to check designated folders for celltype data files & GWAS summary statistics. LDSC annotations are then created, for all celltype data files. All summary statistics within the folder are munged. Heritability is partitioned for each decile of specificity, for each celltype, for each GWAS. Results are then collated and graphed.

What's the order to run scripts to perform the main analysis:
    InitVars.sh
    Run_Step0_PrepTemplate.sh
    Run_Step1a_CheckCelltypeNames.sh
    Run_Step1b_MakeAnnots.sh
    Run_Step2_PrepAnnots.sh
    Run_Step2b_PrepBaselineAnnots.sh
    Run_Step3a_MungeSumStats.sh
    Run_Step4_CellTypeAnalysis.sh
    Run_Step5_GetSummaryStatsOnServer.r

Ignore anything starting with:
    RunGeneList_

What are the main scripts:
1. InitVars.sh
    1. This file needs to be run to set paths
2. create_LDSC_annot.r
    1. This has most of the main R functions used to create the LDSC annotations
3. Run_Step0_PrepTemplate.sh
    1. Needs to be run as a one off... has some explanation inside it's comments

Expected folder structure:
* Creates some folders in $HOME and others on $SCRATCH
* GWAS summary statistics files:
    * export SUMSTATS_FOLDER=$HOME/GWAS_SumStats/Edited
    * export MUNGED_FOLDER=$SCRATCH/GWAS_SumStats/Munged

Known problems:
* InitVars script needs to be modified to create folders as well as just creating variables pointing to them
    * The info in InitVars.sh would tell you which folders are expected to be there
    * This script should also download the LDSC software and associated files, if they are not already there
* InitVars.sh has a hard coded path which needs to be replaced:
    * alias pythonenv="/nas/longleaf/home/nskene/env/bin/python"
* I think it still expects the input files to use the old EWCE celltype_data format, rather than the current CTD format. If so, this needs to be changed... Yes... it does
    * That is, it has celltype_data[[1]]$all_scts instead of ctd[[1]]$specificity
    * To start with, just use an old EWCE format file, you can convert these using the EWCE::convert.... functions
* The code was designed to run on the UCSC Longleaf server, so the shell scripts which submit jobs look like:
    * sbatch -t 1400 -n 1 --mem=8000 -o $LOGOUT -e $LOGERR --wrap="Rscript "$rcode_file" --pre_tI "$TI" --cutoff_thresh "$THRESH" --split_mode "$SPLITMODE" --input_file "$INPUTFILE" --step_size "$STEPSIZE
    * In the first instance these commands should be adjusted to use the PBS submission system used on the Imperial cluster
* The file used to map SNPs to genes is not available at the download link anymore:
    * ftp://ftp.ncbi.nih.gov/snp/database/organism_data/human_9606/b147_SNPContigLocusId_107.bcp.gz
    * This file should be a direct replacement:
        * ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/organism_data/b151_SNPContigLocusId_108.bcp.gz
    * Note: link was found here:
        * ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/database/organism_data/
    * Description of file contents is here:
        * https://www.ncbi.nlm.nih.gov/projects/SNP/snp_db_table_description.cgi?t=SNPContigLocusId
    * To test it, download it, then run "Run_Step0_PrepTemplate.sh".... check for errors
* However... that way of mapping SNPs to genes is really not ideal
    * It's really slow
    * It doesn't allow for control of how large the window is around a gene
* So instead...
    * Use biomaRt to get the locations of genes (matching genome build with the input GWAS)
    * Draw a window around that gene (i.e. 10kb upstream / 10 kb downstream)
    * Assign SNPs within those windows as within that gene
    * Zijing know's a python function that can do this quickly

Suggested plan:
* Instead of modifying InitVars, instead move towards usage of Docker:
    * Use this dockerfile as the template:
        * https://hub.docker.com/r/manninglab/ldsc/dockerfile
* The manninglab has already created a WDL workflow for partitioned heritability:
    * https://github.com/manning-lab/ldsc_wdl
    * Fork this and build the cell typing aspects onto it
