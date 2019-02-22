#!/bin/sh
cd /nfs/team156/ns9/LDSC
export TEMP=Logs_Step1b_GenerateGORes
mkdir $TEMP
#celltype_data_SCAPT_noLevels_NOthresh celltype_data_MBB_wtLEVELS_thresh0_trim0 celltype_data_allKImouse_noLevels_thresh0_trim0
# celltype_data_allKImouse_level1_thresh0_trim0 celltype_data_allKImouse_level2_thresh0_trim0 celltype_data_allKImouse_noLevels_thresh0_trim0 celltype_data_cortex_and_hippocampus_KI_mouse_level1_thresh0_trim0 celltype_data_cortex_and_hippocampus_KI_mouse_level2_thresh0_trim0 celltype_data_striatal_and_midbrain_KI_mouse_level1_thresh0_trim0 celltype_data_striatal_and_midbrain_KI_mouse_level2_thresh0_trim0 celltype_data_TASIC_level1_thresh1_trim0 celltype_data_TASIC_level2_thresh1_trim0 celltype_data_SCAPT_noLevels_NOthresh
for INPUTFILE in celltype_data_allKImouse_level1_thresh0_trim0 celltype_data_allKImouse_level2_thresh0_trim0 celltype_data_CortexDarmanis_level1_thresh0_trim0 celltype_data_CortexSCAPT_level2_thresh0_trim0 celltype_data_HumanMidbrainKI_level1_thresh0_trim0 celltype_data_HumanMidbrainKI_level2_thresh0_trim0 celltype_data_HumanRPKM_level1_thresh0_trim0 celltype_data_HumanRPKM_level2_thresh0_trim0 celltype_data_StriatumQuake_level1_thresh0_trim0 celltype_data_StriatumQuake_level2_thresh0_trim0 celltype_data_TASIC_level1_thresh0_trim0 celltype_data_TASIC_level2_thresh0_trim0 celltype_data_TasicStriatum_level1_thresh0_trim0 celltype_data_TasicStriatum_level2_thresh0_trim0
do
	# proportion
	for SPLITMODE in percentile 
	do
		for STEPSIZE in 0.1 #0.01
		do
			export THRESH=0
			
										
		    export QUEUE=normal
			bsub -q $QUEUE -R 'select[type==X86_64 && mem>8000] rusage[mem=8000]' -M 8000 -o $TEMP'/test_TI'$TI'_'$INPUTFILE'_'$SPLITMODE'_Thresh'$THRESH'_StepSize'$STEPSIZE'.output' -e $TEMP'/test_TI'$TI'_'$INPUTFILE'_'$SPLITMODE'_Thresh'$THRESH'_StepSize'$STEPSIZE'.error' 'R-3.0.0 --file=create_LDSC_annot.r --args --gen_go --cutoff_thresh '$THRESH' --split_mode '$SPLITMODE' --input_file '$INPUTFILE' --step_size '$STEPSIZE
		done
	done
done