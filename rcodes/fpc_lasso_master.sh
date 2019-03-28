#!/bin/bash

# Base directory
dir_base="/home/erik/Documents/projects/FPC_Lasso/current"
# Code directory
dir_codes=$dir_base"/codes"
# Output directory
dir_output=$dir_base"/output"
# Figure directory
dir_figs=$dir_output/figures
# Dataset directory
dir_datasets=$dir_base"/datasets"
# Output for data directory
dir_output_data=$dir_output"/data"

cd $dir_base

#####################################################
####### ---- (1) Q-Q Plot for rwSNS ----- ###########

# rs_sim_qqplots.R
# input:
# 	(i) dir_codes: loads in rs_dgp.R
# 	(ii) dir_output: folder to save images
# output:
# 	(i) sim_qq_rwsns.csv (csv output of simulation results)
# 	(ii) gg_qq_rwsns.png (qq-plot)

Rscript $dir_codes/rs_sim_qqplots.R --dir_codes $dir_codes --dir_output $dir_output

# NOTE: Script takes apporiximately 1 minute to run

################################################################
####### ---- (2) Survival Dataset construction ----- ###########

# (a) Chandan Reddy datasets: http://dmkd.cs.vt.edu/projects/survival/data/

# rs_dataset_survival.R
# input:
# 	(i) dir_dataset: folder where subfolders containing the different data can be found
# 	(ii) dir_output: folder to output
# outupt:
# 	(i) lst.surv: list of length 81, where each element is a list with three elements: X, So, and cr
# 					where X is the design matrix, So is the Surv() object, and cr is the competing risk information

Rscript $dir_codes/rs_dataset_survival.R --dir_dataset $dir_dataset --dir_output $dir_output

#################################################################################
####### ---- (3) Classification/Regression Dataset construction ----- ###########

# ps_dataset_regclass.py
# input:
# 	(i) dir_output: folder where the CSVs will be temporarily stored
# outupt:
# 	(i) A total of Z X.csv and y.csv files
python $dir_codes/ps_dataset_regclass.py --dir_output $dir_output_data

# rs_dataset_reglass.R
# input:
# 	(i) dir_output: folder where the CSVs have been stored from ps_dataset_regclass.py
# outupt:
# 	(i) A regclass_datasets.RData file
Rscript $dir_codes/rs_dataset_regclass.R --dir_output $dir_output_data


########################################################
####### ---- (4) Regularity conditions ----- ###########

# rs_regularity_conditions.R
# outupt:
# 	(i) 
# 	(ii) /output/gg_lamseq_surv.png || /output/gg_lamseq_class.png || /output/gg_lamseq_reg.png
Rscript $dir_codes/rs_regularity_conditions.R --dir_codes $dir_codes --dir_output $dir_output --dir_data $dir_output_data



#########################################################
####### ---- (5) False Positive Control ----- ###########

# (A) STOCHASTICALLY INDEPENDENT DESIGN

# rs_sim_fp_indep.R
# outupt:
# 	(i) /output/storage_fp_indep.csv
# 	(ii) /output/gg_fp_indep.png
Rscript $dir_codes/rs_sim_fp_indep.R --dir_codes $dir_codes --dir_output $dir_output


# (B) CORRELATED NOSIE COLUMNS IS CONSERVATIVE
# rs_sim_fp_corr.R
# outupt:
# 	(i) /output/storage_fp_corr.csv
# 	(ii) /output/gg_fp_corr.png
Rscript $dir_codes/rs_sim_fp_corr.R --dir_codes $dir_codes --dir_output $dir_output

##############################################################
####### ---- (6) RUN SOME INFERENCE EXAMPLES ----- ###########

# INFERENCE COMPARISONS: selectiveInference + knockoffs + G'Sell
Rscript $dir_codes/rs_inference_comp.R --dir_codes $dir_codes --dir_output $dir_output --dir_data $dir_output_data



######################################################
####### ---- (7) PLOT OUTPUT RESULTS ----- ###########

# rs_fig_create.R
# outupt:
# 	(i)
Rscript $dir_codes/rs_fig_create.R --dir_output $dir_output --dir_figs $dir_figs





