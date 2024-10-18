#!/bin/sh
#Slurm commands
source home/.bashrc
module load Anaconda3
source activate home/.conda/envs/rbb


# Argument 1: relative path to K_matrix config file
export K_MATRIX_CONFIG='K_matrices/K_matrix_config/TestingBased_Kmat_config_hg38_mode.csv'

# Argument 2: relative path to the config file containing the names of the cfDNA files to deconvolute
export CFDNA_FILES='config/example_project/neu_percentage_deconv_comparisons/deconvolution_config.csv'

#Argument 3: relative path to desired results output folder (should contain deconvolution_results/, evaluation_results/, and logs/)
export OUTPUT_FOLDER='results/example_project/deconvolution_results'

#Argument 4: the name you want to give your deconvolution results file
export DECONV_RES_NAME='example_results.csv'

python deconvolution.py $K_MATRIX_CONFIG $CFDNA_FILES $OUTPUT_FOLDER $DECONV_RES_NAME