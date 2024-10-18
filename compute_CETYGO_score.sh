#!/bin/sh
#Slurm commands

module load Anaconda3
source activate home/.conda/envs/rbb

#Argment 1: Software used to perform deconvolution
export DECONVOLUTION_SOFTWARE='RBB'
#Path to folder containing bulk sequencing cfDNA files used in deconvolution
export CFDNA_FILES='bulk_seq_data/example_project'
#Path to deconvolution results
export DECONVOLUTION_RESULTS='results/example_project/deconvolution_results/example_results.csv'
#Path to reference matrix file 
export REFERENCE_MATRIX='K_matrices/K_matrix_hg38_small_mode.nc'

Rscript compute_CETYGO_score_test.R $DECONVOLUTION_SOFTWARE $CFDNA_FILES $DECONVOLUTION_RESULTS $REFERENCE_MATRIX