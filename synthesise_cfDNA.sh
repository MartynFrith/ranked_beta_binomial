#!/bin/sh
#Slurm commands
source home/.bashrc
module load Anaconda3
source activate home/.conda/envs/rbb

# Argument 1: Can choose 'synthesise_poisson' or 'synthesise_zero_inflated_poisson'. NOTE! Baseline percentage of zeros is ~5% due to missing sites in reference matrix data
export SYNTHESISER='synthesise_poisson'

# Argument 2: Relative path to the K matrix to use to synthesise cfDNA (.nc). Recommended that a K matrix covering ALL CpG sites is used.
export K_MATRIX='K_matrices/K_matrix_hg38_full.nc'

# Argument 3: Relative path to the ground truth file to use to synthesise cfDNA (.csv)
export GROUND_TRUTH='ground_truths/emperical_ground_truth_1percent_neuron_hg38.csv'

# Argument 4: Mean read depth across the synthetic cfDNA to use 
export MEAN_READ_DEPTH=30

# Argument 5: Number of cfDNA samples to synthesise
export N_INDIVIDUALS=10

python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS
