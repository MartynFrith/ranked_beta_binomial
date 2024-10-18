#!/bin/sh
#Slurm commands
source home/.bashrc
module load Anaconda3
source activate home/.conda/envs/rbb

export RMN='reference_data/rmn/rmn_hg38_small.nc'
export TARGET='Neuron' #Neuron for hg38, neuronal for hg19
export MIN_SAMPLES=0
export MIN_RD=0
export MEDIAN_OR_MODE='median'

export N_SITES=100000
python generate_K_matrix.py $RMN $TARGET $MIN_SAMPLES $MIN_RD $N_SITES $MEDIAN_OR_MODE

