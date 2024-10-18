#!/bin/sh
#Slurm commands
source home/.bashrc
module load Anaconda3
source activate home/.conda/envs/rbb

python generate_rmn.py 