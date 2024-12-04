#!/bin/bash 

#SBATCH -D /lustre/home/mjf221/ranked_beta_binomial
#SBATCH -p mrcq
#SBATCH --time=06:00:00
#SBATCH -A Research_Project-MRC190311
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100G
#SBATCH --error=log_files/generate_rmn.err
#SBATCH --output=log_files/generate_rmn.out
#SBATCH --mail-type=END
#SBATCH --mail-user=mjf221@exeter.ac.uk

source activate /lustre/home/mjf221/.conda/envs/pandas-test
python generate_rmn.py