#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D /lustre/home/mjf221/entropy_deconv
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=06:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC190311 # research project to submit under
#SBATCH --nodes=8 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --error=log_files/synthesise_cfDNA.err
#SBATCH --output=log_files/synthesise_cfDNA.out
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=mjf221@exeter.ac.uk # email address

source /lustre/home/mjf221/.bashrc
module load Anaconda3
source activate /lustre/home/mjf221/.conda/envs/entropy_deconv

#Can choose 'synthesise_poisson' or 'synthesise_zero_inflated_poisson'
export SYNTHESISER='synthesise_zero_inflated_poisson'

#Relative path to the K matrix to use to synthesise cfDNA (.nc). Recommended that a K matrix covering ALL CpG sites is used.
export K_MATRIX='K_matrices/K_matrix_hg19_full.nc'

#Relative path to the ground truth file to use to synthesise cfDNA (.csv)
export GROUND_TRUTH='ground_truths/emperical_ground_truth_hg19.csv'

export MEAN_READ_DEPTH=10
export N_INDIVIDUALS=1

python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS