#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D /lustre/home/mjf221/entropy_deconv
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=06:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC190311 # research project to submit under
#SBATCH --nodes=4 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --error=log_files/deconvolution.err
#SBATCH --output=log_files/deconvolution.out
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=mjf221@exeter.ac.uk # email address

source /lustre/home/mjf221/.bashrc
module load Anaconda3
source activate /lustre/home/mjf221/.conda/envs/entropy_deconv


# Argument 1: relative path to K_matrix config file
export K_MATRIX_CONFIG='K_matrices/K_matrix_config/cedric_0inflatedPoisResultsBased_Kmat_config_hg19_mode.csv'

# Argument 2: relative path to file containing the names of the cfDNA files to deconvolute
export CFDNA_FILES='config/real_cfDNA_deconvolution/deconvolution_config.csv'

python deconvolution.py $K_MATRIX_CONFIG $CFDNA_FILES