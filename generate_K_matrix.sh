#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D /lustre/home/mjf221/entropy_deconv
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=06:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC190311 # research project to submit under
#SBATCH --nodes=4# specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --error=log_files/generate_k_matrix.err
#SBATCH --output=log_files/generate_k_matrix.out
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=mjf221@exeter.ac.uk # email address

source /lustre/home/mjf221/.bashrc
module load Anaconda3
source activate /lustre/home/mjf221/.conda/envs/entropy_deconv

export RMN='data/rmn/rmn_hg38.nc'
export TARGET='Neuron' #Neuron for hg38, neuronal for hg19
export MIN_SAMPLES=0
export MIN_RD=10
export MEDIAN_OR_MODE='median'

export N_SITES=100000
python generate_K_matrix.py $RMN $TARGET $MIN_SAMPLES $MIN_RD $N_SITES $MEDIAN_OR_MODE

export N_SITES=129154
python generate_K_matrix.py $RMN $TARGET $MIN_SAMPLES $MIN_RD $N_SITES $MEDIAN_OR_MODE

export N_SITES=166810
python generate_K_matrix.py $RMN $TARGET $MIN_SAMPLES $MIN_RD $N_SITES $MEDIAN_OR_MODE

export N_SITES=215443
python generate_K_matrix.py $RMN $TARGET $MIN_SAMPLES $MIN_RD $N_SITES $MEDIAN_OR_MODE

export N_SITES=278255
python generate_K_matrix.py $RMN $TARGET $MIN_SAMPLES $MIN_RD $N_SITES $MEDIAN_OR_MODE