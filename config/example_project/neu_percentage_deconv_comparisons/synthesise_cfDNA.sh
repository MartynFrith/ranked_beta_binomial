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

#Can choose 'synthesise_poisson' or 'synthesise_zero_inflated_poisson'. NOTE! Baseline percentage of zeros is ~5% due to missing sites in loyfer reference matrix data
export SYNTHESISER='synthesise_poisson'

#Relative path to the K matrix to use to synthesise cfDNA (.nc). Recommended that a K matrix covering ALL CpG sites is used.
export K_MATRIX='K_matrices/K_matrix_hg38_full_mode.nc'

#Relative path to the ground truth file to use to synthesise cfDNA (.csv)
export GROUND_TRUTH='ground_truths/emperical_ground_truth_1percent_neuron_hg38.csv'
export N_INDIVIDUALS=10

#Need for 10, 20, and 30 reads: 5%, 1%, 0.5%, 0.1%, 0.05%, 0.01%. All with 5% zeros in sequencing data. 

export MEAN_READ_DEPTH=10
#10 rds 5%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_5percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#10 rds 1%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_1percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#10 rds 0.5%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_0.5percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#10 rds 0.1%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_0.1percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#10 rds 0.05%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_0.05percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#10 rds 0.01%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_0.01percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS



export MEAN_READ_DEPTH=20
#20 rds 5%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_5percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#20 rds 1%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_1percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#20 rds 0.5%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_0.5percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#20 rds 0.1%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_0.1percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#20 rds 0.05%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_0.05percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#20 rds 0.01%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_0.01percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS



export MEAN_READ_DEPTH=30
#30 rds 5%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_5percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#30 rds 1%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_1percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#30 rds 0.5%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_0.5percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#30 rds 0.1%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_0.1percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#30 rds 0.05%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_0.05percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS

#30 rds 0.01%
export GROUND_TRUTH='ground_truths/emperical_ground_truth_0.01percent_neuron_hg38.csv'
python synthesise_cfDNA.py $SYNTHESISER $K_MATRIX $GROUND_TRUTH $MEAN_READ_DEPTH $N_INDIVIDUALS
