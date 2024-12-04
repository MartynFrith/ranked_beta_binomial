#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D /lustre/home/mjf221/entropy_deconv
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=06:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC190311 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --error=log_files/evaluate.err
#SBATCH --output=log_files/evaluate.out
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=mjf221@exeter.ac.uk # email address

source /lustre/home/mjf221/.bashrc
module load Anaconda3
source activate /lustre/home/mjf221/.conda/envs/rbb

#Argument 1
#Relative path to deconvolution results file 
export DECONVOLUTION_RESULTS='results/mjf221/deconvolution_results'

#Argument 2
#Relative path to evaluation config file. Should contain 2 columns:
#   'target_ct_proportion' - containing each of the target cell type proportions used in ground truth creation. Must contain no decimal points 
#   'ground_truth_file' - containing the relative path to the ground truth file associated with the target_ct_proportion 

#       example:
#       target_ct_proportion, ground_truth_file
#       1%                  , ground_truths/neuron_1%_ground_truth.csv
#       01%                 , ground_truths/neuron_0.1%_ground_truth.csv

export EVALUATION_CONFIG='config/neu_percentage_deconv_comparisons/evaluate_config.csv'

python evaluate.py $DECONVOLUTION_RESULTS $EVALUATION_CONFIG