#!/bin/sh
#Slurm commands
source home/.bashrc
module load Anaconda3
source activate home/.conda/envs/rbb


#Argument 1
#Relative path to deconvolution results file 
export DECONVOLUTION_RESULTS='results/example_project/deconvolution_results'

#Argument 2
#Relative path to evaluation config file. Should contain 2 columns:
#   'target_ct_proportion' - containing each of the target cell type proportions used in ground truth creation. Must contain no decimal points 
#   'ground_truth_file' - containing the relative path to the ground truth file associated with the target_ct_proportion 

#       example:
#       target_ct_proportion, ground_truth_file
#       1%                  , ground_truths/neuron_1%_ground_truth.csv
#       01%                 , ground_truths/neuron_0.1%_ground_truth.csv

export EVALUATION_CONFIG='config/example_project/neu_percentage_deconv_comparisons/evaluate_config.csv'

python evaluate.py $DECONVOLUTION_RESULTS $EVALUATION_CONFIG