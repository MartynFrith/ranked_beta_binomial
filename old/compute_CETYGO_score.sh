#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D /lustre/home/mjf221/entropy_deconv
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=06:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC190311 # research project to submit under
#SBATCH --nodes=1# specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mem=100G # specify bytes memory to reserve
#SBATCH --error=log_files/compute_CETYGO_score.err
#SBATCH --output=log_files/compute_CETYGO_score.out
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=mjf221@exeter.ac.uk # email address


module load Anaconda3
source activate /lustre/home/mjf221/.conda/envs/rbb

#Argment 1: Software used to perform deconvolution
export DECONVOLUTION_SOFTWARE='RBB'
#Path to folder containing bulk sequencing cfDNA files used in deconvolution
export CFDNA_FILES='bulk_cfDNA_files/loyfer_cfDNA'
#Path to deconvolution results
export DECONVOLUTION_RESULTS='results/mjf221/deconvolution_results/loyfer_real_cfDNA_hg38_cedric_0inflatedPoisResultsBased_Kmat_config_hg38_mode_2024-10-04_1215-38.csv'
#Path to reference matrix file 
export REFERENCE_MATRIX=''

Rscript compute_CETYGO_score_test.R $DECONVOLUTION_SOFTWARE $CFDNA_FILES $DECONVOLUTION_RESULTS $REFERENCE_MATRIX