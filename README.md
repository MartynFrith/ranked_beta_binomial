Ranked Beta Binomial (RBB)

Reference-based deconvolution of cell type composition from cfDNA samples using methylation sequencing data. This software utilises a Baysesian approach to account for the missingness of sequencing data used for the construction of cell type deconvolution reference atlases. Detection of cell type(s) of interest is optimised by calculating relative entropy (Kullback-Leibler divergence) between the target cell type and all other cell types at each CpG site in the genome, and then constructing the reference atlas using the n sites with the highest relative entropy.

Requirements:

Python 3+
numpy
pandas
scipy
scikit-learn
xarray
netCDF4
Workflow: A reference matrix of all tissue types to be included in the deconvolution, aligned to the same genome build as your cfDNA sequencing data must be created before starting. This file can be created using the 'generate_rmn' scripts using cell type specific sequencning .beta files. These .beta files can generated using the wgbstools software found here: https://github.com/nloyfer/wgbs_tools.

To run deconvolution, replace the arguments on lines 118 to 123 with your chosen values:

rmn_path: the local path to your reference matrix.
target: the tissue type you want to optimise the deconvolution for.
min_samples: exclude tissue types with fewer than this number of cell type-specific sequencing samples used to construct the reference matrix.
min_rd: exclude CpG sites with fewer than this number of reads in the cfDNA sequencing data from the deconvolution.
n_sites: select this many sites with the highest relative entropy to use in the deconvolution.
results_filename: the file that you want to save the deconvolution results to.
Then run the python script.
