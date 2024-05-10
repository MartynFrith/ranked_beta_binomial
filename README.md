BAysian CEll TYpe DEconvolution (BACTYDE)

Reference-based deconvolution of cell type composition from cfDNA samples using methylation sequencing data. This software uitlises a Baysesian approach to account for missingness of 
sequencing data used for training data sets when constructing reference altases for cell type deconvolution. Detection of cell type(s) of interest is optimised by calculating relative 
entropy (Kullback-Leibler divergence) between the target cell type and all other cell types at each CpG site in the genome and selecting the n sites with the highest relative entropy. 
