# **Ranked Beta Binomial (RBB)**
(Work in progress!)

Reference-based deconvolution of cell type composition from cfDNA samples using methylation sequencing data. This software utilises a Baysesian approach to account for the missingness of sequencing data used for the construction of cell type deconvolution reference matrices. Detection of the cell type of interest is optimised by calculating relative entropy (Kullback-Leibler divergence) between the target cell type and all other cell types at each CpG site in the genome, and then constructing the reference matrix using the n sites with the highest relative entropy.


**Requirements:**
- Python 3+
- numpy
- pandas
- scipy
- scikit-learn
- xarray
- netCDF4
- multiprocessing
<br/>

A reference matrix of all tissue types to be included in the deconvolution, aligned to the same genome build as your cfDNA sequencing data, must be created before starting. This file can be created using the generate_rmn.py script using cell type specific sequencing .beta files. These .beta files can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5652177 or can be generated from your own sequencing data using the wgbstools software found here: https://github.com/nloyfer/wgbs_tools. Cell type specific sequencing data should be structred as follows:

├── betas  
│   ├── Cell_type1  
│   │   ├── Cell_type1_sample1.beta  
│   │   ├── Cell_type1_sample2.beta  
│   │   ├── Cell_type1_sampleN.beta  
│   ├── Cell_type2  
│   ├── Cell_typeN  
<br/>

To generate the reference matrix, open the generate_rmn.py script and replace the arguments on lines 10 to 12 with the relevant values:
- **source**: file path to the parent folder containing all .beta files.
- **target**: file path to the folder that you want to save ouputs to.
- **cpgFile**: file path to the .bed file containing CpG index information aligned to the same genome build as your cfDNA sequencing data.*
- **rmn_filename**: name of the file to save the rmn xarray to.
- **k_matrix_filename**: name of the file to save the K matrix (reference matrix) to.
*Currently there is no implementation within Ranked Beta Binomial to create this file. Generate this file using the wgbstools init_genome function:  https://github.com/nloyfer/wgbs_tools/blob/master/docs/init_genome_ref_wgbs.md 
<br/>

To run deconvolution, open the rbb.py script and replace the arguments on lines 118 to 123 with the relevant values:
- **config**: file path to a csv formatted file that is structured as follows:

| filenames  | output_tag| 
| ------------- | ------------- |
| path/to/cfDNA_sequencing_data/sample1.beta  | sample1 |
| path/to/cfDNA_sequencing_data/sampleN.beta   | sampleN  |

- **rmn_path**: file path to the rmn xarray.
- **k_matrix_path**: file path to the K matrix.
- **target**: tissue type you want to optimise the deconvolution for.
- **min_samples**: exclude tissue types with fewer than this number of cell type-specific sequencing samples used to construct the reference matrix.
- **min_rd**: exclude CpG sites with fewer than this number of reads in the cfDNA sequencing data from the deconvolution.
- **n_sites**: select this many sites with the highest relative entropy to use in the deconvolution.
- **results_filename**: the file that you want to save the deconvolution results to.

Then run the rbb.py script.

<br/>

Work in progress features that will soon be added include:
- Synthetic cfDNA generator
- K-fold cross validation 
