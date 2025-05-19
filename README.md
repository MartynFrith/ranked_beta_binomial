# **Ranked Beta Binomial (RBB)**
(Work in progress!)

Reference-based deconvolution of cell type composition from cfDNA samples using methylation sequencing data. This software utilises a Baysesian approach to account for the missingness of sequencing data used for the construction of cell type deconvolution reference matrices. Detection of the cell type of interest is optimised by calculating the Jensen-Shannon divergence between the target cell type and all other cell types at each CpG site in the genome, and then constructing the reference matrix using the n sites with the highest relative entropy.


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

A reference matrix, or K matrix, of all tissue types to be included in the deconvolution, aligned to the same genome build as your cfDNA sequencing data, must be created before starting. There are two steps towards generating the K matrix.
 
**1**. Generate the rmn.nc file: This file can be created using the generate_rmn.py script using cell type specific sequencing .beta files. These .beta files can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5652177 or can be generated from your own sequencing data using the wgbstools software found here: https://github.com/nloyfer/wgbs_tools. Cell type specific sequencing data should be structred as follows:

├── betas  
│   ├── Cell_type1  
│   │   ├── Cell_type1_sample1.beta  
│   │   ├── Cell_type1_sample2.beta  
│   │   ├── Cell_type1_sampleN.beta  
│   ├── Cell_type2  
│   ├── Cell_typeN  
<br/>

To generate the rmn.nc file, open the generate_rmn.py script and replace the arguments on lines 32 to 35 with the relevant values, then run:
- **source**: file path to the parent folder containing all .beta files.
- **target**: file path to the folder that you want to save ouputs to.
- **cpgFile**: file path to the .bed file containing CpG index information aligned to the same genome build as your cfDNA sequencing data.*
- **rmn_filename**: name of the file to save the rmn xarray to.
*Currently there is no implementation within Ranked Beta Binomial to create this file. Generate this file using the wgbstools init_genome function:  https://github.com/nloyfer/wgbs_tools/blob/master/docs/init_genome_ref_wgbs.md 

**2.** Generate the k_matrix.nc file: This file can be generated using the generate_kMat.py script. The returned file K matrix file features CpG sites (coordinate 'gen') ranked by highest Jensen-Shannon divergence to lowest. To generate the k_matrix.nc file, open the generate_kMat.py script and replace the arguments on lines 46 to 48 with the relevant values:
- **rmn_path**: file path to the previously generated rmn.nc file, then run.
- **target_ct**: name of the cell type you wish base Jensen-Shannon divergence ranking on.
- **ct_min_samples**: the minimum number of samples a cell type must have in rmn for the cell type to be included in the K matrix.
<br/>
<br/>

Finally, to run deconvolution, open the rbb.py script and replace the arguments on lines 89 to 91 with the relevant values:
- **config_path**: file path to a csv formatted file that is structured as follows:

| filenames  | output_tag| 
| ------------- | ------------- |
| path/to/cfDNA_sequencing_data/sample1.beta  | sample1 |
| path/to/cfDNA_sequencing_data/sampleN.beta   | sampleN  |

- **k_matrix_path**: file path to the K matrix.
- **results_filename**: the file that you want to save the deconvolution results to.

Then run the rbb.py script.

<br/>

Work in progress features that will soon be added include:
- Synthetic cfDNA generator
- K-fold cross validation 
