#%%
import sys
import xarray as xr
import numpy as np
import pandas as pd
from scipy.stats import beta as beta_distribution
from scipy.stats import entropy

#%% OPTIONS
rmn = sys.argv[1]
target = sys.argv[2]
min_samples = int(sys.argv[3])
min_rd = int(sys.argv[4])
n_sites = int(sys.argv[5])


#%% GENERATE MASK FOR READS
rmn = xr.open_dataarray(rmn)
n = rmn.loc[:,:,'n']
ct_mask = np.any(np.greater(n, min_samples), axis=1)
ct_selection = rmn.loc[ct_mask,:]
mask = np.all(np.greater(ct_selection.loc[:,:,'r'], min_rd), axis=0)
selection = ct_selection.loc[:, mask] #This is the xarray that we need to take forward
muDNAm = np.divide(selection.loc[:,:,'m'], selection.loc[:,:,'r'])
muDNAm = muDNAm.fillna(0)


#%% GENERATE K MATRIX 
selection = xr.open_dataarray('data/rmn/rmn_hg19.nc')
m = selection.loc[:,:,'m']
r = selection.loc[:,:,'r']
pseudo_K =  xr.DataArray(beta_distribution.median(1 + m, 1 + (r - m)), #I changed to median so each version of K matrix has only one probability at each site in
                    dims = ('ct', 'gen'), coords = {'ct': selection.coords['ct'], # each cell type
                                                    'gen':selection.coords['gen']})
# pseudo_K.to_netcdf("K_matrices/K_matrix_hg19_full.nc") 


# #%% RANK SITES BY RELATIVE ENTROPY
chance_target_methylated = np.array(selection.loc[target,:,'m'] / selection.loc[target,:,'r'])
chance_target_not_methylated = 1 - chance_target_methylated
p_target = np.array([chance_target_methylated, chance_target_not_methylated])
others = list(selection.coords['ct'].values)
others.remove(target)
p_others = selection.loc[others, :, 'm'].sum(axis=0) / selection.loc[others, :, 'r'].sum(axis=0)
p_others = np.array([p_others, 1-p_others])
kl = entropy(p_target, p_others)
entropy_rank = np.argsort(kl) 


#%% FILTER PSEUDO K MATRIX BY ENTROPY RANKED SITES TO GET SAMPLED K MATRIX
sel = entropy_rank[-n_sites:]
entropy_mask = np.isin(np.arange(entropy_rank.shape[0]), sel) 
entropy_filtered = selection[:, entropy_mask, :]
sampled_K = pseudo_K[:, entropy_mask]

sampled_K.to_netcdf("K_matrices/K_matrix_hg19_{}_{}rds.nc".format((int(n_sites)), min_rd)) #Sampled K saved
# mask.to_netcdf("masks/gen_mask_hg38_{}k_{}rds.nc".format((int(sites/1000)), min_rd)) #gen_mask saved 
# np.save("masks/entropy_mask_hg38_{}k_{}rds.npy".format((int(sites/1000)), min_rd), entropy_mask) #entropy mask saved
