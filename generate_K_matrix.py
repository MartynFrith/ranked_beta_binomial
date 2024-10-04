import numpy as np
from scipy.stats import beta as beta_distribution
import xarray as xr
from scipy.stats import entropy
import sys 

rmn = sys.argv[1]
target = sys.argv[2]
min_samples = int(sys.argv[3])
min_rd = int(sys.argv[4]) 
n_sites = int(sys.argv[5])
median_or_mode = sys.argv[6]
#CREATE ARGUMENT FOR FULL K MATRIX

#%% FILTER SITES BY READ DEPTH 
rmn = xr.open_dataarray(rmn)
n = rmn.loc[:,:,'n']
ct_mask = np.any(np.greater(n, min_samples), axis=1)
ct_selection = rmn.loc[ct_mask,:]
gen_mask = np.all(np.greater(ct_selection.loc[:,:,'r'], min_rd), axis=0) #comment out for full matrix
selection = ct_selection.loc[:, gen_mask] #This is the xarray that we need to take forward!comment out for full matrix


#%% RANK SITES BY ENTROPY

chance_target_methylated = np.array(selection.loc[target,:,'m'] / selection.loc[target,:,'r'])
chance_target_not_methylated = 1 - chance_target_methylated
p_target = np.array([chance_target_methylated, chance_target_not_methylated])
others = list(selection.coords['ct'].values)
others.remove(target)
p_others = selection.loc[others, :, 'm'].sum(axis=0) / selection.loc[others, :, 'r'].sum(axis=0)
p_others = np.array([p_others, 1-p_others])
kl = entropy(p_target, p_others)
entropy_rank = np.argsort(kl) #Entropy ranked sites


#%% CREATE FILTERNING MASK BASED ON ENTROPY RANKED SITES AND N_SITES DEFINED BY USER
sel = entropy_rank[-n_sites:] #comment out for full matrix 
entropy_mask = np.isin(np.arange(entropy_rank.shape[0]), sel) #Entropy masked. comment out for full matrix


#%% GENERATE K MATRIX 
m = selection.loc[:,:,'m']
r = selection.loc[:,:,'r']
if median_or_mode == 'median':
    K_matrix =  xr.DataArray(beta_distribution.mode(1 + m, 1 + (r - m)), #I changed to median so each version of K matrix has only one probability at each site in
                             dims = ('ct', 'gen'), 
                             coords = {'ct': selection.coords['ct'], # each cell type
                                       'gen':selection.coords['gen']})
    
    sampled_K = K_matrix[:, entropy_mask] #Read depth and entropy ranked site filtered K matrix (sampled K)
    sampled_K.to_netcdf("K_matrices/K_matrix_median_beta_selection/K_matrix_hg38_{}_{}rds_median.nc".format(n_sites, min_rd)) #Sampled K saved
elif median_or_mode == 'mode':
    a = 1+m
    b = 1+(r-m)
    beta_mode = (a-1)/(a+b-2)
    K_matrix =  xr.DataArray(beta_mode, 
                             dims = ('ct', 'gen'), 
                             coords = {'ct': selection.coords['ct'],
                                       'gen':selection.coords['gen']})
    
    sampled_K = K_matrix[:, entropy_mask] 
    sampled_K.to_netcdf("K_matrices/K_matrix_mode_beta_selection/K_matrix_hg38_{}_{}rds_mode.nc".format(n_sites, min_rd)) 
else:
    print('Invalid CpG methylation probability calculation function selected. Please chose "median" or "mode".')




### FOR FULL K MATRIX ###
# selection = xr.open_dataarray('data/rmn/rmn_hg38.nc')
# m = selection.loc[:,:,'m']
# r = selection.loc[:,:,'r']
# a = 1+m
# b = 1+(r-m)
# beta_mode = (a-1)/(a+b-2)
# K_matrix =  xr.DataArray(beta_mode, 
#                             dims = ('ct', 'gen'), 
#                             coords = {'ct': selection.coords['ct'],
#                                     'gen':selection.coords['gen']})
# K_matrix.to_netcdf("K_matrices/K_matrix_hg38_full.nc") 