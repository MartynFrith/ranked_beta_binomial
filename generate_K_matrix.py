#%%
import os
import time
import numpy as np
from scipy.stats import beta as beta_distribution
import xarray as xr
from scipy.stats import entropy

#%% OPTIONS
rmn = xr.open_dataarray('data/rmn/rmn_hg38.nc')
min_rd = 10
min_samples = 0
target = 'Neuron'
n = 100000


#%%
#FILTER SITES BY READ DEPTH 
n = rmn.loc[:,:,'n']
ct_mask = np.any(np.greater(n, min_samples), axis=1)
ct_selection = rmn.loc[ct_mask,:]
gen_mask = np.all(np.greater(ct_selection.loc[:,:,'r'], min_rd), axis=0)
gen_selection = ct_selection.loc[:, gen_mask] #This is the xarray that we need to take forward 
muDNAm_masked = np.divide(gen_selection.loc[:,:,'m'], gen_selection.loc[:,:,'r'])
muDNAm_masked = muDNAm_masked.fillna(0)
  
 


#%% RANK SITES BY ENTROPY
selection = gen_selection

chance_target_methylated = np.array(selection.loc[target,:,'m'] / selection.loc[target,:,'r'])
chance_target_not_methylated = 1 - chance_target_methylated
p_target = np.array([chance_target_methylated, chance_target_not_methylated])
others = list(selection.coords['ct'].values)
others.remove(target)
p_others = selection.loc[others, :, 'm'].sum(axis=0) / selection.loc[others, :, 'r'].sum(axis=0)
p_others = np.array([p_others, 1-p_others])
kl = entropy(p_target, p_others)
entropy_rank = np.argsort(kl)




#%% FILTER THE SITES TO BE USED IN DECONVOLUTION BASED ON THE ENTROPY-RANKED SITES GENERATED PREVIOUSLY 
sel = entropy_rank[-n:]
entropy_mask = np.isin(np.arange(entropy_rank.shape[0]), sel)
selection_filtered = selection[:, entropy_mask, :]

#filename = 'rmn_hg38_{}_sites.nc'.format(n)
#selection_filtered.to_netcdf(filename)




#%% GENERATE K MATRIX 

#THE DATA GIVEN (rmn) IS THE OUTPUT OF muDNAM_mask, MEANING THAT INPUT DATA IS ALREADY FILTERED BY READ DEPTH

epoch_time = int(time.time())
np.random.seed(os.getpid() + epoch_time)

m = selection_filtered.loc[:,:,'m']
r = selection_filtered.loc[:,:,'r']
K_matrix =  xr.DataArray(beta_distribution.median(1 + m, 1 + (r - m)),
                    dims = ('ct', 'gen'), coords = {'ct': selection_filtered.coords['ct'],
                                                    'gen':selection_filtered.coords['gen']})

K_matrix.to_netcdf("K_matrices/K_matrix_hg38_{}k_10rds.nc".format(n))
 