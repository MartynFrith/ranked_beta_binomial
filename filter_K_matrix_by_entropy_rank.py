import xarray as xr
import numpy as np
from scipy.stats import entropy

#%% COMPUTE THE KL-DIVERGENCE BETWEEN THE METHYLATED READS AT EACH SITE IN TARGET CELL TYPE (NEURONS) AND ALL OTHER CELL TYPES, THEN RANK BY HIGHEST KL-DIVERGENCE
# We're doing this to hopefully find a selection of sites that are most different between neurons and all other tissues. These sites are likely to be the greatest 
# contributers to the neuronal signal in bulk cfDNA and so we can reduce the sites used in deconvolution to just these sites to improve speed and without impacting
# the accuracy of the deconvolution 

rmn = xr.open_dataarray('rmn_hg38.nc')
target = 'Neuron'
n = 100000

chance_target_methylated = np.array(rmn.loc[target,:,'m'] / rmn.loc[target,:,'r'])
chance_target_not_methylated = 1 - chance_target_methylated
p_target = np.array([chance_target_methylated, chance_target_not_methylated])
others = list(rmn.coords['ct'].values)
others.remove(target)
p_others = rmn.loc[others, :, 'm'].sum(axis=0) / rmn.loc[others, :, 'r'].sum(axis=0)
p_others = np.array([p_others, 1-p_others])
kl = entropy(p_target, p_others)
ranked_sites = np.argsort(kl)


#%% FILTER THE SITES TO BE USED IN DECONVOLUTION BASED ON THE ENTROPY-RANKED SITES GENERATED PREVIOUSLY 
entropy_rank = ranked_sites
sel = entropy_rank[-n:]
entropy_mask = np.isin(np.arange(entropy_rank.shape[0]), sel)
rmn_filtered = rmn[:, entropy_mask, :]

np.savetxt('rmn_hg38_{}_sites.txt'.format(n), rmn_filtered)
