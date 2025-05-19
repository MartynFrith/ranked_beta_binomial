import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import beta as beta_distribution
from scipy.spatial.distance import jensenshannon


def read_file(filename):
    content = np.fromfile(filename, dtype=np.uint8).reshape((-1, 2)) 
    dnam = pd.DataFrame({'dnam': np.where(content[:, 1] >= 1,  content[:, 0], 0), 'rd': content[:, 1]})
    return(dnam)

def generate_jsRanked_kMatrix(rmn, target_ct):
    n = rmn.loc[:,:,'n']
    ct_mask = np.any(np.greater(n, ct_min_samples), axis=1)
    selection = rmn.loc[ct_mask,:]

    chance_target_methylated = np.array(selection.loc[target_ct,:,'m'] / selection.loc[target_ct,:,'r']) #0/0 here results in nans when working out kl divergence
    chance_target_not_methylated = 1 - chance_target_methylated
    p_target = np.array([chance_target_methylated, chance_target_not_methylated])
    others = list(selection.coords['ct'].values)
    others.remove(target_ct)
    p_others = selection.loc[others, :, 'm'].sum(axis=0) / selection.loc[others, :, 'r'].sum(axis=0)
    p_others = np.array([p_others, 1-p_others])
    js = jensenshannon(p_target, p_others, base=2)

    nan_mask = np.isnan(js)
    js_nans_replaced = xr.where(nan_mask, 0.0, js) # NaN values must be replaced for the k matrix to be functional
    sites_ranked_by_entropy = np.argsort(js_nans_replaced) 
    sel = sites_ranked_by_entropy

    m = selection.loc[:,:,'m']
    r = selection.loc[:,:,'r']
    k_matrix =  xr.DataArray(beta_distribution.median(1 + m, 1 + (r - m)), 
                                dims = ('ct', 'gen'), 
                                coords = {'ct': selection.coords['ct'], 
                                        'gen':selection.coords['gen']})
    nan_mask = np.isnan(k_matrix)
    k_matrix_nans_replaced = xr.where(nan_mask, 0.5, k_matrix) # for some reason if a and b params are both then nan is produced. Replace these instances with 0.5

    sampled_K = k_matrix_nans_replaced.isel(gen=sel[::-1])
    return sampled_K.to_netcdf("kMatrix_JSRanked_20ct_hg38.nc")


# USER ARGUMENTS
rmn_path = 'reference_data/rmn/rmn_20ct_hg38.nc'
target_ct = 'Neuron'
ct_min_samples = 0

# GENERATE K MATRIX
rmn = xr.open_dataarray(rmn_path)
generate_jsRanked_kMatrix(rmn, target_ct)

