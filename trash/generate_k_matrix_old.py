import os
import time
import numpy as np
from scipy.stats import beta as beta_distribution
import xarray as xr

"""
Generate K matrix by sampling from beta distributions.
Parameters:
    rmn (xarray.DataArray): Methylation data.
Returns:
    xarray.DataArray: Generated K matrix.
"""

#IDEAS FOR FUTURE IMPROVEMENTS:
#generate rmn and K matrix scripts together 
#give option in batch submission script to filter by n number of entropy-ranked sites 
rmn = xr.open_dataarray('rmn_hg38_200000_sites.nc')  

epoch_time = int(time.time())
np.random.seed(os.getpid() + epoch_time)

m = rmn.loc[:,:,'m']
r = rmn.loc[:,:,'r']
K_matrix =  xr.DataArray(beta_distribution.median(1 + m, 1 + (r - m)),
                    dims = ('ct', 'gen'), coords = {'ct': rmn.coords['ct'],
                                                    'gen':rmn.coords['gen']})
K_matrix.to_netcdf("K_matrix_hg38_200000_sites.nc")
 