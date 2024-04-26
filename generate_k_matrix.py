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

rmn = xr.open_dataarray('log_files/rmn_hg38.nc')  

epoch_time = int(time.time())
np.random.seed(os.getpid() + epoch_time)

m = rmn.loc[:,:,'m']
r = rmn.loc[:,:,'r']
K_matrix =  xr.DataArray(beta_distribution.median(1 + m, 1 + (r - m)),
                    dims = ('ct', 'gen'), coords = {'ct': rmn.coords['ct'],
                                                    'gen':rmn.coords['gen']})
K_matrix.to_netcdf("K_matrix.nc")
