import xarray as xr
from scipy.stats import beta as beta_distribution

rmn = xr.open_dataarray('reference_data/LOOCV_rmn/rmn_hg38_LOOCV_3rdFiles_removed.nc')
m = rmn.loc[:,:,'m']
r = rmn.loc[:,:,'r']
k_matrix =  xr.DataArray(beta_distribution.median(1 + m, 1 + (r - m)), 
                        dims = ('ct', 'gen'), 
                        coords = {'ct': rmn.coords['ct'],
                                    'gen':rmn.coords['gen']})
k_matrix.to_netcdf('reference_k_matrix_LOOCV_3rdFiles_removed.nc')