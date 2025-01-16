import os
import pandas as pd
import xarray as xr
from scipy.stats import beta as beta_distribution
import multiprocessing as mp
import re

def generate_rmn(data_set):
    target = data_set['folder']
    target_filename =  re.search(r"[^/]+(?=/$)", target)[0]

    files = os.listdir(target)
    rmn_d = {filename.split('.')[0]: pd.read_csv('{}{}'.format(target, filename))for filename in files}
    rmn = xr.DataArray(dims=("ct", "gen", "v"), coords={"ct" : list(rmn_d.keys()), 
                                                        "gen" : rmn_d[list(rmn_d.keys())[0]]['gen_ind'].values, 
                                                        "v" : ["r", "m", "n"]})
    
    for ct in list(rmn_d.keys()):
        print(ct)
        rmn.loc[ct] = rmn_d[ct].loc[:, ('r', 'm', 'n')]

    new_ct_list = list()
    for i in rmn_d.keys():
        new_ct = re.search(r"[A-Z].*$", i)[0]
        new_ct_list.append(new_ct)

    rmn = rmn.assign_coords(ct=new_ct_list)
    rmn.to_netcdf(f'testing/trainingData_rmn/{target_filename}_rmn_hg38.nc', engine='netcdf4') 
    print('Generated rmn.nc')
    
    m = rmn.loc[:,:,'m']
    r = rmn.loc[:,:,'r']
    k_matrix =  xr.DataArray(beta_distribution.median(1 + m, 1 + (r - m)), 
                            dims = ('ct', 'gen'), 
                            coords = {'ct': rmn.coords['ct'], 
                                        'gen':rmn.coords['gen']})
    k_matrix.to_netcdf(f'testing/trainingData_kMatrices/{target_filename}_kMatrix_hg38.nc') 
    print('Generated k_matrix_hg38.nc')


## RUN ##

config = pd.read_csv("/lustre/home/mjf221/ranked_beta_binomial/testing/trainingData_rmn/config.csv")
data_set = config.to_dict(orient='records')

#with mp.Pool(mp.cpu_count()) as pool:
#    all_rmn = pool.map(generate_rmn, data_set)



#################################################################################################################
target = config.iloc[0,0]
target_filename =  re.search(r"[^/]+(?=/$)", target)[0]

files = os.listdir(target)
rmn_d = {filename.split('.')[0]: pd.read_csv('{}{}'.format(target, filename))for filename in files}
rmn = xr.DataArray(dims=("ct", "gen", "v"), coords={"ct" : list(rmn_d.keys()), 
                                                    "gen" : rmn_d[list(rmn_d.keys())[0]]['gen_ind'].values, 
                                                    "v" : ["r", "m", "n"]})

for ct in list(rmn_d.keys()):
    print(ct)
    rmn.loc[ct] = rmn_d[ct].loc[:, ('r', 'm', 'n')]

new_ct_list = list()
for i in rmn_d.keys():
    new_ct = re.search(r"[A-Z].*$", i)[0]
    new_ct_list.append(new_ct)

rmn = rmn.assign_coords(ct=new_ct_list)
rmn.to_netcdf(f'testing/trainingData_rmn/{target_filename}_rmn_hg38.nc', engine='netcdf4') 