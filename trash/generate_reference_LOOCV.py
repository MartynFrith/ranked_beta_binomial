# LOAD LIBRARIES AND DATA
import os
import re
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import beta as beta_distribution
import netCDF4

def read_file(betaPath, filename):
    content = np.fromfile(os.path.join(betaPath, filename), dtype=np.uint8).reshape((-1, 2)) #UINT8 MEANS MAX VALUE 255, DO WE WANT THIS? WOULD HIGHER VALUES MAKE A DIFFERENCE?
    dnam = pd.DataFrame({'dnam': np.where(content[:, 1] >= 1,  content[:, 0], 0), 'rd': content[:, 1]})
    return(dnam)


left_out = []
def process_ct(ct):
    betaPath = '{}{}/'.format(source, ct)
    betaFiles = os.listdir(betaPath)
    if len(betaFiles) >= 3:
        leaveoneout = betaFiles[1] #CHANGE FOR DIFFERENT FILE EXCLUDED 
        betaFiles = [file for file in betaFiles if file != leaveoneout]
        betaMat = pd.concat([read_file(betaPath, x) for x in betaFiles], axis = 1)
        rd = pd.DataFrame(betaMat["rd"])
        dnam = pd.DataFrame(betaMat["dnam"])
        rd.columns = np.arange(len(rd.columns))
        dnam.columns = np.arange(len(dnam.columns))
        ct_K = cpgLoci.copy()# sum up over all samples
        n = np.sum(rd != 0, axis=1)
        rd = rd.sum(axis=1)
        dnam = dnam.sum(axis=1)
        ct_K.columns = ['chr', 'chr_ind', 'gen_ind'] #gen / gen_ind is cedric's term for cpg index 
        ct_K['r'] = rd
        ct_K['n'] = n #n is the number of cell specific sequencing files for each cell type 
        ct_K['m'] = dnam
        ct_K.to_csv("{}{}.csv".format(target, ct), index=False)
        left_out.append(leaveoneout)


def process_test_data(excluded_files):
    betaFiles = pd.read_csv(excluded_files, header=None)
    betaFiles = list(betaFiles.iloc[:,0].values)
    for file in betaFiles:
        extract_ct = re.search(r'_(.*?)\-', file)
        ct = extract_ct.group(1)  
        ct = ct.replace('.', '_')    
        current_file = '{}/{}'.format(ct, file)
        betaMat = read_file(source, current_file)
        rd = pd.DataFrame(betaMat["rd"])
        dnam = pd.DataFrame(betaMat["dnam"])
        rd.columns = np.arange(len(rd.columns))
        dnam.columns = np.arange(len(dnam.columns))
        ct_K = cpgLoci.copy()# sum up over all samples
        n = np.sum(rd != 0, axis=1)
        ct_K.columns = ['chr', 'chr_ind', 'gen_ind'] #gen / gen_ind is cedric's term for cpg index 
        ct_K['r'] = rd
        ct_K['n'] = n #n is the number of cell specific sequencing files for each cell type 
        ct_K['m'] = dnam
        ct_K.to_csv("{}{}.csv".format(target, ct), index=False)
        print(current_file)



path = "/lustre/home/mjf221/ranked_beta_binomial/" #path to the project folder 
source = "/lustre/home/mjf221/ranked_beta_binomial/reference_data/betas/" #path to the parent  containing all .beta files
target = "/lustre/home/mjf221/ranked_beta_binomial/reference_data/LOOCV_processed_genome_rmn/" #path to output script to
create_test_data = True 
excluded_files = 'reference_data/LOOCV_rmn/filenames_LOOCV_2ndFiles_removed.csv'
cpgFile = "reference_data/CpG.bed" 

try:
    for file in os.listdir(target):
        os.remove("{}{}".format(target, file))  
    os.rmdir(target) 
except:
    pass 
os.mkdir(target) 

cpgLoci = pd.read_csv(cpgFile, sep = "\t", header = None, names = ["CHR","START","END"], dtype={"CHR": str, "START":int, "END":int})

if create_test_data == True:
        process_test_data(excluded_files)
else: 
    for ct in os.listdir(source):
        print(ct)
        process_ct(ct)
    

#CREATE XARRAY OF SYNTHESISED DATA TABLE
files = os.listdir(target)
rmn_d = {filename.split('.')[0]: pd.read_csv('{}{}'.format(target, filename))for filename in files}
print('Loaded rmn dictionary')

rmn = xr.DataArray(dims=("ct", "gen", "v"), coords={"ct" : list(rmn_d.keys()), "gen" : rmn_d[list(rmn_d.keys())[0]]['gen_ind'].values, 
                                                    "v" : ["r", "m", "n"]})

for ct in list(rmn_d.keys()):
    print(ct)
    rmn.loc[ct] = rmn_d[ct].loc[:, ('r', 'm', 'n')]
rmn.to_netcdf('reference_data/LOOCV_rmn/TestData_hg38_LOOCV_2ndFiles_removed.nc', engine='netcdf4') 

if create_test_data == False:
    np.savetxt("reference_data/LOOCV_rmn/filenames_LOOCV_2ndFiles_removed.csv", left_out, fmt='%s', delimiter = ',')
print('Generated rmn.nc')

m = rmn.loc[:,:,'m']
r = rmn.loc[:,:,'r']
k_matrix =  xr.DataArray(beta_distribution.median(1 + m, 1 + (r - m)), 
                        dims = ('ct', 'gen'), 
                        coords = {'ct': rmn.coords['ct'], # each cell type
                                    'gen':rmn.coords['gen']})
k_matrix_filename = 'reference_k_matrix_LOOCV_1stFiles_removed'
k_matrix.to_netcdf(f'reference_data/K_matrices/{k_matrix_filename}')
print('Generated K matrix')