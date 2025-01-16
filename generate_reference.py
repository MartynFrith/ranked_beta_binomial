import os
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import beta as beta_distribution

def read_file(betaPath, filename):
    content = np.fromfile(os.path.join(betaPath, filename), dtype=np.uint8).reshape((-1, 2)) 
    dnam = pd.DataFrame({'dnam': np.where(content[:, 1] >= 1,  content[:, 0], 0), 'rd': content[:, 1]})
    return(dnam)
    

def process_ct(ct):
    betaPath = '{}{}/'.format(source, ct)
    betaFiles = os.listdir(betaPath)
    betaMat = pd.concat([read_file(betaPath, x) for x in betaFiles], axis = 1)
    rd = pd.DataFrame(betaMat["rd"])
    dnam = pd.DataFrame(betaMat["dnam"])
    rd.columns = np.arange(len(rd.columns))
    dnam.columns = np.arange(len(dnam.columns))
    ct_K = cpgLoci.copy()
    n = np.sum(rd != 0, axis=1)
    rd = rd.sum(axis=1)
    dnam = dnam.sum(axis=1)
    ct_K.columns = ['chr', 'chr_ind', 'gen_ind'] 
    ct_K['r'] = rd
    ct_K['n'] = n 
    ct_K['m'] = dnam
    ct_K.to_csv("{}{}.csv".format(target, ct), index=False)



# USER ARGUMENTS
source = "/lustre/home/mjf221/ranked_beta_binomial/reference_data/betas/" #path to the parent  containing all .beta files
target = "/lustre/home/mjf221/ranked_beta_binomial/reference_data/processed_genome_rnm/" #path to output script to
cpgFile = "reference_data/CpG.bed" 
rmn_filename = "rmn_hg38"
k_matrix_filename = "k_matrix_hg38"


# CHECK DIRECTORIES 
try:
    for file in os.listdir(target):
        os.remove("{}{}".format(target, file)) 
    os.rmdir(target) 
except:
    pass 
os.mkdir(target)


# LOAD CPG REFERENCE AND CELL TYPE SPECIFIC SEQUENCING DATA 
cpgLoci = pd.read_csv(cpgFile, sep = "\t", header = None, names = ["CHR","START","END"], dtype={"CHR": str, "START":int, "END":int})                                                                                             
for ct in os.listdir(source):
    print(ct)
    process_ct(ct)
files = os.listdir(target)
rmn_d = {filename.split('.')[0]: pd.read_csv('{}{}'.format(target, filename))for filename in files}


# CREATE RMN (SINGLE XARRAY CONTAINING ALL CELL TYPE SEQUENCING DATA)
rmn = xr.DataArray(dims=("ct", "gen", "v"), coords={"ct" : list(rmn_d.keys()), "gen" : rmn_d[list(rmn_d.keys())[0]]['gen_ind'].values, 
                                                    "v" : ["r", "m", "n"]})
for ct in list(rmn_d.keys()):
    print(ct)
    rmn.loc[ct] = rmn_d[ct].loc[:, ('r', 'm', 'n')]
rmn.to_netcdf(f'{rmn_filename}.nc') 


# CREATE K MATRIX (REFERENCE MATRIX)
m = rmn.loc[:,:,'m']
r = rmn.loc[:,:,'r']
k_matrix =  xr.DataArray(beta_distribution.median(1 + m, 1 + (r - m)), 
                        dims = ('ct', 'gen'), 
                        coords = {'ct': rmn.coords['ct'], 
                                    'gen':rmn.coords['gen']})
k_matrix.to_netcdf(f'{k_matrix_filename}.nc')


