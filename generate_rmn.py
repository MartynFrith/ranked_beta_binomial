#%% LOAD LIBRARIES AND DATA
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import netCDF4

path = "/lustre/home/mjf221/ranked_beta_binomial/" #path to the project folder 
source = "/lustre/home/mjf221/ranked_beta_binomial/reference_data/betas/" #path to the parent  containing all .beta files
target = "/lustre/home/mjf221/ranked_beta_binomial/reference_data/processed_genome_rnm/" #path to output script to


try:
    for file in os.listdir(target):
        os.remove("{}{}".format(target, file)) #this block goes through the target directory and removes any files in there 
    os.rmdir(target) #then attempts to remove the directory itself
except:
    pass #but does nothing if any errors are raised (except = if any exceptions (meaning errors), pass = do nothing )

os.mkdir(target) # and then either way makes a new directory using the path given by target

## load CpG location info
cpgFile = "reference_data/CpG.bed" #Cedric originally used hg19 but have changed to hg38 
cpgLoci = pd.read_csv(cpgFile, sep = "\t", header = None, names = ["CHR","START","END"], dtype={"CHR": str, "START":int, "END":int})          #to match hg38.beta files that I am using. Difference 
                                                                                                       #is ~1 million more CpGs in hg38 compared hg19



#%%  CREATE FORMATTED .CSV OF BETAS FOR EACH CELL TYPE


def read_file(betaPath, filename):
    content = np.fromfile(os.path.join(betaPath, filename), dtype=np.uint8).reshape((-1, 2)) #UINT8 MEANS MAX VALUE 255, DO WE WANT THIS? WOULD HIGHER VALUES MAKE A DIFFERENCE?
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
    ct_K = cpgLoci.copy()# sum up over all samples
    n = np.sum(rd != 0, axis=1)
    rd = rd.sum(axis=1)
    dnam = dnam.sum(axis=1)
    ct_K.columns = ['chr', 'chr_ind', 'gen_ind'] #gen / gen_ind is cedric's term for cpg index 
    ct_K['r'] = rd
    ct_K['n'] = n #n is the number of cell specific sequencing files for each cell type 
    ct_K['m'] = dnam
    ct_K.to_csv("{}{}.csv".format(target, ct), index=False)



for ct in os.listdir(source):
    print(ct)
    process_ct(ct)


#%% CREATE XARRAY OF SYNTHESISED DATA TABLE
files = os.listdir(target)


rmn_d = {filename.split('.')[0]: pd.read_csv('{}{}'.format(target, filename))for filename in files}
print('Loaded rmn dictionary')


rmn = xr.DataArray(dims=("ct", "gen", "v"), coords={"ct" : list(rmn_d.keys()), "gen" : rmn_d[list(rmn_d.keys())[0]]['gen_ind'].values, 
                                                    "v" : ["r", "m", "n"]})

for ct in list(rmn_d.keys()):
    print(ct)
    rmn.loc[ct] = rmn_d[ct].loc[:, ('r', 'm', 'n')]
rmn.to_netcdf('rmn_hg38.nc') 
print('Generated rmn.nc')
# %%
