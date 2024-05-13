#%%
import sys
import xarray as xr
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import os
import time
import datetime
from scipy.stats import beta as beta_distribution
from scipy.stats import entropy
import scipy.optimize
import scipy.special

K_matrix_config_filename = sys.argv[1]
cfDNA_data_filename = sys.argv[2]

# K_matrix_filename = 'K_matrices/K_matrix_hg19_100k_10rds.nc'
# cfDNA_data_filename = 'cfDNA_files/loyfer_cfdna_110033750.beta'

#%% READ CFDNA DATA
if cfDNA_data_filename.endswith('.beta'):
    cfDNA_data = np.fromfile(cfDNA_data_filename, dtype=np.uint8).reshape((-1, 2)) #load cfDNA methylation sequencing data as numpy array 
elif cfDNA_data_filename.endswith('.npy'):             #OR USE SYNTHETIC CFDNA DATA
    cfDNA_data = np.load(cfDNA_data_filename)[0,:,:].T #select only the first 2d numpy array in synthesised cfDNA
else:                                                  #as is a 3d array, 3rd dimension is individuals
    cfDNA_data = np.asarray(pd.read_csv(cfDNA_data_filename)) 


cfDNA_data = pd.DataFrame(cfDNA_data)
cfDNA_data.columns = ['dnam', 'rd']
dnam = pd.DataFrame(cfDNA_data['dnam'])
dnam.columns = ['dnam']
dnanm = pd.DataFrame(cfDNA_data['rd'] - dnam['dnam'])
cfdna = np.array([dnam.values.flatten(), dnanm.values.flatten()])


K_matrix_config = pd.read_csv(K_matrix_config_filename)
mu_rd = round(cfDNA_data['rd'].mean()/10)*10


if mu_rd < 100: #currently highest mean read depth K matrix created for is 100
    K_matrix_filename = K_matrix_config.loc[K_matrix_config['mu_rd'] == mu_rd, 'filename'].values[0]
    sampled_K = xr.open_dataarray(K_matrix_filename)
else:
    K_matrix_filename = K_matrix_config.loc[K_matrix_config['mu_rd'] == 100, 'filename'].values[0]
    sampled_K = xr.open_dataarray(K_matrix_filename)
entropy_mask = sampled_K.gen.values - 1
cfdna_masked = cfdna[:, entropy_mask]


data_methylated, data_not_methylated = cfdna_masked
read_depths = np.sum([data_methylated, data_not_methylated], axis=0)



#%% NNLS REGRESSION 
X = sampled_K * read_depths # X is represents the number of reads we would get at eahc site if there was no missingness in data based on our psuedo K
y = data_methylated #Y is then the methylated read data that we actually recorded for in the cfDNA at each site in the methylome


reg = LinearRegression(positive=True).fit(X.T, y) #compute NNLS regression of our X against our Y
w0 = reg.coef_  #then extract the cell type proportion coefficents as w 
w0 = np.array(w0) #and return them as an array to be used in subsequent steps 


w_print = np.array_str(w0, suppress_small=True)
print('Initial coefficient estimate:', w_print)



#%% BINOMIAL MODEL
data_methylated = np.uint64(data_methylated)
read_depths = np.uint64(read_depths) #need to make sure floating point of r and m are equal to floating point of K and w in binomial model?


w0[w0==0] = 10**-100    
w_final = np.zeros(sampled_K.shape[0])
theta = np.log(w0)


def negative_log_lihelihood_theta(theta, r, m, K):
    w = np.exp(theta) #theta originates from the below minimisation step, as these functions are called there
    w /= np.sum(w)
    p = (w[:,np.newaxis] * K).sum(axis=0) #create array p which is the sum, for each cell type, of cell type proportion coeffecient multiplied by 
                                          #the probability of a methylated read at each site across the methylome
    return -np.sum(scipy.special.binom(r, m) + m * np.log(p) + (r - m) * np.log(1 - p)) #return the NLL with previously defined p being used


res = scipy.optimize.minimize(negative_log_lihelihood_theta, theta, 
                                args=(read_depths, data_methylated, sampled_K),
                                method='L-BFGS-B')
w_final = np.exp(res.x) #Normalise result
w_final /= np.sum(w_final) #Normalise result


celltypes = list(sampled_K.ct.values)
w_final = list(zip(celltypes, w_final))
print('Result:', res)
print('Cell type proportion estimate:', w_final)



#%% SAVE RESULTS 
current_datetime = datetime.datetime.now()
formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H%-M-%S")
K_matrix_filename = K_matrix_filename [K_matrix_filename .find('h') : K_matrix_filename .find('.')]
cfDNA_data_filename = cfDNA_data_filename[cfDNA_data_filename.find('/') + 1 : cfDNA_data_filename.find('.')]


filename = f"results/{cfDNA_data_filename}_{K_matrix_filename }_{formatted_datetime}.txt"
with open(filename, 'w') as file:
    file.write('NNLS initial coefficient estimation:\n{}\n'.format(w_print))
    file.write('Minimisation ouput:\n{}\n'.format(res))
    file.write('Final deconvoluted cell type proportions:\n{}\n'.format(w_final))
#%%
