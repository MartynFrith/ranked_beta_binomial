#%%
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.stats import beta as beta_distribution
from collections import defaultdict
from functools import partial
from math import comb

from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error, explained_variance_score
from sklearn.neighbors import KernelDensity

import scipy.optimize
import scipy.special
from scipy.stats import entropy

from multiprocessing import Pool
import os
import time



#%%
K = xr.open_dataarray('K_matrix.nc') #load previously created matrix
cfDNA_data = np.fromfile('cfDNA_beta/11077_cfdnaB_S8__bwameth_sorted.beta', dtype=np.uint8).reshape((-1, 2)) #load cfDNA methylation sequencing data as numpy array 
cfDNA_data =  tuple(map(tuple, cfDNA_data)) #convert cfDNA methylation data to a tuple to be usable in deconvolutions



#%% PERFORM NNLS REGRESSION TO GET INTIAL CELL TYPE PROPORTION ESIMATES

methylated_reads, read_depths = zip(*cfDNA_data) #added this zip command as wouldn't split into 2 lists otherwise

#  A coefficient matrix
X = K * read_depths #By performing element-wise multiplication of K and by total reads at each site in each cell type
                    #X becomes the number of expected methylated reads at each site in each cell type, given the observed number of total reads

#  B Ordinate or “dependent variable” values.
y = methylated_reads #Y is then the methylated read data that we actually recorded for in the cfDNA at each site in the methylome

reg = LinearRegression(positive=True).fit(X.T, y) #compute NNLS regression of our X against our Y
w = reg.coef_  #then extract the cell type proportion coefficents as w 
w = np.array(w) #and return them as an array to be used in subsequent steps 

print('Initial coefficient estimate:', w)



#%% COMPUTE NLL USING A BINOMIAL MODEL AND OUR INITIAL CELL TYPE PROPORTION COEFFICIENTS
#r = np.asarray(methylated_reads) 
#m = np.asarray(read_depths) 

#p = (w[:,np.newaxis] * K).sum(axis=0) #create array p which is the sum, for each cell type, of cell type proportion coeffecient multiplied by 
                                        #the probability of a methylated read at each site across the methylome
#nll = -np.sum(scipy.special.binom(r, m) + m * np.log(p) + (r - m) * np.log(1 - p)) #return the NLL with previously defined p being used

#print('p:', p)
#print('NLL:', nll)


#%% FUNCTION FOR COMPUTING NLL USING BINOMIAL MODEL
def negative_log_lihelihood(w, r, m, K):
    """
    Compute negative log-likelihood.
    Parameters:
        w (array-like): Cell type proportion coeffecients.
        r (array-like): Read depths.
        m (array-like): Methylated counts.
        K (array-like): Probability matrix.
    Returns:
        float: Negative log-likelihood.
    """
    p = (w[:,np.newaxis] * K).sum(axis=0) #create array p which is the sum, for each cell type, of cell type proportion coeffecient multiplied by 
                                          #the probability of a methylated read at each site across the methylome
    return -np.sum(scipy.special.binom(r, m) + m * np.log(p) + (r - m) * np.log(1 - p)) #return the NLL with previously defined p being used
#%% FUNCTION FOR COMPUTING THETA FOR THE PREVIOUSLY COMPUTED NLL
def negative_log_lihelihood_theta(theta, r, m, K):
    """
    Compute negative log-likelihood with theta.
    Parameters:
        theta (array-like): Theta values.
        r (array-like): Read depths.
        m (array-like): Methylated counts.
        K (array-like): Coefficient matrix.
    Returns:
        float: Negative log-likelihood.
    """
    w = np.exp(theta) #theta originates from the below minimisation step, as these functions are called there
    w /= np.sum(w)
    return negative_log_lihelihood(w, r, m, K)



#%% MINIMISE THE NLL BY TO FIND MOST ACCURATE ESTIMATION OF CELL TYPE PROPORTION COEFFICIENTS

read_depths = np.sum([methylated_reads, read_depths], axis=0)
w0 = w
w0[w0==0] = 10**-100    
w_final = np.zeros(K.shape[0])
theta = np.log(w0)
res = scipy.optimize.minimize(negative_log_lihelihood_theta, theta, 
                                args=(read_depths, methylated_reads, K),
                                method='L-BFGS-B')
w_final = np.exp(res.x)
w_final /= np.sum(w_final) #Normalise result
#w_final[w_final == 1.33876992e-100] = 0 

print('Result:', res)
print('Cell type proportion estimate:', w_final)

with open('deconvolution_results.txt', 'w') as file:
    file.write('NNLS initial coefficient estimation: {}\n'.format(w0))
    file.write('Minimisation ouput: {}\n'.format(res))
    file.write('Final deconvoluted cell type proportions: {}\n'.format(w_final))
