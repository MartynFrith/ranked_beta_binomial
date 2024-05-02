#%%
import datetime
import xarray as xr
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import scipy.optimize
import scipy.special

#CURRENTLY PROPORTIONS SEEM REALLY OFF. COULD BE DUE TO A) K-MATRIX?, B)CHANGING ALL READ DEPTH VALUES TO JUST TOTAL READ DEPTHS
#WHEN IT MIGHT ACTAULLY NEED TO BE TOTAL READS - METHYLATED READS, C) SYNTHETC CFDNA NOT ACTUALLY BEING THE PROPORTIONS LAID OUT IN THE GROUND TRUTH

#%% LOAD IN CFDNA SEQUENCING DATA IN .BETA FORMAT
K = xr.open_dataarray('K_matrices/K_matrix_hg38_100k_10rds.nc') #load previously created matrix

cfDNA_data = np.fromfile('cfDNA_files/beta_files/11077_control_S1__bwameth_sorted.beta', dtype=np.uint8).reshape((-1, 2)) #load cfDNA methylation sequencing data as numpy array 
cfDNA_data =  tuple(map(tuple, cfDNA_data)) #convert cfDNA methylation data to a tuple to be usable in deconvolutions
cpg_idx = np.asarray(K.gen - 1 ) #zero based indexing so must subtract 1 from each cpg index value in order to filter cfDNA_data properly
cfDNA_data = [cfDNA_data[i] for i in cpg_idx] #index cfDNA_data tuple based on cpg index retrieved from K matrix
methylated_reads, read_depths = zip(*cfDNA_data) 

#OR USE SYNTHETIC CFDNA DATA
# cfDNA_data = np.asarray(pd.read_csv('poisson_synthesised_cfDNA.txt')) # don't need to filter out cpg sites base on cpg indexes in K matrix
# methylated_reads, read_depths = cfDNA_data

#%% PERFORM NNLS REGRESSION TO GET INTIAL CELL TYPE PROPORTION ESIMATES

#  A coefficient matrix
X = K * read_depths #By performing element-wise multiplication of K and by total reads at each site in each cell type
                    #X becomes the number of expected methylated reads at each site in each cell type, given the observed number of total reads

#  B Ordinate or “dependent variable” values.
y = methylated_reads #Y is then the methylated read data that we actually recorded for in the cfDNA at each site in the methylome

reg = LinearRegression(positive=True).fit(X.T, y) #compute NNLS regression of our X against our Y
w0 = reg.coef_  #then extract the cell type proportion coefficents as w 
w0 = np.array(w0) #and return them as an array to be used in subsequent steps 

print('Initial coefficient estimate:', w0)



#%% MINIMISE THE NLL BY TO FIND MOST ACCURATE ESTIMATION OF CELL TYPE PROPORTION COEFFICIENTS
#read_depths = np.sum([methylated_reads, read_depths], axis=0)
methylated_reads = np.asarray(methylated_reads) #PROBLEM AREA
methylated_reads = np.uint64(methylated_reads)
read_depths = np.asarray(read_depths)  #PROBLEM AREA
read_depths = np.uint64(read_depths) #need to make sure floating point of r and m are equal to floating point of K and w in binomial model?

w0[w0==0] = 10**-100    
w_final = np.zeros(K.shape[0])
theta = np.log(w0)


def negative_log_lihelihood_theta(theta, r, m, K):
    w = np.exp(theta) #theta originates from the below minimisation step, as these functions are called there
    w /= np.sum(w)
    p = (w[:,np.newaxis] * K).sum(axis=0) #create array p which is the sum, for each cell type, of cell type proportion coeffecient multiplied by 
                                          #the probability of a methylated read at each site across the methylome
    return -np.sum(scipy.special.binom(r, m) + m * np.log(p) + (r - m) * np.log(1 - p)) #return the NLL with previously defined p being used


res = scipy.optimize.minimize(negative_log_lihelihood_theta, theta, 
                                args=(read_depths, methylated_reads, K),
                                method='L-BFGS-B')


w_final = np.exp(res.x) #Normalise result
w_final /= np.sum(w_final) #Normalise result

celltypes = list(K.ct.values)
w_final = list(zip(celltypes, w_final))

print('Result:', res)
print('Cell type proportion estimate:', w_final)



current_datetime = datetime.datetime.now()
formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H%-M-%S")

filename = f"results/synthetic_cfDNA_100k_10rds_deconvolution_results_{formatted_datetime}.txt"
with open(filename, 'w') as file:
    file.write('NNLS initial coefficient estimation:\n{}\n'.format(w0))
    file.write('Minimisation ouput:\n{}\n'.format(res))
    file.write('Final deconvoluted cell type proportions:\n{}\n'.format(w_final))

#%%
