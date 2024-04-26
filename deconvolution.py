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
def cell_type_deconvolution_LR(cfDNA_data, K): 
    """
    Perform cell type deconvolution using linear regression.

    Parameters:
        cfDNA_data (tuple): Tuple of methylated and not methylated data.
        K (array-like): Coefficient matrix.

    Returns:
        array: Deconvoluted weights.
    """
    data_methylated, data_not_methylated = zip(*cfDNA_data) #added this zip command as wouldn't split into 2 lists otherwise
    read_depths = data_not_methylated #read depth is actually the just the second column of .beta files, so changed read_depths to be equal to that
   
    #  A coefficient matrix
    X = K * read_depths #By performing element-wise multiplication of K and by total reads at each site in each cell type
                        #X becomes the number of expected methylated reads at each site in each cell type, given the observed number of total reads
   
    #  B Ordinate or “dependent variable” values.
    y = data_methylated #Y is then the methylated read data that we actually recorded for in the cfDNA at each site in the methylome
   
    reg = LinearRegression(positive=True).fit(X.T, y) #compute NNLS regression of our X against our Y
    w = reg.coef_  #then extract the cell type proportion coefficents as w 
    return(np.array(w)) #and return them as an array to be used in subsequent steps 

w = cell_type_deconvolution_LR(cfDNA_data=cfDNA_data, K=K)
print(w)
np.savetxt('w_coeffs.txt', w)



#%% COMPUTE NLL USING A BINOMIAL MODEL AND OUR INITIAL CELL TYPE PROPORTION COEFFICIENTS
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




methylated_reads, read_depth = zip(*cfDNA_data)
methylated_reads = np.asarray(methylated_reads) 
read_depth = np.asarray(read_depth) 

r = read_depth 
m = methylated_reads
p = (w[:,np.newaxis] * K).sum(axis=0)
nll = -np.sum(scipy.special.binom(r, m) + m * np.log(p) + (r - m) * np.log(1 - p))

print(p)
print(nll)


#%% GET THETA FOR THE PREVIOUSLY COMPUTED NLL
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
    w = np.exp(theta)
    w /= np.sum(w)
    return negative_log_lihelihood(w, r, m, K)

nll_theta = negative_log_lihelihood_theta(theta=w, r = read_depth, m = methylated_reads, K = K)
print(nll_theta)


#%% MINIMISE THE NLL BY TO FIND MOST ACCURATE ESTIMATION OF CELL TYPE PROPORTION COEFFICIENTS
def cell_type_deconvolution_binomial(cfDNA, K):
    """
    Perform cell type deconvolution using a binomial approach.

    Parameters:
        cfDNA (tuple): Tuple of methylated and not methylated data.
        K (array-like): Coefficient matrix.

    Returns:
        array: Deconvoluted weights.
    """
    methylated, not_methylated = zip(*cfDNA)
    read_depths = np.sum([methylated, not_methylated], axis=0)
    w0 = cell_type_deconvolution_LR(cfDNA, K)
    w0[w0==0] = 10**-100    
    w = np.zeros(K.shape[0])
    theta = np.log(w0)
    res = scipy.optimize.minimize(negative_log_lihelihood_theta, theta, 
                                    args=(read_depths, methylated, K),
                                    method='L-BFGS-B')
    w = np.exp(res.x)
    w /= np.sum(w)

minimise = cell_type_deconvolution_binomial(cfDNA=cfDNA_data, K=K)
print(minimise)















#%% COMPUTE THE KL-DIVERGENCE BETWEEN THE METHYLATED READS AT EACH SITE IN TARGET CELL TYPE (NEURONS) AND ALL OTHER CELL TYPES, THEN RANK BY HIGHEST KL-DIVERGENCE
# We're doing this to hopefully find a selection of sites that are most different between neurons and all other tissues. These sites are likely to be the greatest 
# contributers to the neuronal signal in bulk cfDNA and so we can reduce the sites used in deconvolution to just these sites to improve speed and without impacting
# the accuracy of the deconvolution 
def relative_entropy_rank(rmn, target = 'neuronal'):
    """
    Compute relative entropy rank.

    Parameters:
        rmn (xarray.DataArray): Methylation data.
        target (str): Target cell type.

    Returns:
        array: Ranked sites.
    """
    chance_target_methylated = np.array(rmn.loc[target,:,'m'] / rmn.loc[target,:,'r'])
    chance_target_not_methylated = 1 - chance_target_methylated
    p_target = np.array([chance_target_methylated, chance_target_not_methylated])
    others = list(rmn.coords['ct'].values)
    others.remove(target)
    p_others = rmn.loc[others, :, 'm'].sum(axis=0) / rmn.loc[others, :, 'r'].sum(axis=0)
    p_others = np.array([p_others, 1-p_others])
    kl = entropy(p_target, p_others)
    ranked_sites = np.argsort(kl)
    return ranked_sites



#%% FILTER THE SITES TO BE USED IN DECONVOLUTION BASED ON THE ENTROPY-RANKED SITES GENERATED PREVIOUSLY 
def entropy_filter_n(n, rmn, entropy_rank):
    """
    Filter data based on entropy rank.

    Parameters:
        n (int): Number of sites.
        rmn (xarray.DataArray): Methylation data.
        entropy_rank (array-like): Entropy rank.

    Returns:
        tuple: Filtered data and mask.
    """
    sel = entropy_rank[-n:]
    entropy_mask = np.isin(np.arange(entropy_rank.shape[0]), sel)
    return (rmn[:, entropy_mask, :], entropy_mask)



#%% GENERATE A TEST CFDNA POOL USING A DIRICHLET DISTRIBUTION TO USE IN DECONVOLUTION 
# We can use this test cfDNA pool to assess the performance and accuracy of the deconvolution, as we know the exact proportions of cell types in the
# bulk cfDNA pool 
def generate_test_dirichlet(target_percentage = 1, target_ct_ind = 30):
    """
    Generate test proportions sampling from a Dirichlet distribution.

    Parameters:
        target_percentage (float): Percentage of target cell type.
        target_ct_ind (int): Index of target cell type.

    Returns:
        array: Generated test proportions.
    """
    proportions_test = np.random.dirichlet(np.ones(37), size=1)[0]*(100-target_percentage)
    proportions_test = np.concatenate([proportions_test[:target_ct_ind], 
                                    np.array([target_percentage]), proportions_test[target_ct_ind:]])
    return proportions_test



#%% RUN AN EXPERIMENT WITH MULTIPLE ITERATIONS OVER ENTROPY-RANKED SITES USED, READ DEPTH USED
def experiment(deconvolution_func, synthesiser, target, rmn, ground_truth, min_samples = 0, min_rd_origin = 20, 
                n_individuals = 1, n_iterations = 10, mu_individual_rds = [10], min_k_rds = [20], 
                entropy_n = [5, 10, 100, 1000, 10000, 100000]):
    """
    Perform an experiment to evaluate cell type deconvolution.

    Parameters:
        deconvolution_func (function): Function for cell type deconvolution.
        synthesiser (function): Function for synthesizing data.
        target (str): Target cell type for deconvolution evaluation.
        rmn (xarray.DataArray): Methylation data.
        ground_truth (dict): Dictionary containing ground truth proportions of cell types.
        min_samples (int): Minimum sample count threshold for masking methylation data.
        min_rd_origin (int): Minimum read depth threshold for masking methylation data.
        n_individuals (int): Number of individuals for synthesizing data.
        n_iterations (int): Number of iterations for repeated experiments.
        mu_individual_rds (list): List of mean read depths for synthesizing data.
        min_k_rds (list): List of minimum read depths for generating K matrix.
        entropy_n (list): List of numbers of top-ranked sites to consider for entropy filtering.

    Returns:
        xarray.DataArray: Result of the experiment, containing deconvolution weights for each iteration,
                          individual, group, and entropy configuration.
    """
        
    entropy_n = [int(n) for n in entropy_n]

    groups = list(ground_truth.keys())

    n_groups = len(groups)
    
    muDNAm, selection, mask = generate_muDNAm_mask(rmn, min_rd_origin, min_samples)

    pseudo_K = generate_K(selection)

    individuals_tests = [synthesiser(pseudo_K, proportions, mu_individual_rds[0], n_individuals) for proportions in ground_truth.values()]
 
    entropy_rank = relative_entropy_rank(selection, target)

    parallelised_arguments = []

    for n in entropy_n:
        entropy_filtered, entropy_mask = entropy_filter_n(n, selection, entropy_rank)
        individuals_masked = [[individual[:, entropy_mask]
                               for individual in individuals] for individuals in individuals_tests]

        sampled_Ks = []
        with Pool(n_iterations) as p:
            sampled_Ks = p.map(generate_K, [entropy_filtered]*n_iterations) 

        for group in individuals_masked:
            for individual in group:
                for i in range(n_iterations):
                    parallelised_arguments.append([individual, sampled_Ks[i]])

    with Pool(maxtasksperchild=1) as p:
        flattened_results = p.starmap(deconvolution_func, parallelised_arguments)

    layered_results = 100*np.array(flattened_results).flatten().reshape(len(entropy_n), len(groups), n_individuals, n_iterations, len(rmn.coords['ct']))

    results = xr.DataArray(layered_results, dims=["entropy_n", "group", "individual", "iteration", "ct"],
                          coords={
                                "entropy_n": entropy_n, 
                                "group": groups, 
                                "individual": np.arange(n_individuals), 
                                "iteration": np.arange(n_iterations),
                                "ct": rmn.coords['ct']})
    return(results)
# %%
