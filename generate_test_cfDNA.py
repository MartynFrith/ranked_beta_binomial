#%%
import numpy as np
import pandas as pd
import xarray as xr
import sys
import re
from scipy.stats import beta as beta_distribution


def generate_k_matrix(rmn):    
    m = rmn.loc[:,:,'m']
    r = rmn.loc[:,:,'r']
    k_matrix =  xr.DataArray(beta_distribution.median(1 + m, 1 + (r - m)), 
                            dims = ('ct', 'gen'), 
                            coords = {'ct': rmn.coords['ct'], 
                                        'gen':rmn.coords['gen']})
    return k_matrix


def synthesise_poisson(k_matrix, cell_types_proportions_list = [], mu_rd = 10, n_individuals = 1):
    n_sites = k_matrix.shape[1]
    weightedMean = sum([k_matrix.loc[ct]*cell_types_proportions_list[i]/100 
                        for i, ct in enumerate(k_matrix.coords['ct'].values)])
    nan_indices = np.where(np.isnan(weightedMean.values))
    weightedMean[nan_indices] = 0
    synthesised = []

    for i in np.arange(n_individuals):
        sample_rd = np.random.poisson(lam=mu_rd, size=n_sites)
        dnam = np.random.binomial(n=sample_rd, p=weightedMean)
        synData = np.array([dnam, sample_rd]) #used to be np.array([dnam, sample_rd - dnam]) but this would produce
                                                #an array with meth and unmeth reads, not an array with meth and read depth
        synData[1,nan_indices]=0
        synthesised.append(synData)  
        percent_zeros = 100*(synData[1,][synData[1,]==0].shape[0] / n_sites)
    return synthesised, percent_zeros


def synthesise_zero_inflated_poisson(muDNAm_, cell_types_proportions = [], mu_rd = 10, n_individuals = 1):
    pi = 5/100 #additional pecercentage of zeros to add

    n_sites = muDNAm_.shape[1]
    weightedMean = sum([muDNAm_.loc[ct]*cell_types_proportions[i]/100 
                    for i, ct in enumerate(muDNAm_.coords['ct'].values)])
    ## set all sites that are NaN in weighted mean to 0 and save their indices
    nan_indices = np.where(np.isnan(weightedMean.values))
    weightedMean[nan_indices] = 0
    synthesised = []

    for i in np.arange(n_individuals):
        bernoulli = np.random.binomial(1, 1-pi, n_sites)
        poisson = np.random.poisson(lam=mu_rd, size=n_sites)
        sample_rd = np.multiply(bernoulli, poisson)
        dnam = np.random.binomial(n=sample_rd, p=weightedMean)
        synData = np.array([dnam, sample_rd])
        synData[1,nan_indices]=0 # by setting the second column of syndata to 0 at each site featuring a NaN in the K matrix we make it so these sites
        synthesised.append(synData) # are represented as 0 total reads and 0 methylated reads in the synthetic cfDNA, reflected the fact these sites were missing 
        percent_zeros = 100*(synData[1,][synData[1,]==0].shape[0] / n_sites)# in the sequencing data that makes up the reference matrix
    return synthesised, percent_zeros             



cell_types_proportions_df = pd.read_csv('ground_truths/emperical_ground_truth_5percent_neuron_hg38.csv')
cell_types_proportions_list = list(cell_types_proportions_df['proportion'])
rmn = xr.open_dataarray('reference_data/LOOCV_rmn/TestData_hg38_LOOCV_1stFiles_removed.nc')
synthesiser = 'synthesise_poisson'





if synthesiser == 'synthesise_poisson' :
    synthesised_cfDNA = synthesise_poisson(muDNAm_= muDNAm_, cell_types_proportions=cell_types_proportions_list, mu_rd=mu_rd, n_individuals=n_individuals)
    synthesiser_for_file = 'Pois'
elif synthesiser == 'synthesise_zero_inflated_poisson':
    synthesised_cfDNA = synthesise_zero_inflated_poisson(muDNAm_= muDNAm_, cell_types_proportions=cell_types_proportions_list, mu_rd=mu_rd, n_individuals=n_individuals)
    synthesiser_for_file = 'zeroInflatedPois'
else:
    print('Invalid synthesiser. Please choose either "synthesise_poisson" or "synthesise_zero_inflated_poisson".')

hg19 = 'hg19'
hg38 = 'hg38'
if hg19 in str(K_matrix):
    alignment = hg19
elif hg38 in K_matrix:
    alignment = hg38
else:
    alignment = ''

pattern = r'(\d+(\.\d+)?)percent'
match = re.search(pattern, sys.argv[3])
target_ct_prop = match.group(1)
target_ct_prop = re.sub('\\.', '' , target_ct_prop)
percent_zeros = round(synthesised_cfDNA[1]) # type: ignore

for i in range(len(synthesised_cfDNA[0])):
    df = synthesised_cfDNA[0][i]
    np.save('S{}_{}_cfDNA_{}_{}%neu_{}rds_{}%zeros.npy'.format(i, synthesiser_for_file, alignment, target_ct_prop, str(int(mu_rd)), percent_zeros), df)