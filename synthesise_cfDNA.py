#%%
import numpy as np
import pandas as pd
import xarray as xr
import sys
import re

synthesiser = sys.argv[1]
K_matrix = sys.argv[2]
cell_types_proportions = sys.argv[3]
mu_rd = float(sys.argv[4])
n_individuals = int(sys.argv[5])

#%% CREATE DIRICHLET DISTRIBUTION OF CELL TYPE PROPORTIONS TO USE AS GROUND TRUTH 
def generate_test_dirichlet(target_percentage, target_ct_ind):
    """
    Generate test proportions sampling from a Dirichlet distribution.

    Parameters:
        target_percentage (float): Percentage of target cell type.
        target_ct_ind (int): Index of target cell type.

    Returns:
        array: Generated test proportions.
    """
    proportions_test = np.random.dirichlet(np.ones(38), size=1)[0]*(100-target_percentage) #changed np.ones from 37 to 38, as my K matrices have 39 cell types total
    proportions_test = np.concatenate([proportions_test[:target_ct_ind], 
                                    np.array([target_percentage]), proportions_test[target_ct_ind:]])
    return proportions_test




#%% CREATE GROUND TRUTHS IF DONT HAVE THEM ALREADY
def generate_ground_truth(target, rmn_xr, mu, sites):
    target_ct_ind = list(rmn_xr.coords['ct'].values).index(target)
    ground_truth = {}

    base_ct_list = list(rmn_xr.coords['ct'].values)

    gt_df = pd.read_csv('ground_truth_0.csv')

    gt_proportions = gt_df.set_index('ct').to_dict()['proportion']

    base = np.array([gt_proportions[ct] for ct in base_ct_list]) #needed if we you want to make multiple ground truths, see experiment scripts

    ground_truth = generate_test_dirichlet(1, target_ct_ind)

    ground_truth_list = list(ground_truth)
    total = sum(ground_truth_list) #not sure why but the original code did not resize the ground truth proportions to total to 1
    ground_truth_list = [i / total for i in ground_truth_list] #the two lines here do just that

    # dgt = xr.DataArray(list(ground_truth.values()), dims=["group", "ct"], 
    #                     coords={ "group": list(ground_truth.keys()), 
    #                             "ct": rmn_xr.coords['ct']})

    # dgt.to_netcdf("ground_truth_%i_%i.nc"%(mu, sites))




#%% SYNTHESIS CFDNA TO BE USED IN DECONVOLUTION TESTING 
#Need to use the a K matrix for muDNAm, gives syntheised cfDNA data the same number of sites
def synthesise_poisson(muDNAm_, cell_types_proportions = [], mu_rd = 10, n_individuals = 1):
    """
    Synthesize Cell Free DNA data following a Poisson distribution based on given methylation data and parameters.

    Parameters:
        muDNAm_ (xarray.DataArray): Methylation data.
        cell_types_proportions (list): Proportions of cell types.
        mu_rd (int): Mean read depth.
        n_individuals (int): Number of cfDNA samples you want to synthesise.

    Returns:
        list: Synthesized Poisson data.
    """
    n_sites = muDNAm_.shape[1]
    weightedMean = sum([muDNAm_.loc[ct]*cell_types_proportions[i]/100 
                        for i, ct in enumerate(muDNAm_.coords['ct'].values)])
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
        percent_zeros = 100*(synData[1,][synData[1,]==0].shape[0] / n_sites)         #so changed this to what it is currently 
    return synthesised, percent_zeros




#%%
def synthesise_zero_inflated_poisson(muDNAm_, cell_types_proportions = [], mu_rd = 10, n_individuals = 1):
    """
    Synthesize  Cell Free DNA data following a Zero inflated Poisson distribution based on given methylation data and parameters.

    Parameters:
        muDNAm_ (xarray.DataArray): Methylation data.
        cell_types_proportions_df (list): Proportions of cell types.
        mu_rd (int): Mean read depth.
        n_individuals (int): Number of individuals.

    Returns:
        list: Synthesized Poisson data with 5% zero counts.
    """

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



muDNAm_ = xr.open_dataarray(K_matrix) 
cell_types_proportions_df = pd.read_csv(cell_types_proportions)
cell_types_proportions_list = list(cell_types_proportions_df['proportion'])

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