#%%
import numpy as np
import pandas as pd
import xarray as xr

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
target = 'Neuron'
rmn_xr = xr.open_dataarray('data/rmn/rmn_hg38.nc')
mu = 10
sites = 100000

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
#Need to use the a K matrix for muDNAm, gives syntheised cfDNA data the same numbe of sites

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
    weightedMean = sum([muDNAm_.loc[ct]*cell_types_proportions[i] 
                        for i, ct in enumerate(muDNAm_.coords['ct'].values)])
    synthesised = []
    for i in np.arange(n_individuals):
        sample_rd = np.random.poisson(lam=mu_rd, size=n_sites)
        dnam = np.random.binomial(n=sample_rd, p=weightedMean)
        synData = np.array([dnam, sample_rd]) #used to be np.array([dnam, sample_rd - dnam]) but this would produce
                                              #an array with meth and unmeth reads, not an array with meth and read depth
        synthesised.append(synData)           #so changed this to what it is currently 
    return synthesised



K_matrix = xr.open_dataarray('K_matrices/K_matrix_hg38_100k_10rds.nc') #data used called rmn in the original syntheise_cfDNA functions, but in the experiment 
                                                                  #function the previosuly generated K matrix is used as the input

muDNAm_ =  xr.open_dataarray('data/rmn/rmn_hg38.nc')
cell_types_proportions = pd.read_csv('ground_truths/emperical_ground_truth.csv')
cell_types_proportions = list(cell_types_proportions['proportion'])
mu_rd = 10
n_individuals = 1

synthesised_cfDNA = synthesise_poisson(muDNAm_= K_matrix, cell_types_proportions=cell_types_proportions, mu_rd=10, n_individuals=1)
synthesised_cfDNA = pd.DataFrame(np.squeeze(synthesised_cfDNA))
synthesised_cfDNA.to_csv('poisson_synthesised_cfDNA.txt', index=False)


# %%
