import pandas as pd 
import numpy as np
import scipy
from scipy.stats import beta as beta_distribution
from sklearn.linear_model import LinearRegression
import xarray as xr
from scipy.stats import entropy
import multiprocessing as mp

def read_cfdna(current_file, min_rd):
    if current_file.endswith('.beta'):
        cfDNA_data = np.fromfile(current_file, dtype=np.uint8).reshape((-1, 2))#load cfDNA methylation sequencing data as numpy array 
    elif current_file.endswith('.lbeta'):
        cfDNA_data = np.fromfile(current_file, dtype=np.uint16).reshape((-1, 2))
    elif current_file.endswith('.npy'):             #OR USE SYNTHETIC CFDNA DATA
        cfDNA_data = np.load(current_file).T        #select only the first 2d numpy array in synthesised cfDNA
    else:                                                  #as is a 3d array, 3rd dimension is individuals
        cfDNA_data = np.asarray(pd.read_csv(current_file)) 
    cfDNA_data = pd.DataFrame(cfDNA_data)
    cfDNA_data.columns = ['dnam', 'rd']
    cfDNA_data['cpg_idx'] = cfDNA_data.index 
    high_reads_idx = cfDNA_data[np.greater(cfDNA_data.loc[:, 'rd'], min_rd)==True] #get locs with 5 or more reads 
    high_reads_idx = np.array(high_reads_idx['cpg_idx']) 
    return cfDNA_data, high_reads_idx


def generate_k_matrix(rmn, target, min_samples, n_sites, high_reads_idx):
    selection = rmn[:,high_reads_idx,:] #remove sites with less than minimum read depth in bulk sequencing data, need to minus -1 due to python 0-based indexing but cpg_idx is 1-based
    #n = rmn.loc[:,:,'n'] 
    #ct_mask = np.any(np.greater(n, min_samples), axis=1) # find cell types with less than the minimum number of cell type specific sequencing samples from reference matrix
    #selection = rmn.loc[ct_mask,:] # remove the above cell types
    chance_target_methylated = np.array(selection.loc[target,:,'m'] / selection.loc[target,:,'r'])
    chance_target_not_methylated = 1 - chance_target_methylated
    p_target = np.array([chance_target_methylated, chance_target_not_methylated])
    others = list(selection.coords['ct'].values)
    others.remove(target)
    p_others = selection.loc[others, :, 'm'].sum(axis=0) / selection.loc[others, :, 'r'].sum(axis=0) #takes long time with many sites, might need batch job
    p_others = np.array([p_others, 1-p_others])
    kl = entropy(p_target, p_others)  
    kl_data_frame = pd.DataFrame({'cpg_idx' : selection.gen.values})
    kl_data_frame['entropy'] =kl 
    entropy_rank = kl_data_frame.sort_values(by='entropy', ascending=False) #Entropy ranked sites
    top_ranked_sites = entropy_rank[0:n_sites] 
    top_ranked_sites_cpg_idx = np.asarray(top_ranked_sites['cpg_idx'])
    sampled_k = selection.sel(gen=top_ranked_sites_cpg_idx)
    m = sampled_k.loc[:,:,'m']
    r = sampled_k.loc[:,:,'r']
    k_matrix =  xr.DataArray(beta_distribution.median(1 + m, 1 + (r - m)), 
                            dims = ('ct', 'gen'), 
                            coords = {'ct': sampled_k.coords['ct'], # each cell type
                                        'gen':sampled_k.coords['gen']})
    return k_matrix


def deconvolute(cfDNA_data, k_matrix):
    cfdna_entropy_filtered = cfDNA_data.loc[k_matrix.gen.values-1, :]
    data_methylated= np.array(cfdna_entropy_filtered.loc[:,'dnam'])
    data_not_methylated = np.array(cfdna_entropy_filtered.loc[:,'rd'] - cfdna_entropy_filtered.loc[:,'dnam'])
    if np.any(np.isnan(k_matrix)==True):
        gen_with_nans = k_matrix.gen.where(k_matrix.isnull().any(dim="ct"), drop=True)
        dropped_gen_values = gen_with_nans.values
        k_matrix = k_matrix.dropna(dim="gen", how="any")

        data_methylated= np.array(cfdna_entropy_filtered[~cfdna_entropy_filtered['cpg_idx'].isin(dropped_gen_values-1)]['dnam'])
        read_depths = np.array(cfdna_entropy_filtered[~cfdna_entropy_filtered['cpg_idx'].isin(dropped_gen_values-1)]['rd'])
    else:
        data_methylated= np.array(cfdna_entropy_filtered.loc[:,'dnam'])
        read_depths = np.sum([data_methylated, data_not_methylated], axis=0) # NNLS REGRESSION 
    X = k_matrix * read_depths # X is represents the number of reads we would get at eahc site if there was no missingness in data based on our psuedo K
    y = data_methylated #Y is then the methylated read data that we actually recorded for in the cfDNA at each site in the methylome
    reg = LinearRegression(positive=True).fit(X.T, y) #compute NNLS regression of our X against our Y
    w0 = reg.coef_  #then extract the cell type proportion coefficents as w 
    w0 = np.array(w0) #and return them as an array to be used in subsequent steps 
    w_print = np.array_str(w0, suppress_small=True)
    print('Initial coefficient estimate:', w_print)# BINOMIAL MODEL
    data_methylated = np.uint64(data_methylated)
    read_depths = np.uint64(read_depths) #need to make sure floating point of r and m are equal to floating point of K and w in binomial model
    w0[w0==0] = 10**-100    
    w_final = np.zeros(k_matrix.shape[0])
    theta = np.log(w0)
    def negative_log_likelihood_theta(theta, r, m, K):
        w = np.exp(theta) #theta originates from the below minimisation step, as these functions are called there
        w /= np.sum(w)
        p = (w[:,np.newaxis] * K).sum(axis=0) #create array p which is the sum, for each cell type, of cell type proportion coeffecient multiplied by the probability of a methylated read at each site across the methylome
        return -np.sum(scipy.special.binom(r, m) + m * np.log(p) + (r - m) * np.log(1 - p)) #return the NLL with previously defined p being used
    res = scipy.optimize.minimize(negative_log_likelihood_theta, theta, 
                            args=(read_depths, data_methylated, k_matrix),
                            method='L-BFGS-B')
    w_final = np.exp(res.x) #Normalise result
    w_final /= np.sum(w_final) #Normalise result # SAVE RESULTS AND RESULTS LOG FILE 
    celltypes = list(k_matrix.ct.values)
    w_df = pd.DataFrame(columns=['ct', 'proportions'])
    w_df.ct = celltypes
    w_df.proportions = w_final
    print(w_df)
    return w_df

##### RUN MULTIPLE DECONVOLUTIONS (E.G. 400, 800, 1600, 3200, 6400, 12800), CALCULATE CETYGO FOR EACH, THEN CHOOSE MOST ACCURATE ESTIMATION ### 
##### POTENTIAL FOR FUTURE HIGH ACCURACY MODE WHERE MOST ACCURATE NUMBER OF SITES FOUND BY GRADIENT DESCENT-APPROACH #####

def multiprocess(args):
    i, cfDNA_data_filenames, min_rd, rmn, target, min_samples, n_sites, output_tags = args
    current_file = cfDNA_data_filenames[i]
    cfDNA_data, high_reads_idx = read_cfdna(current_file, min_rd=min_rd)
    k_matrix = generate_k_matrix(rmn=rmn, target=target, min_samples=min_samples, n_sites=n_sites, high_reads_idx=high_reads_idx)
    w_df = deconvolute(cfDNA_data=cfDNA_data, k_matrix=k_matrix)
    w_df = w_df.rename(columns={'proportions': f"proportions_{output_tags[i]}"})
    return w_df


dataframes = []
if __name__ == '__main__':

    config = pd.read_csv('config/mjf221/loyfer_deconv.csv', index_col=False)
    output_tags = config['output_tag'].tolist()
    cfDNA_data_filenames = config['filenames'].tolist()

    rmn_path = 'reference_data/rmn/rmn_hg38.nc'
    target = 'Neuron'
    min_samples = 1
    min_rd = 200
    n_sites = 2000
    results_filename = 'my_test.csv'
    rmn = xr.open_dataarray(rmn_path)

    args = [(i, cfDNA_data_filenames, min_rd, rmn, target, min_samples, n_sites, output_tags) for i in range(len(cfDNA_data_filenames))]
    with mp.Pool(processes=len(cfDNA_data_filenames)) as pool:
        dataframes = pool.map(multiprocess, args)

    results_df = pd.concat(dataframes, ignore_index=False, axis =1)
    results_df = results_df.loc[:,~results_df.columns.duplicated()]
    results_df.to_csv(results_filename, index = False, float_format='%.4f')



def cetygo_score():
    return