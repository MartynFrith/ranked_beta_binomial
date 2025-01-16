import pandas as pd 
import numpy as np
import scipy
from sklearn.linear_model import LinearRegression
import xarray as xr
from scipy.stats import entropy
import multiprocessing as mp


def read_cfdna(current_file, min_rd):
    if current_file.endswith('.beta'):
        cfDNA_data = np.fromfile(current_file, dtype=np.uint8).reshape((-1, 2)) 
    elif current_file.endswith('.lbeta'):
        cfDNA_data = np.fromfile(current_file, dtype=np.uint16).reshape((-1, 2))
    elif current_file.endswith('.npy'):           
        cfDNA_data = np.load(current_file)[0].T        
    else:                                              
        cfDNA_data = np.asarray(pd.read_csv(current_file)) 
    cfDNA_data = pd.DataFrame(cfDNA_data)
    cfDNA_data.columns = ['dnam', 'rd']
    cfDNA_data['cpg_idx'] = cfDNA_data.index 
    high_reads_idx = cfDNA_data[np.greater_equal(cfDNA_data.loc[:, 'rd'], min_rd)==True] 
    high_reads_idx = np.array(high_reads_idx['cpg_idx']) 
    return cfDNA_data, high_reads_idx


def entropy_rank(rmn, k_matrix, target, n_sites, high_reads_idx):
    selection = rmn[:,high_reads_idx,:] 
    chance_target_methylated = np.array(selection.loc[target,:,'m'] / selection.loc[target,:,'r'])
    chance_target_not_methylated = 1 - chance_target_methylated
    p_target = np.array([chance_target_methylated, chance_target_not_methylated])
    others = list(selection.coords['ct'].values)
    others.remove(target)
    p_others = selection.loc[others, :, 'm'].sum(axis=0) / selection.loc[others, :, 'r'].sum(axis=0) 
    p_others = np.array([p_others, 1-p_others])
    kl = entropy(p_target, p_others)  
    kl_data_frame = pd.DataFrame({'cpg_idx' : selection.gen.values})
    kl_data_frame['entropy'] =kl 
    entropy_rank = kl_data_frame.sort_values(by='entropy', ascending=False) 
    top_ranked_sites = entropy_rank[0:n_sites] 
    top_ranked_sites_cpg_idx = np.asarray(top_ranked_sites['cpg_idx'])
    sampled_k = k_matrix.sel(gen=top_ranked_sites_cpg_idx)
    return sampled_k


def deconvolute(cfDNA_data, sampled_k):
    cfdna_entropy_filtered = cfDNA_data.loc[sampled_k.gen.values-1, :]
    data_methylated= np.array(cfdna_entropy_filtered.loc[:,'dnam'])
    data_not_methylated = np.array(cfdna_entropy_filtered.loc[:,'rd'] - cfdna_entropy_filtered.loc[:,'dnam'])
    if np.any(np.isnan(sampled_k)==True):
        gen_with_nans = sampled_k.gen.where(sampled_k.isnull().any(dim="ct"), drop=True)
        dropped_gen_values = gen_with_nans.values
        sampled_k = sampled_k.dropna(dim="gen", how="any")

        data_methylated= np.array(cfdna_entropy_filtered[~cfdna_entropy_filtered['cpg_idx'].isin(dropped_gen_values-1)]['dnam'])
        read_depths = np.array(cfdna_entropy_filtered[~cfdna_entropy_filtered['cpg_idx'].isin(dropped_gen_values-1)]['rd'])
    else:
        data_methylated= np.array(cfdna_entropy_filtered.loc[:,'dnam'])
        read_depths = np.sum([data_methylated, data_not_methylated], axis=0) 
    X = sampled_k * read_depths 
    y = data_methylated 
    reg = LinearRegression(positive=True).fit(X.T, y) 
    w0 = reg.coef_  
    w0 = np.array(w0) 
    w_print = np.array_str(w0, suppress_small=True)
    print('Initial coefficient estimate:', w_print)
    data_methylated = np.uint64(data_methylated)
    read_depths = np.uint64(read_depths) 
    w0[w0==0] = 10**-100    
    w_final = np.zeros(sampled_k.shape[0])
    theta = np.log(w0)
    def negative_log_likelihood_theta(theta, r, m, K):
        w = np.exp(theta)
        w /= np.sum(w)
        p = (w[:,np.newaxis] * K).sum(axis=0) 
        return -np.sum(scipy.special.binom(r, m) + m * np.log(p) + (r - m) * np.log(1 - p)) 
    res = scipy.optimize.minimize(negative_log_likelihood_theta, theta, 
                            args=(read_depths, data_methylated, sampled_k),
                            method='L-BFGS-B')
    w_final = np.exp(res.x) 
    w_final /= np.sum(w_final)
    celltypes = list(sampled_k.ct.values)
    w_df = pd.DataFrame(columns=['ct', 'proportions'])
    w_df.ct = celltypes
    w_df.proportions = w_final
    print(w_df)
    return w_df


def multiprocess(args):
    i, cfDNA_data_filenames, min_rd, rmn, k_matrix, target, n_sites, output_tags = args
    current_file = cfDNA_data_filenames[i]
    cfDNA_data, high_reads_idx = read_cfdna(current_file, min_rd=min_rd)
    sampled_k = entropy_rank(rmn=rmn, k_matrix=k_matrix, target=target, n_sites=n_sites, high_reads_idx=high_reads_idx)
    w_df = deconvolute(cfDNA_data=cfDNA_data, sampled_k=sampled_k)
    w_df = w_df.rename(columns={'proportions': f"proportions_{output_tags[i]}"})
    return w_df


dataframes = []
if __name__ == '__main__':
    # USER ARGUMENTS 
    config = 'config/mjf221/testing.csv'
    rmn_path = 'reference_data/rmn/rmn_hg38.nc'
    k_matrix_path = 'reference_data/K_matrices/k_matrix_hg38.nc'
    target = 'Neuron'
    min_rd = 240
    n_sites = 500
    results_filename = 'my_test.csv'
    # LOAD REQUIRED FILES
    load_config = pd.read_csv(config, index_col=False)
    output_tags = load_config['output_tag'].tolist()
    cfDNA_data_filenames = load_config['filenames'].tolist()
    rmn = xr.open_dataarray(rmn_path)
    k_matrix = xr.open_dataarray(k_matrix_path)
    # PERFORM DECONVOLUTION
    args = [(i, cfDNA_data_filenames, min_rd, rmn, k_matrix, target, n_sites, output_tags) for i in range(len(cfDNA_data_filenames))]
    with mp.Pool(processes=len(cfDNA_data_filenames)) as pool:
        dataframes = pool.map(multiprocess, args)
    # WRITE RESULTS
    results_df = pd.concat(dataframes, ignore_index=False, axis =1)
    results_df = results_df.loc[:,~results_df.columns.duplicated()]
    results_df.to_csv(results_filename, index = False, float_format='%.4f')
