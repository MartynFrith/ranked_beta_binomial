import pandas as pd 
import numpy as np
import scipy
from sklearn.linear_model import LinearRegression
import xarray as xr
import multiprocessing as mp


def read_cfdna(current_file):
    if current_file.endswith('.beta'):
        cfDNA_data = np.fromfile(current_file, dtype=np.uint8).reshape((-1, 2))
    elif current_file.endswith('.lbeta'):
        cfDNA_data = np.fromfile(current_file, dtype=np.uint16).reshape((-1, 2))
    elif current_file.endswith('.npy'):
        cfDNA_data = np.load(current_file).T
    else:
        cfDNA_data = np.asarray(pd.read_csv(current_file))
    cfDNA_data = pd.DataFrame(cfDNA_data, columns=['dnam', 'rd'])
    dnam = pd.DataFrame(cfDNA_data['dnam'])
    read_depths = pd.DataFrame(cfDNA_data['rd'])
    cfDNA_data = np.array([dnam.values.flatten(), read_depths.values.flatten()])
    mu_rd = read_depths.mean().values
    if mu_rd > 100: 
        mu_rd = 100
    elif mu_rd < 10:
        mu_rd = 10
    else:
        mu_rd = int(round(read_depths.mean()/10)*10)
    return cfDNA_data, mu_rd

def negative_log_likelihood_theta(theta, r, m, K):
    w = np.exp(theta)
    w /= np.sum(w)
    p = (w[:,np.newaxis] * K).sum(axis=0) 
    return -np.sum(scipy.special.binom(r, m) + m * np.log(p) + (r - m) * np.log(1 - p)) 

def compute_normalised_rmse(dnam, read_depths, sampled_k, w_df):
    measured_meth = np.divide(dnam, read_depths)
    measured_meth[np.where(np.isnan(measured_meth))] = 0
    predicted_meth = [w_df[w_df['ct'] == i]['proportions'].values[0] * sampled_k.loc[i,:] for i in sampled_k.ct.values]
    predicted_meth =  np.sum([i.values for i in predicted_meth], axis=0)
    predicted_meth[np.where(np.isnan(predicted_meth))] = 0
    rmse = np.sqrt((sum(np.subtract(measured_meth, predicted_meth)**2))/len(sampled_k.gen))
    normalised_rmse = rmse*np.sqrt(np.nanmean(read_depths))
    return normalised_rmse

def deconvolute(cfDNA_data, k_matrix, n_sites):
    sampled_k = k_matrix.isel(gen=slice(0, n_sites)) 
    sites_mask = sampled_k.gen.values - 1
    sampled_cfDNA = cfDNA_data[:, sites_mask]
    dnam, read_depths = sampled_cfDNA
    X = sampled_k * read_depths 
    y = dnam 
    reg = LinearRegression(positive=True).fit(X.T, y)
    w0 = reg.coef_  
    w0[np.where(w0 < 10e-4)] = 0
    w0[w0==0] = 10**-100    
    w_final = np.zeros(sampled_k.shape[0])
    theta = np.log(w0)
    dnam = np.uint64(dnam)
    read_depths = np.uint64(read_depths) 
    res = scipy.optimize.minimize(negative_log_likelihood_theta, theta, 
                            args=(read_depths, dnam, sampled_k),
                            method='L-BFGS-B')
    w_final = np.exp(res.x) 
    w_final /= np.sum(w_final)
    celltypes = list(sampled_k.ct.values)
    w_df = pd.DataFrame(columns=['ct', 'proportions'])
    w_df.ct = celltypes
    w_df.proportions = w_final
    w_df.loc[w_df['proportions'] < 10e-5, 'proportions'] = 0   
    normalised_rmse = compute_normalised_rmse(dnam, read_depths, sampled_k, w_df)
    return w_df, normalised_rmse 

def multiprocess(args):
    i, cfDNA_data_filenames, site_rd_config, k_matrix, output_tags = args
    current_file = cfDNA_data_filenames[i]
    cfDNA_data, mu_rd = read_cfdna(current_file) 
    n_sites = site_rd_config.loc[site_rd_config['mu_rd'] == mu_rd, 'n_sites'].values[0]
    w_df, score = deconvolute(cfDNA_data, k_matrix, n_sites)
    score = pd.DataFrame({'sample':output_tags[i], 'score':score}, index=[0])
    w_df = w_df.rename(columns={'proportions': f"proportions_{output_tags[i]}"})
    return w_df, score


dataframes = []
if __name__ == '__main__':
# USER ARGUMENTS 
    config_path = 'config/mjf221/loyfer_deconv.csv'
    k_matrix_path = 'reference_data/K_matrices/kMatrix_NeuronJSRanked_full_hg38.nc'
    results_filename = 'my_test.csv'
    # LOAD REQUIRED FILES
    config = pd.read_csv(config_path, index_col=False)
    output_tags = config['output_tag'].tolist()
    cfDNA_data_filenames = config['filenames'].tolist()
    k_matrix = xr.open_dataarray(k_matrix_path)
    mu_rd_arr = [10, 20, 30, 40 ,50, 60, 70, 80, 90, 100]
    n_sites_arr = [278255, 215443, 129154, 166810, 100000, 278255, 129154, 166810, 129154, 129154]
    site_rd_config = pd.DataFrame({'mu_rd':mu_rd_arr, 'n_sites':n_sites_arr})
    # PERFORM DECONVOLUTION
    args = [(i, cfDNA_data_filenames, site_rd_config, k_matrix, output_tags) for i in range(len(cfDNA_data_filenames))]
    with mp.Pool(processes=len(cfDNA_data_filenames)) as pool:
        res_list = pool.map(multiprocess, args)
    dataframes, scores = zip(*res_list)
    # WRITE RESULTS
    results_df = pd.concat(dataframes, ignore_index=False, axis =1)
    results_df = results_df.loc[:,~results_df.columns.duplicated()]
    score_df = pd.concat(scores, ignore_index=False, axis =0)
    results_df.to_csv(results_filename, index = False, float_format='%.4f')
