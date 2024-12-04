#%%
import sys
import xarray as xr
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import datetime
import scipy.optimize
import scipy.special

K_matrix_config_filename = sys.argv[1]
cfDNA_data_filenames = pd.read_csv(sys.argv[2], index_col=False)
cfDNA_data_filenames = cfDNA_data_filenames['filenames'].tolist()
output_folder = sys.argv[3]
deconv_res_file_name = sys.argv[4]
# K_matrix_config_filename = 'K_matrices/K_matrix_config/cedric_0inflatedPoisResultsBased_Kmat_config_hg38_mode.csv'
# cfDNA_data_filenames = 'cfDNA_files/loyfer_cfDNA/GSM6810014_CNVS-NORM-110033633-cfDNA-WGBS-Rep1.beta'
# cfDNA_data_filename = 'cfDNA_files/zero_inflated_poisson_cfDNA_full_hg38_30rds.npy'

#%% READ CFDNA DATA
dataframes = []
for file in cfDNA_data_filenames:
    if file.endswith('.beta'):
        cfDNA_data = np.fromfile(file, dtype=np.uint8).reshape((-1, 2))#load cfDNA methylation sequencing data as numpy array 
    elif file.endswith('.lbeta'):
        cfDNA_data = np.fromfile(file, dtype=np.uint16).reshape((-1, 2))
    elif file.endswith('.npy'):             #OR USE SYNTHETIC CFDNA DATA
        cfDNA_data = np.load(file).T        #select only the first 2d numpy array in synthesised cfDNA
    else:                                                  #as is a 3d array, 3rd dimension is individuals
        cfDNA_data = np.asarray(pd.read_csv(file)) 
    cfDNA_data = pd.DataFrame(cfDNA_data)
    cfDNA_data.columns = ['dnam', 'rd']
    dnam = pd.DataFrame(cfDNA_data['dnam'])
    dnam.columns = ['dnam']
    dnanm = pd.DataFrame(cfDNA_data['rd'] - dnam['dnam'])
    cfdna = np.array([dnam.values.flatten(), dnanm.values.flatten()])
    K_matrix_config = pd.read_csv(K_matrix_config_filename)
    mu_rd = round(cfDNA_data['rd'].mean()/10)*10
    if mu_rd <= 100 and mu_rd >= 10: #currently highest mean read depth K matrix created for is 100
        K_matrix_filename = K_matrix_config.loc[K_matrix_config['mu_rd'] == mu_rd, 'filename'].values[0]
    elif mu_rd < 10:
        K_matrix_filename = K_matrix_config.loc[K_matrix_config['mu_rd'] == 10, 'filename'].values[0]
    else:
        K_matrix_filename = K_matrix_config.loc[K_matrix_config['mu_rd'] == 100, 'filename'].values[0]
    sampled_K = xr.open_dataarray('K_matrices/{}'.format(K_matrix_filename))
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
    read_depths = np.uint64(read_depths) #need to make sure floating point of r and m are equal to floating point of K and w in binomial model
    w0[w0==0] = 10**-100    
    w_final = np.zeros(sampled_K.shape[0])
    theta = np.log(w0)
    def negative_log_likelihood_theta(theta, r, m, K):
        w = np.exp(theta) #theta originates from the below minimisation step, as these functions are called there
        w /= np.sum(w)
        p = (w[:,np.newaxis] * K).sum(axis=0) #create array p which is the sum, for each cell type, of cell type proportion coeffecient multiplied by 
                                        #the probability of a methylated read at each site across the methylome
        return -np.sum(scipy.special.binom(r, m) + m * np.log(p) + (r - m) * np.log(1 - p)) #return the NLL with previously defined p being used
    res = scipy.optimize.minimize(negative_log_likelihood_theta, theta, 
                            args=(read_depths, data_methylated, sampled_K),
                            method='L-BFGS-B')
    w_final = np.exp(res.x) #Normalise result
    w_final /= np.sum(w_final) #Normalise result


    #%% SAVE RESULTS AND RESULTS LOG FILE 
    celltypes = list(sampled_K.ct.values)
    w_df = pd.DataFrame(columns=['ct', 'proportions'])
    w_df.ct = celltypes
    w_df.proportions = w_final
    filename = file[file.find('/') + 1 : file.find('.') +1]
    filename = filename[filename.find('/')+1: filename.find('.')]
    prop_col = f'proportion_{filename}'
    w_df.rename(columns={'proportions':prop_col}, inplace = True)
    print('Result:', res)
    print('Cell type proportion estimate:', w_final)
    current_datetime = datetime.datetime.now()
    formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H%-M-%S")
    K_matrix_filename = K_matrix_filename [K_matrix_filename .find('h') : K_matrix_filename .find('.')]
    log_filename = f"{output_folder}logs/{filename}_{K_matrix_filename }_{formatted_datetime}.txt"
    with open(log_filename, 'w') as log_file:
        log_file.write('NNLS initial coefficient estimation:\n{}\n'.format(w_print))
        log_file.write('Minimisation ouput:\n{}\n'.format(res))
        log_file.write('Final deconvoluted cell type proportions:\n{}\n'.format(w_final))
    dataframes.append(w_df)


#%% WRITE OUT RESULTS
results_df = pd.concat(dataframes, ignore_index=False, axis =1)
results_df = results_df.loc[:,~results_df.columns.duplicated()]
K_matrix_config_filename = K_matrix_config_filename [K_matrix_config_filename.find('g/')+2: K_matrix_config_filename.find('.')]
results_filename = f"{output_folder}{deconv_res_file_name}_{K_matrix_config_filename}_{formatted_datetime}.csv"
results_df.to_csv(results_filename, index = False, float_format='%.4f')

#%%