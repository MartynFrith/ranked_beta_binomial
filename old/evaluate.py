from xml.etree.ElementTree import tostring
import numpy as np
import pandas as pd
import sys
import re
from sklearn.metrics import explained_variance_score, mean_squared_error, r2_score
np.set_printoptions(suppress=True)


res = []
cols = ['sample', 'r2', 'mse', '1- xvar', 'xmse']

results_filename = sys.argv[1]
results = pd.read_csv(results_filename)
# results_filename  = 'results/10_sample_variable_neuronal_proportion_ZeroInflatedPois_cfDNA.csv'
# results = pd.read_csv(results_filename)

ground_truth_files = pd.read_csv(sys.argv[2])
# ground_truth_files = pd.read_csv('config/neu_percentage_deconv_comparisons/evaluate_config.csv')


col_names = results.columns[1:results.shape[1]]
for column in col_names:
    pattern = r'hg38_([\d.]+%)'
    experiment_target_ct_percentage = re.search(pattern, column).group(1)
    ground_truths = pd.read_csv(ground_truth_files.loc[ground_truth_files.target_ct_proportion == experiment_target_ct_percentage].iloc[0,1])
    truths = ground_truths.proportion
    truths = np.array(truths)/100
    predictions = results[column]
    predictions = np.array(predictions) 
    res.append([column,
            r2_score(truths, predictions),
            mean_squared_error(truths, predictions),
            1-explained_variance_score(truths, predictions),
            mean_squared_error(truths, predictions)*explained_variance_score(truths, predictions)])

res_full = pd.DataFrame(res, columns=cols)

res_full_filename = results_filename[results_filename.find('/')+1:results_filename.find('.')]
res_full.to_csv('{}_deconvolution_metrics.csv'.format(res_full_filename), index=False)


