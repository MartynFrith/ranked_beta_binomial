import pandas as pd
pd.options.mode.chained_assignment = None
import os
import re
import numpy as np
import multiprocessing as mp


## 1. Load aggregated tissue file 
## 2. Divide the values in the 'r', 'm', and 'n' column by number of folds desired
## 3. For each row and each of the r, m, and n columns select a random number between 0 and the value in the column from
##    a normal distribution (mean = value) to use as the value for that CpG in the present fold. This is the test data.
## 4. Subtract the values in each r, n, and m column  for each CpG in the present fold from the values in the aggregated 
##     tissue file. This is the training data. 
## 5. Repeat k times to get the k folds needed for cross validation? 
##
## want to select from normal distribution but have no way of getting a standard deviation
## for now will do reasonable approximation of standard deviation, 1/4th of range is 1 SD


def chunk_data(training_data):
    chunk_size = np.floor(len(training_data) / mp.cpu_count())
    chunks = []
    for i in range(0, mp.cpu_count() - 1):
        chunk_i_start = int(i * chunk_size)
        chunk_i_end = int((i + 1) * chunk_size)
        chunk_i = training_data.iloc[chunk_i_start:chunk_i_end, :]
        chunks.append(chunk_i)
    final_chunk_start = int(len(chunks) * chunk_size)
    final_chunk_end = len(training_data)
    final_chunk = training_data.iloc[final_chunk_start:final_chunk_end, :]
    chunks.append(final_chunk)
    return chunks


def generate_test_data(chunk):
    chunk_values = chunk.to_numpy()  
    r_vals = chunk['r'].to_numpy()
    n_vals = chunk['n'].to_numpy()
    m_vals = chunk['m'].to_numpy()

    chunk_values[:, 4] = n_vals * 10
    chunk_values[:, 3] = np.round(np.random.normal(loc=r_vals, scale=r_vals / 4))
    chunk_values[:, 5] = np.round(np.random.normal(loc=m_vals, scale=m_vals / 4))
    return pd.DataFrame(chunk_values, columns=chunk.columns)


def process_chunk(chunk):
    return generate_test_data(chunk)



if __name__ == "__main__":
    processed_genome_filenames = os.listdir('reference_data/processed_genome_rnm')
    for file in processed_genome_filenames:
        training_data = pd.read_csv(f'reference_data/processed_genome_rnm/{file}')
        tissue_type  = re.match(r'^(.*?)\.', file).group(1)
        num_folds = 10
        training_data.loc[:, 'r':'m'] = training_data.loc[:, 'r':'m'] / num_folds
        chunks = chunk_data(training_data=training_data)
        for n in range(0, num_folds):
            results = []
            with mp.Pool(mp.cpu_count()) as pool:
                results = pool.map(process_chunk, chunks)
            test_data_set_n = pd.concat(results, ignore_index=True)
            training_data_set_n = training_data.copy()
            training_data_set_n.loc[:,'r':'m'] = training_data.loc[:,'r':'m']*10 - test_data_set_n.loc[:,'r':'m']
            print(training_data_set_n)
            print(test_data_set_n)
            test_data_set_n.to_csv(f'testing/testData_processed_genome/test_data_set{n}_{tissue_type}.csv')
            training_data_set_n.to_csv(f'testing/trainingData_processed_genome/training_data_set{n}_{tissue_type}.csv')
