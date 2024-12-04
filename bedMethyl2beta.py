import numpy as np
import pandas as pd 
import multiprocessing as mp

bedMethyl = pd.read_csv('bulk_seq_data/ONT_cfDNA/P20018_N27_primary_alignments_hg38v2.bed', sep='\t', header=None)
header = ['chrom', 'chromStart', 'chromEnd', 'mod', 'score', 'strand', 'startPos', 'endPos', 'color', 
          'nCov', 'percentMod', 'nMod', 'nCanon', 'nOtherMod', 'nDel', 'nFail', 'nDiff', 'nNoCall']
bedMethyl.columns = header[:len(bedMethyl.columns)]


for i in range(0, len(bedMethyl)):
    if bedMethyl.loc[i, 'strand'] == '-':
        bedMethyl.loc[i, 'chromStart'] = bedMethyl.loc[i, 'chromStart'] - 1 

bedMethyl.loc[:, 'id']  = (bedMethyl.iloc[:, 0].astype(str) + ':' + bedMethyl.iloc[:, 1].astype(str)).tolist()

####cpg idx on to end ####

def retrieve_beta(args):
    i, bedMethyl, rowsPerTask = args 
    chunk_start = round((i-1)*rowsPerTask)
    chunk_end = round(i * rowsPerTask)
    chunk = bedMethyl.iloc[chunk_start:chunk_end, :]
    wgbs_list = []
    for current_idx, group in chunk.groupby('id'):
        nMod_sum = group['nMod'].sum()
        nCov_sum = group['nCov'].sum()
        wgbs_list.append([current_idx, nMod_sum, nCov_sum])
    return pd.DataFrame(wgbs_list, columns=['id', 'nMod_sum', 'nCov_sum'])


if __name__ == '__main__':
    avail_cores = mp.cpu_count()
    rowsPerTask = len(bedMethyl) / avail_cores
    args = [(i, bedMethyl, rowsPerTask) for i in range(1, avail_cores + 1)]
    with mp.Pool(processes=avail_cores) as pool:
        wgbs_chunks = pool.map(retrieve_beta, args)
    wgbs_df = pd.concat(wgbs_chunks, ignore_index=True)
    wgbs_df = wgbs_df.drop_duplicates()
    print(wgbs_df)

wgbs_df.to_csv('P20018_N27_primary_alignments_hg38v2.beta', sep='\t', header=None)