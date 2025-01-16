import numpy as np
import pandas as pd 


cpg_idx = pd.read_csv('reference_data/CpG.bed', sep='\t', header=None)
cpg_idx.loc[:, 'id']  = (cpg_idx.iloc[:, 0].astype(str) + ':' + cpg_idx.iloc[:, 1].astype(str)).tolist()
cpg_idx = cpg_idx.iloc[:,2:4]

bedMethyl = pd.read_csv('P20018_N27_primary_alignments_hg38v2_trad.bed', sep='\t', header=None, usecols = [0,1,2,5,9,11])
header = ['chrom', 'chromStart', 'chromEnd','strand', 'nCov', 'nMod']
bedMethyl.columns = header[:len(bedMethyl.columns)]
reverse_mask = bedMethyl['strand'] == '-'  
bedMethyl.loc[reverse_mask, 'chromEnd'] -= 1 
bedMethyl.loc[:, 'id']  = (bedMethyl.iloc[:, 0].astype(str) + ':' + bedMethyl.iloc[:, 2].astype(str)).tolist()


bedMethyl_agg = bedMethyl.groupby('id', as_index=False).agg({'nMod': 'sum', 'nCov': 'sum'})
bedMethyl_agg = pd.merge(bedMethyl_agg, cpg_idx, how='left', on='id')

bedMethyl_agg.to_csv('P20018_N27_primary_alignments_hg38v2.beta', sep='\t', header=None)




bedMethyl_merge = pd.merge(bedMethyl, cpg_idx, how='left', on='id')
bedMethyl_idxs_only = bedMethyl_merge[bedMethyl_merge[2].notna()]
bedMethyl_agg = bedMethyl_idxs_only.groupby('id', as_index=False).agg({'nMod': 'sum', 'nCov': 'sum'})


bedMethyl[bedMethyl['id'] =='chr1:10577']
cpg_idx['id'].isin(bedMethyl['id']).all()
not_in_bedMethyl = cpg_idx.loc[~cpg_idx['id'].isin(bedMethyl['id']), 'id'].tolist()
