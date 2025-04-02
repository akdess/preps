import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import loompy
from geneformer import TranscriptomeTokenizer


os.environ['PYTHONHASHSEED'] = '0'
random.seed(0)
np.random.seed(0)
torch.manual_seed(0)


ref_directory = './'
name = 'WHB-10Xv3-Neurons-raw'
adata = sc.read_h5ad(f'{ref_directory}{name}.h5ad')
print('original', adata.shape)

sc.pp.subsample(adata, fraction=0.1)
print('subsampled', adata.shape)

adata.obs['group'] = adata.obs['anatomical_division_label'].values # Designate one column as "group" that contains groups information for model finetuning
adata.obs['isTumor'] = 0 # A trick related to finetuning: only those cells labeled with "isTumor = 0" are to be used for model finetuning
adata.obs = adata.obs[['group', 'isTumor']].copy()

adata.var['ensembl_id'] = adata.var.index.tolist()
adata.obs['n_counts'] = np.sum(adata.X.toarray(), axis=1)
adata.obs['filter_pass'] = 1
adata.obs['individual'] = adata.obs.index.tolist()

vmax = 2 ** 30
adata = adata[adata.obs['n_counts'].map(lambda x: x < vmax)].copy()
print('filtered', adata.shape)
print('max of n_counts', adata.obs['n_counts'].max())

data = adata.X.toarray().T
df_row_metadata = adata.var.copy()
df_col_metadata = adata.obs.copy()
loompy.create(f'{ref_directory}{name}.loom', data, df_row_metadata.to_dict('list'), df_col_metadata.to_dict('list'))

del adata, data

tk = TranscriptomeTokenizer({'individual': 'individual', 'isTumor': 'isTumor', 'group': 'group', 'n_counts': 'n_counts'}, nproc=1)
tk.tokenize_data(ref_directory, ref_directory, name)
os.remove(f'{ref_directory}{name}.loom')
