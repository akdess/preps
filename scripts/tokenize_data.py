import argparse
parser = argparse.ArgumentParser(description='scRNA-seq data tokenization.')
parser.add_argument('test_name', help='Input the directory name of the dataset to be tokenized (e.g., mouse, glioma).')
parser.add_argument('-s', '--species', choices=['human', 'mouse'], default='human', help='Input -s human or -s mouse to designate species (default human).')
args = parser.parse_args()
test_name = args.test_name
species = args.species


import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
import torch
import loompy
from gprofiler import GProfiler
from geneformer import TranscriptomeTokenizer
from scipy.sparse import issparse


os.environ['PYTHONHASHSEED'] = '0'
random.seed(0)
np.random.seed(0)
torch.manual_seed(0)


output_directory = f'{test_name}/'
os.makedirs(output_directory, exist_ok=True)

input_directory = './'
file_path = f'{input_directory}{test_name}.h5ad'
adata = sc.read_h5ad(file_path)

# else:
#     # Load genes, barcodes, matrix, and metadata if h5ad is unavailable
#     with open(input_directory + 'genes.tsv') as f:
#         genes = f.read().rstrip().split('\n')

#     with open(input_directory + 'barcodes.tsv') as f:
#         barcodes = f.read().rstrip().split('\n')

#     adata = sc.read_mtx(input_directory + 'matrix.mtx').T 
#     adata.obs_names = barcodes
#     adata.var_names = genes

#     adata.write(output_directory + 'adata.h5ad')
#     print(f'{output_directory}adata.h5ad saved')


# Preserve 'group' if it exists (critical for finetuning), otherwise set a placeholder
if 'group' not in adata.obs.columns:
    adata.obs['group'] = '_' 

# Preserve 'isTumor' if it exists, otherwise set to 0 (so cells are used by default)
if 'isTumor' not in adata.obs.columns:
    adata.obs['isTumor'] = 0 

# All other columns of adata.obs are excluded as they may disturb tokenization
adata.obs = adata.obs[['group', 'isTumor']].copy()

gp = GProfiler(return_dataframe=True)
if species == 'human':
    df_converted = gp.convert(organism='hsapiens', query=adata.var_names.tolist(), target_namespace='ENSG')
    converted_col = 'converted'

else:
    df_converted = gp.orth(organism='mmusculus', query=adata.var_names.tolist(), target='hsapiens')
    converted_col = 'ortholog_ensg'

df_converted = df_converted[~df_converted[converted_col].isin([None, np.nan, 'None', 'N/A'])]
df_converted = df_converted.drop_duplicates(subset=['incoming'])
df_converted = df_converted.drop_duplicates(subset=[converted_col])
df_converted[['incoming', converted_col, 'name', 'description']].to_excel(output_directory + f'{test_name}_convertedGenes.xlsx', index=False)

# Filter out those genes with no ENSG IDs
adata = adata[:, df_converted['incoming'].tolist()].copy()

# Add metadata required by tokenizer, don't change the feature names
adata.var['ensembl_id'] = df_converted[converted_col].tolist()
adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))
adata.obs['filter_pass'] = 1
adata.obs['individual'] = adata.obs.index.tolist()

# Save as [name].loom in [ref_directory]
loom_path = f'{output_directory}{test_name}.loom'
adata.write_loom(loom_path, write_obsm_varm=False)
del adata, df_converted

# Tokenize [name].loom
# Ensure [name].loom is the only loom file in [ref_directory]
# Output is a folder [name].dataset in [ref_directory]
tk = TranscriptomeTokenizer({'individual': 'individual', 'isTumor': 'isTumor', 'group': 'group', 'n_counts': 'n_counts'}, nproc=1)
tk.tokenize_data(output_directory, output_directory, test_name)
os.remove(loom_path)
