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

os.environ['PYTHONHASHSEED'] = '0'
random.seed(0)
np.random.seed(0)
torch.manual_seed(0)


token_ref_directory = f'{test_name}/'
if os.path.isfile(token_ref_directory + 'adata.h5ad'):
    print(f'Loading {token_ref_directory}adata.h5ad')
    adata = sc.read_h5ad(token_ref_directory + 'adata.h5ad')
    print(f'{token_ref_directory}adata.h5ad loaded')

else:
    # Load genes, barcodes, matrix, and metadata if h5ad is unavailable
    with open(token_ref_directory + 'genes.tsv') as f:
        genes = f.read().rstrip().split('\n')

    with open(token_ref_directory + 'barcodes.tsv') as f:
        barcodes = f.read().rstrip().split('\n')

    print(f'Loading {token_ref_directory}matrix.mtx')
    mat = mmread(token_ref_directory + 'matrix.mtx')
    df = pd.DataFrame.sparse.from_spmatrix(mat, index=genes, columns=barcodes).fillna(0)
    adata = sc.AnnData(df.T)
    print(f'{token_ref_directory}matrix.mtx loaded')
    del mat, df, genes, barcodes

    index_col = 'CellID'
    df_ref_meta = pd.read_csv(token_ref_directory + 'meta.tsv', sep='\t', index_col=index_col)
    df_ref_meta = df_ref_meta.loc[adata.obs_names, :]
    adata.obs = df_ref_meta.copy()

    adata.write(token_ref_directory + 'adata.h5ad')
    print(f'{token_ref_directory}adata.h5ad saved')

adata.obs['group'] = '_' # If adata is for GPT model finetuning, designate one column of adata.obs as "group" that contains group information
adata.obs['isTumor'] = 0 # A trick related to finetuning: only those cells labeled with "isTumor = 0" are to be used for model finetuning
adata.obs = adata.obs[['group', 'isTumor']].copy() # All other columns of adata.obs are excluded as they may disturb tokenization

gp = GProfiler(return_dataframe=True)
if species == 'human':
    df_genes_converted = gp.convert(organism='hsapiens', query=adata.var_names.tolist(), target_namespace='ENSG')
    df_genes_converted = df_genes_converted[~df_genes_converted['incoming'].duplicated()]
    df_genes_converted = df_genes_converted[~df_genes_converted['converted'].isin([None, np.nan, 'None', 'N/A'])]
    df_genes_converted[['incoming', 'converted', 'name', 'description']].to_excel(token_ref_directory + f'{test_name}_convertedGenes.xlsx', index=False)

else:
    df_genes_converted = gp.orth(organism='mmusculus', query=adata.var_names.tolist(), target='hsapiens')
    df_genes_converted = df_genes_converted[~df_genes_converted['incoming'].duplicated()]
    df_genes_converted = df_genes_converted[~df_genes_converted['ortholog_ensg'].isin([None, np.nan, 'None', 'N/A'])]
    df_genes_converted[['incoming', 'converted', 'ortholog_ensg', 'name', 'description']].to_excel(token_ref_directory + f'{test_name}_convertedGenes.xlsx', index=False)

# Filter out those genes with no ENSG IDs
adata = adata[:, df_genes_converted['incoming'].tolist()].copy()

# Add metadata required by tokenizer, don't change the feature names
if species == 'human':
    adata.var['ensembl_id'] = df_genes_converted['converted'].tolist()
else:
    adata.var['ensembl_id'] = df_genes_converted['ortholog_ensg'].tolist()

adata.obs['n_counts'] = np.sum(adata.X.toarray(), axis=1) # total read counts in each cell
adata.obs['filter_pass'] = 1
adata.obs['individual'] = adata.obs.index.tolist() # cell IDs

# Save as [name].loom in [ref_directory]
data = adata.X.toarray().T
df_row_metadata = adata.var.copy()
df_col_metadata = adata.obs.copy()
loompy.create(f'{token_ref_directory}{test_name}.loom', data, df_row_metadata.to_dict('list'), df_col_metadata.to_dict('list'))
del adata, data, df_row_metadata, df_col_metadata

# Tokenize [name].loom
# Ensure [name].loom is the only loom file in [ref_directory]
# Output is a folder [name].dataset in [ref_directory]
tk = TranscriptomeTokenizer({'individual': 'individual', 'isTumor': 'isTumor', 'group': 'group', 'n_counts': 'n_counts'}, nproc=1)
tk.tokenize_data(token_ref_directory, token_ref_directory, test_name)
os.remove(f'{token_ref_directory}{test_name}.loom')
