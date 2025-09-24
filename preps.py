import argparse
parser = argparse.ArgumentParser(description='scRNA-seq data tokenization, annotation, and prediction.')
parser.add_argument('test_name', help='Input the directory name of the dataset to be tokenized (e.g., mouse, glioma).')
parser.add_argument('-s', '--species', choices=['human', 'mouse'], default='human', help='Input -s human or -s mouse to designate species (default human).')
parser.add_argument('-g', '--gpu_name', choices=list(map(str, range(1000))), default='0', help='Input the idle GPU on which to run the code (e.g., 0, 1, 2).')
parser.add_argument('-m', '--models', choices=['patchseq', 'allen', 'celltype'], default='patchseq', help='Input -m patchseq or -m celltype to designate models (default patchseq).')
args = parser.parse_args()
test_name = args.test_name
species = args.species
gpu_name = args.gpu_name
models = args.models

import os
os.environ['CUDA_VISIBLE_DEVICES'] = gpu_name
os.environ['NCCL_DEBUG'] = 'INFO'

import random
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
import torch
import loompy
from gprofiler import GProfiler
from geneformer import TranscriptomeTokenizer

from datasets import load_from_disk
from sklearn.metrics import accuracy_score, f1_score
from transformers import BertForSequenceClassification
from transformers import Trainer
from transformers.training_args import TrainingArguments
from geneformer import DataCollatorForCellClassification
from scipy.special import softmax
import shutil
from joblib import load
from tqdm import tqdm

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

ann_output_directory = f'{test_name}_preds/'
os.mkdir(ann_output_directory)
shutil.copytree(f'{test_name}/{test_name}.dataset', ann_output_directory + 'tokenized_copy.dataset')

ref_num_tups = [('aldinger_2000perCellType', 21), 
                ('allen_2000perCellType', 20), 
                ('bhaduri_3000perCellType', 10), 
                ('bhaduri_d2_4000perCellType', 10), 
                ('codex_1000perCellType', 16), 
                ('devbrain_3000perCellType', 10), 
                ('dirks_primary_gbm_combined_2000perCellType', 13), 
                ('primary_gbm_2000perCellType', 8), 
                ('recurrent_gbm_1000perCellType', 14), 
                ('TissueImmune_2000perCellType', 45)]

for ref_name, num_classes in tqdm(ref_num_tups):
    output_prefix = f'preds_by_{ref_name}_num_classes_{num_classes}'
    model_directory=f'{ref_name}/finetune/240605_geneformer_CellClassifier_0_L2048_B12_LR5e-05_LSlinear_WU500_E10_Oadamw_F0/'


    # load data, labels, model
    tokenized_dataset = load_from_disk(ann_output_directory + 'tokenized_copy.dataset')
    # tokenized_dataset = tokenized_dataset.filter(lambda x: x['isTumor'] == 1, num_proc=1)
    labels = [0] * tokenized_dataset.num_rows
    tokenized_dataset = tokenized_dataset.add_column('label', labels)

    df_target_names = pd.read_excel(f'{ref_name}/finetune/target_names.xlsx', header=None)


    def compute_metrics(pred):
        labels = pred.label_ids
        preds = pred.predictions.argmax(-1)
        # calculate accuracy and macro f1 using sklearn's function
        acc = accuracy_score(labels, preds)
        macro_f1 = f1_score(labels, preds, average='macro')
        return {
        'accuracy': acc,
        'macro_f1': macro_f1
        }


    # set model parameters
    # max input size
    max_input_size = 2 ** 11  # 2048
    # number gpus
    num_gpus = 1
    # batch size for training and eval
    geneformer_batch_size = 12

    model = BertForSequenceClassification.from_pretrained(model_directory, 
                                                        num_labels=num_classes,
                                                        output_attentions=False,
                                                        output_hidden_states=False).to('cuda')

    # predict
    training_args = {
        "do_train": False,
        "do_eval": False,
        "evaluation_strategy": "epoch",
        "group_by_length": True,
        "length_column_name": "length",
        "disable_tqdm": False,
        "per_device_train_batch_size": geneformer_batch_size,
        "per_device_eval_batch_size": geneformer_batch_size,
        "output_dir": ann_output_directory
    }
    training_args_init = TrainingArguments(**training_args)
    trainer = Trainer(
        model=model,
        args=training_args_init,
        data_collator=DataCollatorForCellClassification(),
        train_dataset=tokenized_dataset,
        eval_dataset=tokenized_dataset
    )
    predictions = trainer.predict(tokenized_dataset)

    df_preds = pd.DataFrame(predictions.predictions, index=tokenized_dataset['individual'])
    df_preds.index.name = 'individual'
    df_preds.to_csv(ann_output_directory + f'{output_prefix}_preds.csv', sep=',')

    y_predict_prob = softmax(predictions.predictions, axis=1)
    y_predict_int = np.argmax(y_predict_prob, axis=1)
    y_predict_class = [df_target_names.iat[i, 0] for i in y_predict_int]
    y_predict_score = np.amax(y_predict_prob, axis=1)

    # save results
    col_ann = f'{ref_name}_ann'
    col_score = f'{ref_name}_score'

    df_predict_prob = pd.DataFrame(y_predict_prob, index=tokenized_dataset['individual'])
    df_predict_prob.index.name = 'individual'
    df_predict_prob.columns = df_target_names.iloc[:, 0].tolist()
    df_predict_prob[col_ann] = y_predict_class
    df_predict_prob[col_score] = y_predict_score
    df_predict_prob.to_csv(ann_output_directory + f'{output_prefix}_scores.csv', sep=',')

pred_output_directory = f'{test_name}_{models}/'
os.mkdir(pred_output_directory)

if models == 'allen':
    ref_embs_directory = 'allen_preds/'
    model_tups = [['all ephys ElasticNet_emb_layer_preds/', 'prediction of AP width (ms) (Allen model)_log by pcs alpha 0.15 l1_ratio 0.05 positive False MAE 0.098.joblib'], 
                ['all ephys ElasticNet_emb_layer_preds/', 'prediction of Fitted MP (mV) (Allen model) by embs alpha 0.45 l1_ratio 0.3 positive False MAE 4.608.joblib'], 
                ['all ephys ElasticNet_emb_layer_preds/', 'prediction of Upstroke-to-downstroke ratio (Allen model) by embs alpha 0.05 l1_ratio 0.7 positive False MAE 0.51.joblib'], 
                ['all ephys ElasticNet_emb_layer_preds/', 'prediction of AP threshold (mV) (Allen model) by pcs alpha 0.95 l1_ratio 0.95 positive True MAE 3.808.joblib'], 
                ['all ephys ElasticNet_emb_layer_scores/', 'prediction of Membrane time constant (ms) (Allen model)_log by embs alpha 0.05 l1_ratio 0.0 positive False MAE 5.995.joblib'], 
                ['all ephys ElasticNet_emb_layer_preds/', 'prediction of Sag ratio (Allen model) by pcs alpha 0.25 l1_ratio 0.0 positive False MAE 0.041.joblib'], 
                ['all ephys ElasticNet_emb_layer_scores/', 'prediction of Rheobase (pA) (Allen model)_log by embs alpha 0.05 l1_ratio 0.0 positive False MAE 52.222.joblib'], 
                ['all ephys ElasticNet_emb_layer_preds/', 'prediction of AP amplitude (mV) (Allen model) by embs alpha 0.3 l1_ratio 0.95 positive False MAE 5.892.joblib'], 
                ['all ephys ElasticNet_emb_layer_preds/', 'prediction of Latency (ms) (Allen model)_log by embs alpha 0.8 l1_ratio 0.05 positive False MAE 59.476.joblib'], 
                ['all ephys ElasticNet_emb_layer_scores/', 'prediction of Input resistance (MOhm) (Allen model)_log by embs alpha 0.05 l1_ratio 0.0 positive False MAE 29.347.joblib']]
elif models == 'patchseq':
    ref_embs_directory = 'combined_patchseq_all_preds/'
    model_tups = [['Fitted MP (mV) AP threshold (mV) Afterhyperpolarization (mV) ElasticNet_emb_layer_preds/', 'prediction of Fitted MP (mV) by embs alpha 0.95 l1_ratio 0.0 MAE 17.014.joblib'], 
                ['Fitted MP (mV) AP threshold (mV) Afterhyperpolarization (mV) ElasticNet_emb_layer_preds/', 'prediction of AP threshold (mV) by pcs alpha 0.95 l1_ratio 0.7 MAE 8.599.joblib'], 
                ['Rheobase (pA) Sag ratio Membrane time constant (ms) ElasticNet_emb_layer_scores/', 'prediction of Sag ratio_log by embs alpha 0.2 l1_ratio 0.0 MAE 0.086.joblib'], 
                ['AP width (ms) Upstroke-to-downstroke ratio Latency (ms) ElasticNet_emb_layer_preds/', 'prediction of Latency (ms)_log by embs alpha 0.05 l1_ratio 0.85 MAE 58.346.joblib'], 
                ['AP width (ms) Upstroke-to-downstroke ratio Latency (ms) ElasticNet_emb_layer_preds/', 'prediction of Upstroke-to-downstroke ratio by embs alpha 0.85 l1_ratio 0.05 MAE 1.596.joblib'], 
                ['AP width (ms) Upstroke-to-downstroke ratio Latency (ms) ElasticNet_emb_layer_preds/', 'prediction of AP width (ms)_log by embs alpha 0.3 l1_ratio 0.05 MAE 0.915.joblib'], 
                ['Input resistance (MOhm) AP amplitude (mV) Max number of APs ElasticNet_emb_layer_preds/', 'prediction of Input resistance (MOhm)_log by embs alpha 0.95 l1_ratio 0.0 MAE 470.228.joblib'], 
                ['Rheobase (pA) Sag ratio Membrane time constant (ms) ElasticNet_emb_layer_scores/', 'prediction of Rheobase (pA)_log by embs alpha 0.05 l1_ratio 0.0 MAE 34.342.joblib'], 
                ['Input resistance (MOhm) AP amplitude (mV) Max number of APs ElasticNet_emb_layer_preds/', 'prediction of AP amplitude (mV) by embs alpha 0.95 l1_ratio 0.7 MAE 17.56.joblib'], 
                ['ISI adaptation index ElasticNet_emb_layer_preds/', 'prediction of ISI adaptation index by embs alpha 0.15 l1_ratio 0.75 MAE 0.457.joblib'], 
                ['Rheobase (pA) Sag ratio Membrane time constant (ms) ElasticNet_emb_layer_scores/', 'prediction of Membrane time constant (ms)_log by embs alpha 0.1 l1_ratio 0.05 MAE 16.943.joblib']]

else:
    ref_embs_directory = 'combined_patchseq_all_preds/'
    model_tups = [['CellTypeLogisticRegression_emb_layer_scores/', 'prediction of Cell Type by embs C 3 l1_ratio 0.95 acc 0.7058823529411765.joblib']]

for model_directory, model_name in tqdm(model_tups):
    assert len(model_name.split(' by ')) == 2
    y_name = model_name.split(' by ')[0].replace('prediction of ', '')
    X_name = model_name.split(' by ')[1].split(' ')[0]
    assert X_name in ['embs', 'pcs']
    print(y_name)

    scaler = load(ref_embs_directory + model_directory + 'scaler.joblib')
    pca = load(ref_embs_directory + model_directory + 'pca.joblib')
    preps_model = load(ref_embs_directory + model_directory + model_name)
    if models == 'celltype':
        label_encoder = load(ref_embs_directory + model_directory + 'label_encoder.joblib')

    embs_directory = f'{test_name}_preds/'
    emb_layer = model_directory.split('_')[-1].replace('/', '')
    assert emb_layer in ['preds', 'scores']
    embs_files = [x for x in os.listdir(embs_directory) if x.endswith(f'{emb_layer}.csv')]
    embs_files = sorted(embs_files)

    # ensure glioma embeddings are concatenated in the same order as patchseq training data
    ref_embs_files = [x for x in os.listdir(ref_embs_directory) if x.endswith(f'{emb_layer}.csv')]
    ref_embs_files = sorted(ref_embs_files)
    assert all([x == y for x, y in zip(ref_embs_files, embs_files)])

    df_merged = None
    for file_name in embs_files:
        df = pd.read_csv(embs_directory + file_name, sep=',', index_col='individual')
        
        if emb_layer in [-1, 0]:
            df = df.drop(columns=['Unnamed: 0', 'group'])
            assert df.shape[1] == 256
        elif emb_layer == 'features':
            df = df.iloc[:, :-2]
            assert df.shape[1] == 32
        else:
            if emb_layer == 'scores':
                df = df.iloc[:, :-2]
            assert df.shape[1] == int(file_name.split('_')[-2])
        
        df.columns = [f"{file_name.replace('.csv', '')}_dim_{x}" for x in df.columns.tolist()]
        
        if df_merged is None:
            df_merged = df.copy()
        else:
            df_merged = pd.merge(df_merged, df, how='inner', left_index=True, right_index=True)

    if X_name == 'embs':
        X_test = df_merged.values
    else:
        X_test = pca.transform(scaler.transform(df_merged.values))

    y_predict = preps_model.predict(X_test)

    if y_name.endswith('_log'):
        print('exp back')
        y_predict = np.exp(y_predict) - 1
        y_name = y_name.replace('_log', '')

    if models == 'celltype':
        y_predict = label_encoder.inverse_transform(y_predict)
        y_score = np.amax(preps_model.predict_proba(X_test), axis=1)
    
    df = pd.DataFrame({y_name: y_predict}, index=df_merged.index.tolist())
    if models == 'celltype':
        df['prob'] = y_score
    df.index.name = 'cells'
    df.to_excel(pred_output_directory + f'{y_name}.xlsx')
