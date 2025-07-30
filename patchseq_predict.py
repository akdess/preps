import os
import random
import numpy as np
import pandas as pd
from joblib import load
from tqdm import tqdm
import argparse


os.environ['PYTHONHASHSEED'] = '0'
random.seed(0)
np.random.seed(0)


parser = argparse.ArgumentParser(description='Electrophysiological feature prediction.')
parser.add_argument('name', help='Input the name of the dataset to be predicted (e.g., mouse, glioma).')
parser.add_argument('-m', '--models', choices=['patchseq', 'allen', 'celltype'], default='patchseq', help='Input -m patchseq or -m allen to designate models (default patchseq).')

args = parser.parse_args()
name = args.name
models = args.models

output_directory = f'{name}_{models}/'
os.mkdir(output_directory)

# ref_embs_directory = 'allen_preds/'
# model_tups = [['v_baseline threshold_v_long_square trough_v_long_square_rel ElasticNet_emb_layer_preds/', 'prediction of v_baseline by embs alpha 0.25 l1_ratio 0.95 MAE 4.523.joblib'], 
#               ['v_baseline threshold_v_long_square trough_v_long_square_rel ElasticNet_emb_layer_scores/', 'prediction of threshold_v_long_square by pcs alpha 0.95 l1_ratio 0.4 MAE 3.806.joblib'], 
#               ['sag ElasticNet_emb_layer_preds/', 'prediction of sag by embs alpha 0.6 l1_ratio 0.0 MAE 0.042.joblib'], 
#               ['latency_rheo upstroke_downstroke_ratio_long_square width_long_square ElasticNet_emb_layer_preds/', 'prediction of latency_rheo_log by embs alpha 0.1 l1_ratio 0.1 MAE 0.06.joblib'], 
#               ['latency_rheo upstroke_downstroke_ratio_long_square width_long_square ElasticNet_emb_layer_preds/', 'prediction of upstroke_downstroke_ratio_long_square_log by embs alpha 0.05 l1_ratio 0.15 MAE 0.512.joblib'], 
#               ['latency_rheo upstroke_downstroke_ratio_long_square width_long_square ElasticNet_emb_layer_scores/', 'prediction of width_long_square by embs alpha 0.05 l1_ratio 0.0 MAE 0.00010707.joblib'], 
#               ['input_resistance rheobase_i peak_v_long_square_rel ElasticNet_emb_layer_scores/', 'prediction of input_resistance_log by embs alpha 0.05 l1_ratio 0.0 MAE 29.919.joblib'], 
#               ['input_resistance rheobase_i peak_v_long_square_rel ElasticNet_emb_layer_scores/', 'prediction of rheobase_i_log by embs alpha 0.05 l1_ratio 0.0 MAE 52.934.joblib'], 
#               ['input_resistance rheobase_i peak_v_long_square_rel ElasticNet_emb_layer_preds/', 'prediction of peak_v_long_square_rel by embs alpha 0.95 l1_ratio 0.0 MAE 5.957.joblib'], 
#               ['adapt_mean tau ElasticNet_emb_layer_preds/', 'prediction of adapt_mean_log by embs alpha 0.35 l1_ratio 0.0 MAE 0.139.joblib'], 
#               ['adapt_mean tau ElasticNet_emb_layer_preds/', 'prediction of tau_log by embs alpha 0.2 l1_ratio 0.0 MAE 0.00745101.joblib']]

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
    model_tups = [['CellTypeLogisticRegression_emb_layer_preds/', 'prediction of Cell Type by embs C 0.04 l1_ratio 0.05 acc 0.7058823529411765.joblib']]

for model_directory, model_name in tqdm(model_tups):
    assert len(model_name.split(' by ')) == 2
    y_name = model_name.split(' by ')[0].replace('prediction of ', '')
    X_name = model_name.split(' by ')[1].split(' ')[0]
    assert X_name in ['embs', 'pcs']
    print(y_name)

    scaler = load(ref_embs_directory + model_directory + 'scaler.joblib')
    pca = load(ref_embs_directory + model_directory + 'pca.joblib')
    model = load(ref_embs_directory + model_directory + model_name)
    if models == 'celltype':
        label_encoder = load(ref_embs_directory + model_directory + 'label_encoder.joblib')

    embs_directory = f'{name}_preds/'
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

    y_predict = model.predict(X_test)

    if y_name.endswith('_log'):
        print('exp back')
        y_predict = np.exp(y_predict) - 1
        y_name = y_name.replace('_log', '')

    if models == 'celltype':
        y_predict = label_encoder.inverse_transform(y_predict)
    
    df = pd.DataFrame({y_name: y_predict}, index=df_merged.index.tolist())
    df.index.name = 'cells'
    df.to_excel(output_directory + f'{y_name}.xlsx')
