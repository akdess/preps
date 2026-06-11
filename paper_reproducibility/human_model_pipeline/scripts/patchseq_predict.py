import argparse
parser = argparse.ArgumentParser(description='Electrophysiological feature prediction.')
parser.add_argument('test_name', help='Input the name of the dataset to be predicted (e.g., mouse, glioma).')
parser.add_argument('-m', '--models', choices=['patchseq', 'm1_patchseq_total', 'allen', 'celltype', 'm1_celltype_total', 'm1_celltype_exons', 'rf', 'rfhvg', 'mlp', 'mlphvg', 'eln', 'elnhvg'], default='m1_celltype_exons', help='Input -m patchseq or -m celltype to designate models (default patchseq).')
args = parser.parse_args()
test_name = args.test_name
models = args.models

# test_name = 'glioma2'
# models = 'rf'

import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
from joblib import load
from tqdm import tqdm
from scipy.sparse import issparse

os.environ['PYTHONHASHSEED'] = '0'
random.seed(0)
np.random.seed(0)


pred_output_directory = f'{test_name}_{models}/'
os.mkdir(pred_output_directory)

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
elif models == 'm1_patchseq_total':
    ref_embs_directory = 'm1_patchseq_total_preds/'
    model_tups = [['ElasticNet_emb_layer_preds/', 'prediction of AP width (ms) by embs n 1056 alpha 0.01 l1 0.3 corr 0.766 p 2.7039197303039005e-52 MAE 0.207.joblib'], 
                  ['ElasticNet_emb_layer_preds/', 'prediction of ISI adaptation index_log by embs n 1056 alpha 0.1 l1 0.3 corr 0.378 p 2.2723436696849558e-10 MAE 0.643.joblib'], 
                  ['ElasticNet_emb_layer_preds/', 'prediction of Upstroke-to-downstroke ratio by embs n 1056 alpha 0.01 l1 1.0 corr 0.781 p 1.903381747446319e-55 MAE 0.598.joblib'], 
                  ['ElasticNet_emb_layer_preds/', 'prediction of AP threshold (mV) by pcs n 1056 alpha 1.0 l1 0.01 corr 0.157 p 0.010636072782549555 MAE 5.455.joblib'], 
                  ['ElasticNet_emb_layer_preds/', 'prediction of Membrane time constant (ms)_log by embs n 1056 alpha 0.1 l1 0.01 corr 0.33 p 4.1451304498668005e-08 MAE 4.412.joblib'], 
                  ['ElasticNet_emb_layer_preds/', 'prediction of Sag ratio_log by embs n 1056 alpha 0.001 l1 0.9 corr 0.284 p 2.7605680508399486e-06 MAE 0.072.joblib'], 
                  ['ElasticNet_emb_layer_preds/', 'prediction of Rheobase (pA)_log by embs n 1056 alpha 0.1 l1 0.1 corr 0.484 p 6.418915574492753e-17 MAE 42.214.joblib'], 
                  ['ElasticNet_emb_layer_preds/', 'prediction of AP amplitude (mV) by embs n 1056 alpha 0.1 l1 0.7 corr 0.58 p 3.9612841243281993e-25 MAE 6.519.joblib'], 
                  ['ElasticNet_emb_layer_preds/', 'prediction of Latency (ms)_log by embs n 1056 alpha 0.1 l1 0.1 corr 0.081 p 0.18727974493101324 MAE 64.512.joblib'], 
                  ['ElasticNet_emb_layer_preds/', 'prediction of Input resistance (MOhm) by embs n 1056 alpha 1.0 l1 0.9 corr 0.58 p 3.5174679368441846e-25 MAE 70.739.joblib']]
elif models == 'rf':
    ref_embs_directory = 'patchseq/'
    model_tups = [['RandomForest_emb_layer_genes/', 'prediction of AP width (ms) by embs max_depth 5 min_samples_leaf 2 corr 0.569 pval 0.017203570129644845 MAE 0.943.joblib'], 
                  ['RandomForest_emb_layer_genes/', 'prediction of Fitted MP (mV) by embs max_depth 5 min_samples_leaf 2 corr 0.491 pval 0.04510657863986557 MAE 18.025.joblib'], 
                  ['RandomForest_emb_layer_genes/', 'prediction of Upstroke-to-downstroke ratio by embs max_depth 5 min_samples_leaf 2 corr 0.616 pval 0.008399248652029633 MAE 1.774.joblib'], 
                  ['RandomForest_emb_layer_genes/', 'prediction of AP threshold (mV) by embs max_depth 5 min_samples_leaf 2 corr 0.799 pval 0.00011962451309732247 MAE 9.318.joblib'], 
                  ['RandomForest_emb_layer_genes/', 'prediction of Membrane time constant (ms)_log by embs max_depth 5 min_samples_leaf 4 corr 0.618 pval 0.008168682415829196 MAE 17.436.joblib'], 
                  ['RandomForest_emb_layer_genes/', 'prediction of Sag ratio_log by embs max_depth 5 min_samples_leaf 4 corr 0.235 pval 0.36431644441289895 MAE 0.111.joblib'], 
                  ['RandomForest_emb_layer_genes/', 'prediction of Rheobase (pA)_log by embs max_depth 5 min_samples_leaf 4 corr 0.23 pval 0.3755172847639521 MAE 34.353.joblib'], 
                  ['RandomForest_emb_layer_genes/', 'prediction of AP amplitude (mV) by embs max_depth 5 min_samples_leaf 2 corr 0.873 pval 4.804896750919207e-06 MAE 20.091.joblib'], 
                  ['RandomForest_emb_layer_genes/', 'prediction of Latency (ms)_log by embs max_depth 5 min_samples_leaf 2 corr 0.18 pval 0.4905910567986553 MAE 65.306.joblib'], 
                  ['RandomForest_emb_layer_genes/', 'prediction of Input resistance (MOhm)_log by embs max_depth 5 min_samples_leaf 2 corr 0.892 pval 1.5285200402965027e-06 MAE 664.917.joblib']]
elif models == 'rfhvg':
    ref_embs_directory = 'patchseq/'
    model_tups = [['RandomForest_emb_layer_hvg/', 'prediction of AP width (ms) by embs max_depth 5 min_samples_leaf 2 corr 0.696 pval 0.0019254247852281484 MAE 0.929.joblib'], 
                  ['RandomForest_emb_layer_hvg/', 'prediction of Fitted MP (mV) by embs max_depth 5 min_samples_leaf 4 corr 0.148 pval 0.5703171585586623 MAE 18.247.joblib'], 
                  ['RandomForest_emb_layer_hvg/', 'prediction of Upstroke-to-downstroke ratio by embs max_depth 5 min_samples_leaf 2 corr 0.644 pval 0.005238767406837353 MAE 1.83.joblib'], 
                  ['RandomForest_emb_layer_hvg/', 'prediction of AP threshold (mV) by embs max_depth 5 min_samples_leaf 2 corr 0.824 pval 4.704855372264906e-05 MAE 9.272.joblib'], 
                  ['RandomForest_emb_layer_hvg/', 'prediction of Membrane time constant (ms)_log by embs max_depth 5 min_samples_leaf 4 corr 0.564 pval 0.018409032356925044 MAE 16.885.joblib'], 
                  ['RandomForest_emb_layer_hvg/', 'prediction of Sag ratio_log by embs max_depth 5 min_samples_leaf 4 corr 0.115 pval 0.6616334377721342 MAE 0.105.joblib'], 
                  ['RandomForest_emb_layer_hvg/', 'prediction of Rheobase (pA)_log by embs max_depth 5 min_samples_leaf 2 corr 0.381 pval 0.13137750333080375 MAE 34.438.joblib'], 
                  ['RandomForest_emb_layer_hvg/', 'prediction of AP amplitude (mV) by embs max_depth 5 min_samples_leaf 2 corr 0.87 pval 5.661057732146265e-06 MAE 20.412.joblib'], 
                  ['RandomForest_emb_layer_hvg/', 'prediction of Latency (ms)_log by embs max_depth 5 min_samples_leaf 2 corr 0.233 pval 0.36811657431039574 MAE 66.431.joblib'], 
                  ['RandomForest_emb_layer_hvg/', 'prediction of Input resistance (MOhm) by embs max_depth 5 min_samples_leaf 4 corr 0.833 pval 3.328101411199949e-05 MAE 645.234.joblib']]
elif models == 'mlp':
    ref_embs_directory = 'patchseq/'
    model_tups = [['MLP_emb_layer_genes/', 'prediction of AP width (ms) by pcs hidden_layer_sizes 16_8 alpha 10.0 learning_rate constant corr 0.173 pval 0.5059178526351913 MAE 1.215.joblib'], 
                  ['MLP_emb_layer_genes/', 'prediction of Fitted MP (mV) by pcs hidden_layer_sizes 16 alpha 10.0 learning_rate constant corr 0.25 pval 0.3334967331231563 MAE 24.318.joblib'], 
                  ['MLP_emb_layer_genes/', 'prediction of Upstroke-to-downstroke ratio by pcs hidden_layer_sizes 16_8 alpha 10.0 learning_rate constant corr 0.629 pval 0.006850364516967192 MAE 1.63.joblib'], 
                  ['MLP_emb_layer_genes/', 'prediction of AP threshold (mV) by pcs hidden_layer_sizes 32 alpha 10.0 learning_rate constant corr 0.798 pval 0.00012170528103144876 MAE 11.65.joblib'], 
                  ['MLP_emb_layer_genes/', 'prediction of Membrane time constant (ms) by pcs hidden_layer_sizes 16 alpha 10.0 learning_rate constant corr -0.08 pval 0.7610613463273146 MAE 24.646.joblib'], 
                  ['MLP_emb_layer_genes/', 'prediction of Sag ratio by pcs hidden_layer_sizes 16_8 alpha 10.0 learning_rate constant corr 0.227 pval 0.3814177944382938 MAE 0.656.joblib'], 
                  ['MLP_emb_layer_genes/', 'prediction of Rheobase (pA)_log by pcs hidden_layer_sizes 16_8 alpha 10.0 learning_rate constant corr 0.263 pval 0.3081931797159237 MAE 35.864.joblib'], 
                  ['MLP_emb_layer_genes/', 'prediction of AP amplitude (mV) by pcs hidden_layer_sizes 16 alpha 0.1 learning_rate constant corr 0.832 pval 3.4458933591159005e-05 MAE 34.339.joblib'], 
                  ['MLP_emb_layer_genes/', 'prediction of Latency (ms) by pcs hidden_layer_sizes 16 alpha 10.0 learning_rate constant corr 0.322 pval 0.2074062938717503 MAE 62.398.joblib'], 
                  ['MLP_emb_layer_genes/', 'prediction of Input resistance (MOhm) by pcs hidden_layer_sizes 32 alpha 10.0 learning_rate constant corr 0.758 pval 0.000423149214784301 MAE 582.093.joblib']]
elif models == 'mlphvg':
    ref_embs_directory = 'patchseq/'
    model_tups = [['MLP_emb_layer_hvg/', 'prediction of AP width (ms)_log by pcs hidden_layer_sizes 16_8 alpha 10.0 learning_rate constant corr 0.273 pval 0.2898859028878062 MAE 1.019.joblib'], 
                  ['MLP_emb_layer_hvg/', 'prediction of Fitted MP (mV) by pcs hidden_layer_sizes 16 alpha 0.01 learning_rate constant corr 0.206 pval 0.42672124196430467 MAE 22.652.joblib'], 
                  ['MLP_emb_layer_hvg/', 'prediction of Upstroke-to-downstroke ratio_log by pcs hidden_layer_sizes 16_8 alpha 10.0 learning_rate constant corr 0.459 pval 0.06355053761788025 MAE 1.834.joblib'], 
                  ['MLP_emb_layer_hvg/', 'prediction of AP threshold (mV) by pcs hidden_layer_sizes 32 alpha 10.0 learning_rate constant corr 0.912 pval 3.360669885054962e-07 MAE 15.658.joblib'], 
                  ['MLP_emb_layer_hvg/', 'prediction of Membrane time constant (ms)_log by pcs hidden_layer_sizes 32 alpha 10.0 learning_rate constant corr 0.137 pval 0.6011820468581192 MAE 22.01.joblib'], 
                  ['MLP_emb_layer_hvg/', 'prediction of Sag ratio by pcs hidden_layer_sizes 16_8 alpha 10.0 learning_rate constant corr 0.145 pval 0.5794644997646701 MAE 0.264.joblib'], 
                  ['MLP_emb_layer_hvg/', 'prediction of Rheobase (pA)_log by pcs hidden_layer_sizes 16_8 alpha 10.0 learning_rate constant corr 0.083 pval 0.7516216326580749 MAE 34.872.joblib'], 
                  ['MLP_emb_layer_hvg/', 'prediction of AP amplitude (mV) by pcs hidden_layer_sizes 16 alpha 1.0 learning_rate constant corr 0.781 pval 0.0002167168396748687 MAE 35.426.joblib'], 
                  ['MLP_emb_layer_hvg/', 'prediction of Latency (ms) by pcs hidden_layer_sizes 16 alpha 0.1 learning_rate constant corr -0.016 pval 0.9501878989536846 MAE 76.593.joblib'], 
                  ['MLP_emb_layer_hvg/', 'prediction of Input resistance (MOhm) by pcs hidden_layer_sizes 32 alpha 0.1 learning_rate constant corr 0.698 pval 0.0018186877975065322 MAE 680.36.joblib']]
elif models == 'eln':
    ref_embs_directory = 'patchseq/'
    model_tups = [['ElasticNet_emb_layer_genes/', 'prediction of AP width (ms) by embs alpha 0.95 l1_ratio 0.5 corr 0.537 pval 0.026116841562111415 MAE 0.84.joblib'], 
                  ['ElasticNet_emb_layer_genes/', 'prediction of Fitted MP (mV) by embs alpha 0.95 l1_ratio 0.7 corr -0.375 pval 0.13810617989830642 MAE 26.525.joblib'], 
                  ['ElasticNet_emb_layer_genes/', 'prediction of Upstroke-to-downstroke ratio_log by embs alpha 0.25 l1_ratio 0.7 corr 0.634 pval 0.006276619927717832 MAE 1.478.joblib'], 
                  ['ElasticNet_emb_layer_genes/', 'prediction of AP threshold (mV) by embs alpha 0.95 l1_ratio 0.95 corr 0.665 pval 0.003608681238085664 MAE 10.762.joblib'], 
                  ['ElasticNet_emb_layer_genes/', 'prediction of Membrane time constant (ms)_log by embs alpha 0.95 l1_ratio 0.95 corr 0.659 pval 0.004021259623302045 MAE 18.018.joblib'], 
                  ['ElasticNet_emb_layer_genes/', 'prediction of Sag ratio_log by embs alpha 0.1 l1_ratio 0.95 corr 0.029 pval 0.9118594271781185 MAE 0.09.joblib'], 
                  ['ElasticNet_emb_layer_genes/', 'prediction of Rheobase (pA)_log by embs alpha 0.95 l1_ratio 0.7 corr 0.013 pval 0.9602692610968125 MAE 33.714.joblib'], 
                  ['ElasticNet_emb_layer_genes/', 'prediction of AP amplitude (mV) by embs alpha 0.95 l1_ratio 0.7 corr 0.796 pval 0.00013235009471175056 MAE 22.85.joblib'], 
                  ['ElasticNet_emb_layer_genes/', 'prediction of Latency (ms)_log by embs alpha 0.95 l1_ratio 0.5 corr 0.177 pval 0.49615833649609375 MAE 60.399.joblib'], 
                  ['ElasticNet_emb_layer_genes/', 'prediction of Input resistance (MOhm)_log by embs alpha 0.95 l1_ratio 0.95 corr 0.598 pval 0.011312956689847237 MAE 739.182.joblib']]
elif models == 'elnhvg':
    ref_embs_directory = 'patchseq/'
    model_tups = [['ElasticNet_emb_layer_hvg/', 'prediction of AP width (ms)_log by embs alpha 0.25 l1_ratio 0.95 corr 0.533 pval 0.02767092845817332 MAE 0.819.joblib'], 
                  ['ElasticNet_emb_layer_hvg/', 'prediction of Fitted MP (mV) by embs alpha 0.95 l1_ratio 0.95 corr -0.249 pval 0.3344651383348409 MAE 27.128.joblib'], 
                  ['ElasticNet_emb_layer_hvg/', 'prediction of Upstroke-to-downstroke ratio_log by embs alpha 0.25 l1_ratio 0.95 corr 0.623 pval 0.007494255869836565 MAE 1.523.joblib'], 
                  ['ElasticNet_emb_layer_hvg/', 'prediction of AP threshold (mV) by embs alpha 0.95 l1_ratio 0.7 corr 0.508 pval 0.037145856641156896 MAE 13.013.joblib'], 
                  ['ElasticNet_emb_layer_hvg/', 'prediction of Membrane time constant (ms)_log by embs alpha 0.95 l1_ratio 0.5 corr 0.497 pval 0.04221669340594632 MAE 16.335.joblib'], 
                  ['ElasticNet_emb_layer_hvg/', 'prediction of Sag ratio_log by embs alpha 0.1 l1_ratio 0.7 corr 0.067 pval 0.7990503708642127 MAE 0.097.joblib'], 
                  ['ElasticNet_emb_layer_hvg/', 'prediction of Rheobase (pA)_log by embs alpha 0.95 l1_ratio 0.7 corr 0.097 pval 0.7120951384508717 MAE 33.82.joblib'], 
                  ['ElasticNet_emb_layer_hvg/', 'prediction of AP amplitude (mV) by embs alpha 0.1 l1_ratio 0.5 corr 0.616 pval 0.008476006898181735 MAE 32.682.joblib'], 
                  ['ElasticNet_emb_layer_hvg/', 'prediction of Latency (ms)_log by embs alpha 0.95 l1_ratio 0.5 corr 0.366 pval 0.14863584840568128 MAE 68.125.joblib'], 
                  ['ElasticNet_emb_layer_hvg/', 'prediction of Input resistance (MOhm)_log by embs alpha 0.95 l1_ratio 0.95 corr 0.323 pval 0.20671755590153554 MAE 779.105.joblib']]
elif models == 'm1_celltype_total':
    ref_embs_directory = 'm1_patchseq_total_preds/'
    model_tups = [['LogisticRegression_emb_layer_preds/', 'prediction of RNA family by embs n 981 C 0.1 L1 0.01 F1 0.796 Acc 0.805.joblib']]
elif models == 'm1_celltype_exons':
    ref_embs_directory = 'm1_patchseq_exons_preds/'
    model_tups = [['LogisticRegression_emb_layer_preds/', 'prediction of RNA family by embs C 1.0 L1 0.1 F1 0.393 Acc 0.378.joblib']]
else:
    assert models == 'celltype'
    ref_embs_directory = 'combined_patchseq_all_preds/'
    model_tups = [['CellTypeLogisticRegression_emb_layer_scores/', 'prediction of Cell Type by embs C 3 l1_ratio 0.95 acc 0.7058823529411765.joblib']]

for model_directory, model_name in tqdm(model_tups):
    assert len(model_name.split(' by ')) == 2
    y_name = model_name.split(' by ')[0].replace('prediction of ', '')
    X_name = model_name.split(' by ')[1].split(' ')[0]
    assert X_name in ['embs', 'pcs']
    print(y_name)

    model_dir_path = ref_embs_directory + model_directory
    # scaler = load(model_dir_path + 'scaler.joblib')
    scaler_embs = load(model_dir_path + 'scaler_embs.joblib')
    scaler_pcs = load(model_dir_path + 'scaler_pcs.joblib')
    pca_name = [x for x in os.listdir(model_dir_path) if x.startswith('pca') and x.endswith('.joblib')][0]
    pca = load(model_dir_path + pca_name)
    preps_model = load(model_dir_path + model_name)
    if models in ['celltype', 'm1_celltype_total', 'm1_celltype_exons']:
        label_encoder = load(model_dir_path + 'label_encoder.joblib')

    emb_layer = model_directory.split('_')[-1].replace('/', '')
    if models in ['rf', 'rfhvg', 'mlp', 'mlphvg', 'eln', 'elnhvg']:
        assert emb_layer in ['genes', 'hvg']
        if emb_layer == 'genes':
            adata = sc.read_h5ad('glioma2/glioma2_pp_with_patchseq.h5ad')
        else:
            adata = sc.read_h5ad('glioma2/glioma2_hvg_with_patchseq.h5ad')
        adata_X = adata.X.toarray() if issparse(adata.X) else adata.X
        df_merged = pd.DataFrame(adata_X, index=adata.obs_names, columns=adata.var_names)

    else:
        assert emb_layer in ['preds', 'scores']
        embs_directory = f'{test_name}_preds/'
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

    X_test = scaler_embs.transform(df_merged.values)
    if X_name == 'pcs':
        X_test = scaler_pcs.transform(pca.transform(X_test))

    # if X_name == 'embs':
    #     X_test = df_merged.values
    # else:
    #     X_test = pca.transform(scaler.transform(df_merged.values))

    y_predict = preps_model.predict(X_test)

    if y_name.endswith('_log'):
        print('exp back')
        y_predict = np.exp(y_predict) - 1
        y_name = y_name.replace('_log', '')

    if models in ['celltype', 'm1_celltype_total', 'm1_celltype_exons']:
        y_predict = label_encoder.inverse_transform(y_predict)
        y_score = np.amax(preps_model.predict_proba(X_test), axis=1)
    
    df = pd.DataFrame({y_name: y_predict}, index=df_merged.index.tolist())
    if models in ['celltype', 'm1_celltype_total', 'm1_celltype_exons']:
        df['prob'] = y_score
    df.index.name = 'cells'
    df.to_excel(pred_output_directory + f'{y_name}.xlsx')
