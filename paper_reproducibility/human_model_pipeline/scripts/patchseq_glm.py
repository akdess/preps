import os
import random
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold, train_test_split
from sklearn.linear_model import ElasticNet, LogisticRegression
# from sklearn.ensemble import RandomForestRegressor
# from sklearn.neural_network import MLPRegressor
from scipy.sparse import issparse
from scipy.stats import pearsonr
from sklearn.metrics import confusion_matrix, f1_score, accuracy_score
from joblib import dump
from matplotlib import pyplot as plt
from tqdm import tqdm


os.environ['PYTHONHASHSEED'] = '0'
random.seed(0)
np.random.seed(0)

plt.rc('figure', figsize=(6, 6))
plt.rc('font', size=10)

# Load ephys features
# df_meta = pd.read_csv('m1_patchseq/m1_patchseq_ephys_features.csv', index_col='cell id')
df_meta = pd.read_csv('m1_patchseq_meta_data.txt', sep='\t', index_col='Cell')
# ephys_list = ['Input resistance (MOhm)', 'Latency (ms)', 'AP amplitude (mV)', 'Rheobase (pA)', 'Sag ratio', 
#               'Membrane time constant (ms)', 'AP threshold (mV)', 'Upstroke-to-downstroke ratio', 'ISI adaptation index', 'AP width (ms)']
celltype_col = 'RNA family'
ephys_list = [celltype_col]
df_meta = df_meta[ephys_list].dropna()
df_meta = df_meta[df_meta[celltype_col].map(lambda x: x != 'low quality')]
df_meta = df_meta.groupby(celltype_col).filter(lambda x: len(x) >= 10)
label_encoder = LabelEncoder()
df_meta[celltype_col] = label_encoder.fit_transform(df_meta[celltype_col].values)

# Select embedding files
embs_directory = 'm1_patchseq_exons_preds/'
emb_layer = ['preds', 'scores'][1]
embs_files = [x for x in os.listdir(embs_directory) if x.endswith(f'{emb_layer}.csv')]
embs_files = sorted(embs_files)

# Create output folder
output_directory = f'{embs_directory}LogisticRegression_emb_layer_{emb_layer}/'
os.mkdir(output_directory)

# Load embeddings
df_merged = None
for file_name in embs_files:
    df = pd.read_csv(embs_directory + file_name, sep=',', index_col='individual')
    df = df[df.index.isin(df_meta.index.tolist())]
    
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

# Align ephys features BEFORE splitting so indices match perfectly
df_meta = df_meta.loc[df_merged.index.tolist(), :]
X_all = df_merged.values
y_stratify = df_meta[celltype_col].values

X_train, X_test, cells_train, cells_test = train_test_split(
    X_all, df_meta.index.values, 
    test_size=0.2, 
    random_state=0, 
    stratify=y_stratify 
)

X_dict = {}
y_dict = {}
for ephy in ephys_list:
    y_train = df_meta.loc[cells_train, ephy].values
    y_test = df_meta.loc[cells_test, ephy].values
    y_dict[ephy] = [y_train, y_test] # original y_train and y_test
    # if np.amin(y_train) >= 0:
    #     y_dict[f'{ephy}_log'] = [np.log(y_train + 1), y_test] # log y_train, keep y_test

# Fit Scalers and PCA on the new perfectly stratified X_train
scaler_embs = StandardScaler()
embs_train = scaler_embs.fit_transform(X_train)
embs_test = scaler_embs.transform(X_test)
X_dict['embs'] = [embs_train, embs_test]

pca = PCA(n_components=0.85, random_state=0)
pcs_train = pca.fit_transform(embs_train)
pcs_test = pca.transform(embs_test)
num_pcs_used = pca.n_components_

scaler_pcs = StandardScaler()
pcs_train = scaler_pcs.fit_transform(pcs_train)
pcs_test = scaler_pcs.transform(pcs_test)
X_dict['pcs'] = [pcs_train, pcs_test]

# Save scaler and pca models
dump(scaler_embs, output_directory + 'scaler_embs.joblib')
dump(scaler_pcs, output_directory + 'scaler_pcs.joblib')
dump(pca, output_directory + f'pca_{num_pcs_used}.joblib')
dump(label_encoder, output_directory + 'label_encoder.joblib')


# # Fit elastic net model for each ephys feature
# for y_name, y_tup in tqdm(y_dict.items()):
#     y_train, y_test = y_tup

#     for X_name, X_tup in X_dict.items():
#         X_train, X_test = X_tup

#         # if X_name == 'pcs':
#         #     continue

#         model = ElasticNet(max_iter=10000, selection='random', random_state=0)
#         grid = {'alpha': [0.0001, 0.001, 0.01, 0.1, 1.0], 'l1_ratio': [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]}
        
#         # cv = RepeatedKFold(n_splits=10, n_repeats=10, random_state=0)
#         cv = KFold(n_splits=10, shuffle=True, random_state=0)
#         search = GridSearchCV(model, grid, scoring='neg_mean_absolute_error', cv=cv, n_jobs=-1)
#         results = search.fit(X_train, y_train)
#         alpha = results.best_params_['alpha']
#         l1_ratio = results.best_params_['l1_ratio']

#         model = results.best_estimator_
#         y_predict = model.predict(X_test)

#         if y_name.endswith('_log'):
#             y_predict = np.exp(y_predict) - 1

#         try:
#             corr, pval = pearsonr(y_predict, y_test)
#             corr = np.round(corr, 3)
#         except:
#             corr, pval = 'NA', 'NA'
        
#         mae = np.mean(np.abs(y_predict - y_test))
#         mae = np.round(mae, 3) if np.round(mae, 3) > 0.01 else np.round(mae, 8)
        
#         file_prefix = f'prediction of {y_name} by {X_name} n {n_train} alpha {alpha} l1 {l1_ratio} corr {corr} p {pval} MAE {mae}'
#         dump(model, output_directory + f'{file_prefix}.joblib')
        
#         plot_title = (
#             f"Prediction: {y_name} by {X_name}\n"
#             f"Params: alpha={alpha}, l1_ratio={l1_ratio}\n"
#             f"Performance: r={corr} (p={pval}), MAE={mae}"
#         )
        
#         plt.figure()
#         plt.scatter(y_test, y_predict)
#         plt.title(plot_title, fontsize=10)
#         plt.xlabel(f'{y_name} for test')
#         plt.ylabel(f'{y_name} by prediction')
#         figure = plt.gcf()
#         figure.patch.set_facecolor('white')
#         figure.savefig(output_directory + f'{file_prefix}.pdf', bbox_inches='tight', dpi=300)
#         plt.close('all')

#         df_plot = pd.DataFrame({'cells_test': cells_test, 'y_test': y_test, 'y_predict': y_predict}).set_index('cells_test')
#         df_plot.to_csv(output_directory + f'{file_prefix}.csv', sep=',')

# Fit logistic regression model for cell type classification
for y_name, y_tup in tqdm(y_dict.items()):
    y_train, y_test = y_tup

    for X_name, X_tup in X_dict.items():
        X_train, X_test = X_tup

        model = LogisticRegression(
            penalty='elasticnet', 
            solver='saga', 
            multi_class='multinomial', 
            max_iter=5000, 
            random_state=0, 
            class_weight='balanced'
        )

        grid = {
            'C': [0.01, 0.1, 1.0, 10.0], 
            'l1_ratio': [0.01, 0.1, 0.5, 0.9, 1.0]
        }

        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
        search = GridSearchCV(model, grid, scoring='f1_weighted', cv=cv, n_jobs=-1)
        results = search.fit(X_train, y_train)
        
        c = results.best_params_['C']
        l1_ratio = results.best_params_['l1_ratio']
        model = results.best_estimator_

        y_predict = model.predict(X_test)
        f1 = f1_score(y_test, y_predict, average='weighted')
        acc = accuracy_score(y_test, y_predict)
        
        f1 = np.round(f1, 3)
        acc = np.round(acc, 3)

        file_prefix = f'prediction of {y_name} by {X_name} C {c} L1 {l1_ratio} F1 {f1} Acc {acc}'
        dump(model, output_directory + f'{file_prefix}.joblib')

        y_predict_label = label_encoder.inverse_transform(y_predict)
        y_test_label = label_encoder.inverse_transform(y_test)

        df_plot = pd.DataFrame({'cells_test': cells_test, 'y_test': y_test_label, 'y_predict': y_predict_label}).set_index('cells_test')
        df_plot.to_csv(output_directory + f'{file_prefix}.csv', sep=',')

        try:
            labels = label_encoder.classes_
            conf_array = confusion_matrix(y_true=y_test_label, y_pred=y_predict_label, labels=labels)
            assert np.sum(conf_array) == len(y_test_label)
            conf_df = pd.DataFrame(conf_array, index=[f'{x}_true' for x in labels], columns=[f'{x}_predicted' for x in labels])
            conf_df.index.name = 'confusion_matrix'
            conf_df.to_csv(output_directory + f'{file_prefix} confusion_matrix.csv', sep=',')
        except:
            print('confusion_matrix not generated')
