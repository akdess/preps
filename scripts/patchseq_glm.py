import argparse
parser = argparse.ArgumentParser(description='Train prediction models on Patch-seq embeddings.')
parser.add_argument('dataset_name', help='Input the training dataset prefix (e.g., m1_patchseq).')
parser.add_argument('--meta_file', default=None, help='Custom metadata file name')
parser.add_argument('--ephys_file', default=None, help='Custom ephys features file name')
parser.add_argument('--cell_id_col', default='0', help='Name or index of the cell ID column (default: 0 for first column)')
parser.add_argument('--celltype_col', default='RNA family', help='Column name for cell types in metadata')
args = parser.parse_args()

dataset_name = args.dataset_name
meta_file = args.meta_file if args.meta_file else f'{dataset_name}_meta_data.txt'
ephys_file = args.ephys_file if args.ephys_file else f'{dataset_name}_ephys_features.csv'
cell_id_col = args.cell_id_col
celltype_col = args.celltype_col

import os
import random
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold, train_test_split
from sklearn.linear_model import ElasticNet, LogisticRegression
from scipy.sparse import issparse
from scipy.stats import pearsonr
from sklearn.metrics import confusion_matrix, f1_score, accuracy_score
from joblib import dump
from matplotlib import pyplot as plt
from tqdm import tqdm
import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)

os.environ['PYTHONHASHSEED'] = '0'
random.seed(0)
np.random.seed(0)

plt.rc('figure', figsize=(6, 6))
plt.rc('font', size=10)

# ---------------------------------------------------------
# Helper: Robustly load DataFrames regardless of format
# ---------------------------------------------------------
def load_df(file_path, cell_id_col):
    # engine='python' and sep=None auto-detects delimiter (comma, tab, etc.)
    df = pd.read_csv(file_path, sep=None, engine='python')
    if cell_id_col.isdigit():
        col_idx = int(cell_id_col)
        if col_idx < len(df.columns):
            df.set_index(df.columns[col_idx], inplace=True)
    else:
        if cell_id_col in df.columns:
            df.set_index(cell_id_col, inplace=True)
        else:
            print(f"⚠️ Warning: Column '{cell_id_col}' not found in {file_path}. Defaulting to first column.")
            df.set_index(df.columns[0], inplace=True)
    return df

# ---------------------------------------------------------
# Load and preprocess Data
# ---------------------------------------------------------
# Meta data
df_meta = load_df(meta_file, cell_id_col)

if celltype_col not in df_meta.columns:
    raise ValueError(f"Cell type column '{celltype_col}' not found in metadata.")
    
df_meta = df_meta[[celltype_col]].dropna()
df_meta = df_meta[df_meta[celltype_col].map(lambda x: str(x).lower() != 'low quality')]
df_meta = df_meta.groupby(celltype_col).filter(lambda x: len(x) >= 10)
label_encoder = LabelEncoder()
df_meta[celltype_col] = label_encoder.fit_transform(df_meta[celltype_col].values)

# Ephys features
df_ephys = load_df(ephys_file, cell_id_col)
ephys_list = ['Input resistance (MOhm)', 'Latency (ms)', 'AP amplitude (mV)', 'Rheobase (pA)', 'Sag ratio', 
              'Membrane time constant (ms)', 'AP threshold (mV)', 'Upstroke-to-downstroke ratio', 'ISI adaptation index', 'AP width (ms)']

# Dynamically keep ephys features that actually exist in the dataframe to prevent crashes
available_ephys = [col for col in ephys_list if col in df_ephys.columns]
if not available_ephys:
    available_ephys = df_ephys.select_dtypes(include=[np.number]).columns.tolist()
    print(f"Standard ephys columns not found. Using available numeric columns: {available_ephys}")
    
df_ephys = df_ephys[available_ephys].dropna()
print(f"Cells remaining after strict ephys dropna: {df_ephys.shape[0]}")
df_ephys = df_ephys[df_ephys.index.isin(df_meta.index.tolist())]

# ---------------------------------------------------------
# Process over both embedding layers ('preds' and 'scores')
# ---------------------------------------------------------
embs_directory = f'{dataset_name}_preds/'

for emb_layer in ['preds', 'scores']:
    print(f"\n{'='*60}\n🚀 Processing embedding layer: {emb_layer}\n{'='*60}")
    embs_files = [x for x in os.listdir(embs_directory) if x.endswith(f'{emb_layer}.csv')]
    embs_files = sorted(embs_files)
    
    if not embs_files:
        print(f"⚠️ No '{emb_layer}.csv' files found in {embs_directory}. Skipping.")
        continue

    # Output folders
    elnet_output_dir = f'{embs_directory}ElasticNet_emb_layer_{emb_layer}/'
    logreg_output_dir = f'{embs_directory}LogisticRegression_emb_layer_{emb_layer}/'
    os.makedirs(elnet_output_dir, exist_ok=True)
    os.makedirs(logreg_output_dir, exist_ok=True)

    # Load embeddings
    df_merged = None
    for file_name in embs_files:
        df = pd.read_csv(embs_directory + file_name, sep=',', index_col='individual')
        df = df[df.index.isin(df_meta.index.tolist())]
        df = df[df.index.isin(df_ephys.index.tolist())]
        
        # Safely drop metadata columns from embeddings if present
        cols_to_drop = [c for c in ['Unnamed: 0', 'group'] if c in df.columns]
        if cols_to_drop:
            df = df.drop(columns=cols_to_drop)
            
        if emb_layer == 'features' or emb_layer == 'scores':
            df = df.iloc[:, :-2]
        
        df.columns = [f"{file_name.replace('.csv', '')}_dim_{x}" for x in df.columns.tolist()]
        
        if df_merged is None:
            df_merged = df.copy()
        else:
            df_merged = pd.merge(df_merged, df, how='inner', left_index=True, right_index=True)
            
    if df_merged is None or df_merged.empty:
        print("⚠️ No valid intersecting cells found. Skipping this layer.")
        continue

    # Align features BEFORE splitting so indices match perfectly
    common_cells = df_merged.index.intersection(df_meta.index).intersection(df_ephys.index)
    
    if len(common_cells) < 10:
         print(f"⚠️ Not enough intersecting cells ({len(common_cells)}). Skipping this layer.")
         continue
         
    df_meta_sub = df_meta.loc[common_cells, :]
    df_ephys_sub = df_ephys.loc[common_cells, :]
    df_merged_sub = df_merged.loc[common_cells, :]
    
    X_all = df_merged_sub.values
    y_stratify = df_meta_sub[celltype_col].values

    X_train, X_test, cells_train, cells_test = train_test_split(
        X_all, common_cells.values, 
        test_size=0.2, 
        random_state=0, 
        stratify=y_stratify 
    )

    # Split targets for Regression (Ephys) vs Classification (Cell type)
    y_dict_ephys = {}
    for ephy in available_ephys:
        y_train = df_ephys_sub.loc[cells_train, ephy].values
        y_test = df_ephys_sub.loc[cells_test, ephy].values
        y_dict_ephys[ephy] = [y_train, y_test]
        if np.amin(y_train) >= 0:
            y_dict_ephys[f'{ephy}_log'] = [np.log(y_train + 1), y_test]

    y_dict_class = {}
    y_train_class = df_meta_sub.loc[cells_train, celltype_col].values
    y_test_class = df_meta_sub.loc[cells_test, celltype_col].values
    y_dict_class[celltype_col] = [y_train_class, y_test_class]

    # Fit Scalers and PCA
    scaler_embs = StandardScaler()
    embs_train = scaler_embs.fit_transform(X_train)
    embs_test = scaler_embs.transform(X_test)
    
    X_dict = {}
    X_dict['embs'] = [embs_train, embs_test]

    n_components = min(0.85, len(embs_train) - 1, embs_train.shape[1]) if isinstance(0.85, float) else 0.85
    pca = PCA(n_components=n_components, random_state=0)
    pcs_train = pca.fit_transform(embs_train)
    pcs_test = pca.transform(embs_test)
    num_pcs_used = pca.n_components_

    scaler_pcs = StandardScaler()
    pcs_train = scaler_pcs.fit_transform(pcs_train)
    pcs_test = scaler_pcs.transform(pcs_test)
    X_dict['pcs'] = [pcs_train, pcs_test]

    # Save scaler and pca models to BOTH directories
    for out_dir in [elnet_output_dir, logreg_output_dir]:
        dump(scaler_embs, out_dir + 'scaler_embs.joblib')
        dump(scaler_pcs, out_dir + 'scaler_pcs.joblib')
        dump(pca, out_dir + f'pca_{num_pcs_used}.joblib')
        dump(label_encoder, out_dir + 'label_encoder.joblib')


    # ---------------------------------------------------------
    # 1. Fit Elastic Net model for each ephys feature
    # ---------------------------------------------------------
    print("📈 Training Elastic Net Regressors for continuous Ephys Features...")
    for y_name, y_tup in tqdm(y_dict_ephys.items()):
        y_train, y_test = y_tup

        for X_name, X_tup in X_dict.items():
            X_train_fold, X_test_fold = X_tup

            model = ElasticNet(max_iter=10000, selection='random', random_state=0)
            grid = {'alpha': [0.0001, 0.001, 0.01, 0.1, 1.0], 'l1_ratio': [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]}
            
            cv = KFold(n_splits=min(10, len(X_train_fold)), shuffle=True, random_state=0)
            search = GridSearchCV(model, grid, scoring='neg_mean_absolute_error', cv=cv, n_jobs=-1)
            results = search.fit(X_train_fold, y_train)
            alpha = results.best_params_['alpha']
            l1_ratio = results.best_params_['l1_ratio']

            model = results.best_estimator_
            y_predict = model.predict(X_test_fold)

            if y_name.endswith('_log'):
                y_predict = np.exp(y_predict) - 1

            try:
                corr, pval = pearsonr(y_predict, y_test)
                corr = np.round(corr, 3)
            except:
                corr, pval = 'NA', 'NA'
            
            mae = np.mean(np.abs(y_predict - y_test))
            mae = np.round(mae, 3) if np.round(mae, 3) > 0.01 else np.round(mae, 8)
            
            file_prefix = f'prediction of {y_name} by {X_name} alpha {alpha} l1 {l1_ratio} corr {corr} p {pval} MAE {mae}'
            dump(model, elnet_output_dir + f'{file_prefix}.joblib')
            
            plot_title = (
                f"Prediction: {y_name} by {X_name}\n"
                f"Params: alpha={alpha}, l1_ratio={l1_ratio}\n"
                f"Performance: r={corr} (p={pval}), MAE={mae}"
            )
            
            plt.figure()
            plt.scatter(y_test, y_predict)
            plt.title(plot_title, fontsize=10)
            clean_y_name = y_name.replace('_log', '')
            plt.xlabel(f'{clean_y_name} for test')
            plt.ylabel(f'{clean_y_name} by prediction')
            figure = plt.gcf()
            figure.patch.set_facecolor('white')
            figure.savefig(elnet_output_dir + f'{file_prefix}.pdf', bbox_inches='tight', dpi=300)
            plt.close('all')

            df_plot = pd.DataFrame({'cells_test': cells_test, 'y_test': y_test, 'y_predict': y_predict}).set_index('cells_test')
            df_plot.to_csv(elnet_output_dir + f'{file_prefix}.csv', sep=',')


    # ---------------------------------------------------------
    # 2. Fit Logistic Regression model for cell type classification
    # ---------------------------------------------------------
    print("📊 Training Logistic Regression Classifiers for categorical Cell Types...")
    for y_name, y_tup in tqdm(y_dict_class.items()):
        y_train, y_test = y_tup

        for X_name, X_tup in X_dict.items():
            X_train_fold, X_test_fold = X_tup

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

            cv = StratifiedKFold(n_splits=min(5, len(np.unique(y_train))), shuffle=True, random_state=0)
            search = GridSearchCV(model, grid, scoring='f1_weighted', cv=cv, n_jobs=-1)
            results = search.fit(X_train_fold, y_train)
            
            c = results.best_params_['C']
            l1_ratio = results.best_params_['l1_ratio']
            model = results.best_estimator_

            y_predict = model.predict(X_test_fold)
            f1 = f1_score(y_test, y_predict, average='weighted')
            acc = accuracy_score(y_test, y_predict)
            
            f1 = np.round(f1, 3)
            acc = np.round(acc, 3)

            file_prefix = f'prediction of {y_name} by {X_name} C {c} L1 {l1_ratio} F1 {f1} Acc {acc}'
            dump(model, logreg_output_dir + f'{file_prefix}.joblib')

            y_predict_label = label_encoder.inverse_transform(y_predict)
            y_test_label = label_encoder.inverse_transform(y_test)

            df_plot = pd.DataFrame({'cells_test': cells_test, 'y_test': y_test_label, 'y_predict': y_predict_label}).set_index('cells_test')
            df_plot.to_csv(logreg_output_dir + f'{file_prefix}.csv', sep=',')

            try:
                labels = label_encoder.classes_
                conf_array = confusion_matrix(y_true=y_test_label, y_pred=y_predict_label, labels=labels)
                assert np.sum(conf_array) == len(y_test_label)
                conf_df = pd.DataFrame(conf_array, index=[f'{x}_true' for x in labels], columns=[f'{x}_predicted' for x in labels])
                conf_df.index.name = 'confusion_matrix'
                conf_df.to_csv(logreg_output_dir + f'{file_prefix}_confusion_matrix.csv', sep=',')
            except Exception as e:
                print(f'confusion_matrix not generated: {e}')

print("\n✅ Training Complete for both ElasticNet Regressors and Logistic Regression Classifiers!")