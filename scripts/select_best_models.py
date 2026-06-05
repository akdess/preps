import os
import shutil
import glob

def select_best_models(model_dir):
    print(f"\nEvaluating models in: {model_dir}")
    joblib_files = glob.glob(os.path.join(model_dir, '*.joblib'))
    
    # Isolate the actual prediction models (exclude scalers, PCA, etc.)
    model_files = [f for f in joblib_files if 'prediction of' in os.path.basename(f)]
    
    feature_dict = {}
    for f in model_files:
        basename = os.path.basename(f).replace('.joblib', '')
        
        # Determine metric based on file name structure
        if 'MAE' in basename:
            metric_type = 'MAE'
            metric_val = float(basename.split('MAE ')[-1])
        elif 'Acc' in basename:
            metric_type = 'Acc'
            metric_val = float(basename.split('Acc ')[-1])
        else:
            continue
        
        # Extract the specific feature or cell type being predicted
        feature_name = basename.split(' by ')[0].replace('prediction of ', '')
        
        if feature_name not in feature_dict:
            feature_dict[feature_name] = []
        feature_dict[feature_name].append((metric_val, f, metric_type))
        
    # Find and copy the best model for each feature
    for feature, models in feature_dict.items():
        metric_type = models[0][2]
        if metric_type == 'MAE':
            best_model = min(models, key=lambda x: x[0]) # Minimize MAE
        else:
            best_model = max(models, key=lambda x: x[0]) # Maximize Accuracy
            
        best_file = best_model[1]
        
        # Create a standardized name for downstream prediction
        clean_feature_name = feature.replace("/", "_").replace(" ", "_")
        dest_file = os.path.join(model_dir, f'best_model_{clean_feature_name}.joblib')
        
        shutil.copy(best_file, dest_file)
        print(f"✅ Best {feature}: {os.path.basename(best_file)} -> Saved as best_model_{clean_feature_name}.joblib")

if __name__ == "__main__":
    # Add all directories where patchseq_glm.py saves models
    directories_to_scan = [
        'm1_patchseq_exons_preds/ElasticNet_emb_layer_preds/',
        'm1_patchseq_exons_preds/LogisticRegression_emb_layer_preds/'
    ]
    for d in directories_to_scan:
        if os.path.exists(d):
            select_best_models(d)