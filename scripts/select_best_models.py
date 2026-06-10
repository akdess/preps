import os
import shutil
import glob
import argparse

def select_best_models(base_dir):
    print(f"\nEvaluating models across all subdirectories in: {base_dir}")
    
    # 1. Glob all joblib files recursively in the base directory
    search_pattern = os.path.join(base_dir, '**', 'prediction of *.joblib')
    model_files = glob.glob(search_pattern, recursive=True)
    
    # Exclude previously prefixed 'best__' files if they exist to avoid double counting
    model_files = [f for f in model_files if 'best__' not in os.path.basename(f)]
    
    if not model_files:
        print("⚠️ No prediction model files found.")
        return

    # 2. Group models by feature
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
        elif 'acc' in basename:  # Handle lower-case 'acc' from older scripts
            metric_type = 'Acc'
            metric_val = float(basename.split('acc ')[-1])
        else:
            continue
        
        # Extract the specific feature or cell type being predicted
        feature_name = basename.split(' by ')[0].replace('prediction of ', '')
        
        if feature_name not in feature_dict:
            feature_dict[feature_name] = []
        feature_dict[feature_name].append((metric_val, f, metric_type))
        
    # Clean up any old 'best__' files before copying new ones
    old_bests = glob.glob(os.path.join(base_dir, '**', 'best__*.joblib'), recursive=True)
    for old in old_bests:
        os.remove(old)

    # 3. Find and duplicate the absolute best model for each feature
    for feature, models in feature_dict.items():
        metric_type = models[0][2]
        if metric_type == 'MAE':
            best_model = min(models, key=lambda x: x[0]) # Minimize MAE
        else:
            best_model = max(models, key=lambda x: x[0]) # Maximize Accuracy
            
        best_file = best_model[1]
        target_dir = os.path.dirname(best_file)
        
        # Prefix the original filename with 'best__' to preserve metadata
        dest_file = os.path.join(target_dir, f'best__{os.path.basename(best_file)}')
        
        shutil.copy(best_file, dest_file)
        print(f"✅ Best {feature} ({metric_type}: {best_model[0]}): -> Saved in {os.path.basename(target_dir)}/{os.path.basename(dest_file)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Select best prediction models across embedding layers.')
    parser.add_argument('embs_dir', help='Base directory containing model subfolders (e.g., m1_patchseq_exons_preds/)')
    args = parser.parse_args()
    
    if os.path.exists(args.embs_dir):
        select_best_models(args.embs_dir)
    else:
        print(f"❌ Directory not found: {args.embs_dir}")