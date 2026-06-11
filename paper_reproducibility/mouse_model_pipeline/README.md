# Mouse Model Pipeline

This folder contains the mouse cell-type prediction and ephys feature prediction model pipeline used for the fine-grained mouse cortical subclass figures (the "previous-models" / 304-dim package used for Figure 5).

## Contents

- `models/celltype_classifier/`
  - Saved PREPS cell-type classifier models.
  - Main model: `logreg__embs.joblib`.

- `models/ephys_regressors/`
  - Saved PREPS ephys feature regressor models.
  - Main models: `*__embs.joblib`.

- `models/mouse_geneformer_finetuned/`
  - The fine-tuned Mouse-Geneformer model that generated the 304-dim bridge embeddings in
    `embeddings/`. Fine-tuned on the Allen WMB Isocortex reference with 24 subclass labels
    (`target_names.csv`), 2,000 cells per class, 10 epochs, seed 0 (`finetune_config.json`).
    Held-out eval: accuracy 0.963, balanced 0.897, macro-F1 0.916 (`eval_metrics.json`).
  - Files: `config.json`, `model.safetensors`, `target_names.csv`, `finetune_config.json`,
    `eval_metrics.json`, `training_args.bin`.
  - Bridge embedding = 256 hidden + 24 logits + 24 probs = 304 dims (matches the saved
    classifier/regressor `embedding_dim`).

- `embeddings/`
  - WMB and GSE185862 bridge embedding files used by the saved models:
    - `m1_patchseq_embeddings.npz`
    - `wmb_isocortex_1_embeddings.npz`
    - `wmb_isocortex_2_embeddings.npz`
    - `gse185862_ssv4_embeddings.npz`

- `scripts/predict_wmb_gse.py`
  - Applies the saved cell-type classifier and ephys regressors to WMB and GSE185862.

- `scripts/train_mouse_models.py`
  - Training counterpart for the saved models. It trains cell-type classifiers and ElasticNetCV ephys regressors from M1 Patch-seq bridge embeddings, writes `.joblib` bundles, and records CV summaries.

- `scripts/plot_gse185862_main_figure.R`
  - Generates the GSE185862 2x2 concordance/ephys figure.

- `scripts/plot_pca_baseline_wmb_concordance.R`
  - Generates the PCA baseline WMB concordance figure.

- `scripts/run_pca_m1_cv_baseline.py`
  - Runs PCA/logistic-regression 5-fold M1 Patch-seq cross-validation from raw M1 expression.

- `figures/`
  - Main generated figure files copied into this folder.

- `results/`
  - Source and summary tables needed for the plots.

## Generate Embeddings From the Fine-Tuned Model

The bridge embeddings in `embeddings/` were produced by `models/mouse_geneformer_finetuned/`:
for each tokenized cell, take the 256-dim mean hidden state and concatenate the 24 class logits
and 24 softmax probabilities from the fine-tuned classification head (256 + 24 + 24 = 304).

## Run WMB and GSE Predictions

From the repository root:

```bash
/usr/bin/python analysis_fine_grained_mouse_cortical_subclasses/mouse_model_pipeline/scripts/predict_wmb_gse.py \
  --dataset both \
  --classifier logreg \
  --mode embs \
  --out analysis_fine_grained_mouse_cortical_subclasses/mouse_model_pipeline/example_outputs/predictions
```

## Retrain the Mouse Cell-Type and Ephys Models

This retrains from `embeddings/m1_patchseq_embeddings.npz` and writes new models under `trained_models_from_code/`.

```bash
/usr/bin/python analysis_fine_grained_mouse_cortical_subclasses/mouse_model_pipeline/scripts/train_mouse_models.py
```

Cell-type model selection trains `logreg`, `randomforest`, and `extratrees` on both `embs` and `pcs` and writes `summary.csv`. Ephys model selection trains one `ElasticNetCV` per feature and mode; each `.joblib` stores the selected `alpha_` and `l1_ratio_`.

## Generate Main Plots

```bash
analysis_fine_grained_mouse_cortical_subclasses/mouse_model_pipeline/scripts/generate_main_plots.sh
```

This runs:

```bash
Rscript analysis_fine_grained_mouse_cortical_subclasses/mouse_model_pipeline/scripts/plot_gse185862_main_figure.R
Rscript analysis_fine_grained_mouse_cortical_subclasses/mouse_model_pipeline/scripts/plot_pca_baseline_wmb_concordance.R
```

The main figure files are included in `figures/`:

```text
figures/Figure5_WMB_concordance_ephys.png
figures/Figure5_WMB_concordance_ephys.pdf
figures/Figure_GSE185862_concordance_ephys.png
figures/Figure_GSE185862_concordance_ephys.pdf
figures/Figure_WMB_PCA_baseline_concordance.png
figures/Figure_WMB_PCA_baseline_concordance.pdf
```

## PCA Baseline

```bash
/usr/bin/python analysis_fine_grained_mouse_cortical_subclasses/mouse_model_pipeline/scripts/run_pca_m1_cv_baseline.py
```

PCA WMB prediction outputs are included under `pca_baseline/results/`.

## Provenance

See `manifest.json` for file counts, sizes, and SHA256 checksums.
