#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/serinharmanci/preps_revision_052726"

cd "${ROOT}"

# GSE185862 main/supplementary 2x2 figure.
Rscript "analysis_fine_grained_mouse_cortical_subclasses/mouse_model_pipeline/scripts/plot_gse185862_main_figure.R"

# PCA baseline WMB concordance supplement.
Rscript "analysis_fine_grained_mouse_cortical_subclasses/mouse_model_pipeline/scripts/plot_pca_baseline_wmb_concordance.R"

echo "[done] Main plot scripts completed."
