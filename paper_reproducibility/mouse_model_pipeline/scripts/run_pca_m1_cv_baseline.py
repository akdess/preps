#!/usr/bin/env python3
"""Run PCA/logistic-regression baseline cross-validation on M1 Patch-seq expression."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


PROJECT_ROOT = Path("/mnt/data/serinharmanci/preps_revision_052726")
TRAINING_DIR = PROJECT_ROOT / "MOUSE_PATCHSEQ/Mouse_training_patchseq/mouse_m1_patchseq"
PACKAGE_DIR = Path(__file__).resolve().parents[1]
DEFAULT_OUT = PACKAGE_DIR / "example_outputs/pca_m1_cv"
KEEP_LABELS = ["IT", "ET", "CT", "Pvalb", "Sst", "Vip", "Lamp5", "Sncg"]
N_PCA = 50
N_FOLDS = 5
RANDOM_STATE = 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--top-genes", type=int, default=3000)
    parser.add_argument("--n-pca", type=int, default=N_PCA)
    parser.add_argument("--min-cells", type=int, default=8)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(TRAINING_DIR / "m1_patchseq_meta_data.txt", sep="\t", index_col="Cell")
    meta = meta[meta["RNA family"].isin(KEEP_LABELS)].copy()
    counts_raw = pd.read_csv(TRAINING_DIR / "m1_patchseq_exon_counts.csv", index_col=0)
    common = meta.index.intersection(counts_raw.columns)
    meta = meta.loc[common]
    labels = meta["RNA family"]
    keep = labels.value_counts()
    labels = labels[labels.isin(keep[keep >= args.min_cells].index)]
    counts = counts_raw[labels.index].T

    cpm = counts.div(counts.sum(axis=1), axis=0) * 1e6
    log_cpm = np.log1p(cpm)
    top_genes = log_cpm.var(axis=0).nlargest(args.top_genes).index
    X = log_cpm[top_genes].values
    y = labels.values

    pipe = Pipeline(
        [
            ("scaler", StandardScaler()),
            ("pca", PCA(n_components=args.n_pca, random_state=RANDOM_STATE)),
            (
                "clf",
                LogisticRegression(
                    max_iter=5000,
                    C=1.0,
                    random_state=RANDOM_STATE,
                    multi_class="multinomial",
                    n_jobs=-1,
                ),
            ),
        ]
    )
    cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    rows = []
    oof = pd.DataFrame(index=labels.index)
    oof["true_label"] = y
    oof["pred_label"] = None
    for fold, (tr, te) in enumerate(cv.split(X, y), start=1):
        pipe.fit(X[tr], y[tr])
        pred = pipe.predict(X[te])
        oof.iloc[te, oof.columns.get_loc("pred_label")] = pred
        rows.append(
            {
                "fold": fold,
                "n_train": len(tr),
                "n_test": len(te),
                "accuracy": float((pred == y[te]).mean()),
                "balanced_accuracy": float(balanced_accuracy_score(y[te], pred)),
                "macro_f1": float(f1_score(y[te], pred, average="macro")),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(args.out / "pca_logreg_m1_5fold_cv_summary.csv", index=False)
    oof.to_csv(args.out / "pca_logreg_m1_5fold_oof_predictions.csv")
    class_counts = labels.value_counts().sort_index()
    class_counts.to_csv(args.out / "pca_logreg_m1_class_counts.csv", header=["n_cells"])
    print(summary.to_string(index=False))
    print("\nmean:")
    print(summary[["accuracy", "balanced_accuracy", "macro_f1"]].mean().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
