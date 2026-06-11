#!/usr/bin/env python3
"""Train mouse PREPS cell-type classifiers and ephys regressors.

This is the training counterpart of `predict_wmb_gse.py`.

Cell-type model selection:
  - trains logreg/randomforest/extratrees on `embs` and `pcs`
  - evaluates by 5-fold stratified CV
  - saves every trained final model plus `summary.csv`
  - the main pipeline uses the best CV model (`logreg__embs` in the saved run)

Ephys model selection:
  - trains one ElasticNetCV per ephys feature and input mode
  - ElasticNetCV chooses alpha and l1_ratio internally
  - evaluates by 5-fold out-of-fold prediction
  - saves final full-data models plus `summary.csv`
  - the main pipeline uses `*__embs.joblib`
"""

from __future__ import annotations

import argparse
import json
import os
import random
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import dump
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.linear_model import ElasticNetCV, LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, mean_absolute_error, r2_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder, StandardScaler


ROOT = Path("/mnt/data/serinharmanci/preps_revision_052726")
PIPELINE = Path(__file__).resolve().parents[1]
DEFAULT_EMB = PIPELINE / "embeddings/m1_patchseq_embeddings.npz"
DEFAULT_OUT = PIPELINE / "trained_models_from_code"
META = ROOT / "MOUSE_PATCHSEQ/Mouse_training_patchseq/mouse_m1_patchseq/m1_patchseq_meta_data.txt"
EPHYS = ROOT / "MOUSE_PATCHSEQ/Mouse_training_patchseq/mouse_m1_patchseq/m1_patchseq_ephys_features.csv"

SEED = 0
N_FOLDS = 5
PCA_DIM = 10
L1_RATIOS = [0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 1.0]
N_ALPHAS = 100
INNER_CV = 5
MAX_ITER = 10000
DROP_LABELS = {"low quality", "NP"}

EPHYS_FEATURES = [
    "Input resistance (MOhm)",
    "Latency (ms)",
    "AP amplitude (mV)",
    "Rheobase (pA)",
    "Sag ratio",
    "Membrane time constant (ms)",
    "AP threshold (mV)",
    "Upstroke-to-downstroke ratio",
    "Resting membrane potential (mV)",
    "AP width (ms)",
    "Afterhyperpolarization (mV)",
    "ISI adaptation index",
    "Max number of APs",
]
LOG_FEATURES = {
    "Input resistance (MOhm)",
    "Latency (ms)",
    "AP width (ms)",
    "Membrane time constant (ms)",
    "Rheobase (pA)",
    "Sag ratio",
    "Max number of APs",
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def seed() -> None:
    os.environ["PYTHONHASHSEED"] = str(SEED)
    random.seed(SEED)
    np.random.seed(SEED)


def sanitize(x: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", x)).strip("_")


def load_embeddings(path: Path) -> pd.DataFrame:
    arr = np.load(path, allow_pickle=True)
    df = pd.DataFrame(np.asarray(arr["embeddings"], dtype=np.float64), index=arr["cell_ids"])
    df.columns = [f"emb_{i:03d}" for i in range(df.shape[1])]
    df.index.name = "cell_id"
    return df


def load_labels(min_cells: int) -> pd.Series:
    meta = pd.read_csv(META, sep="\t")
    meta = meta[["Cell", "RNA family"]].dropna()
    meta = meta[~meta["RNA family"].isin(DROP_LABELS)]
    counts = meta["RNA family"].value_counts()
    meta = meta[meta["RNA family"].isin(counts[counts >= min_cells].index)]
    return meta.set_index("Cell")["RNA family"]


def load_ephys() -> pd.DataFrame:
    df = pd.read_csv(EPHYS, index_col="cell id")
    return df[[f for f in EPHYS_FEATURES if f in df.columns]]


def make_transform(mode: str, X_train: np.ndarray, X_test: np.ndarray | None = None):
    scaler = StandardScaler()
    X_tr = scaler.fit_transform(X_train)
    X_te = scaler.transform(X_test) if X_test is not None else None
    pca = None
    if mode == "pcs":
        pca = PCA(n_components=PCA_DIM, random_state=SEED)
        X_tr = pca.fit_transform(X_tr)
        X_te = pca.transform(X_te) if X_te is not None else None
    return scaler, pca, X_tr, X_te


def clf_factory(name: str):
    if name == "logreg":
        return LogisticRegression(max_iter=5000, solver="lbfgs", n_jobs=-1, random_state=SEED)
    if name == "randomforest":
        return RandomForestClassifier(n_estimators=500, n_jobs=-1, random_state=SEED)
    if name == "extratrees":
        return ExtraTreesClassifier(n_estimators=500, n_jobs=-1, random_state=SEED)
    raise ValueError(name)


def train_celltype(embs: pd.DataFrame, labels: pd.Series, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    common = embs.index.intersection(labels.index)
    X = embs.loc[common].values
    label_values = labels.loc[common].values
    encoder = LabelEncoder()
    y = encoder.fit_transform(label_values)
    cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    rows = []
    for mode in ["embs", "pcs"]:
        for clf_name in ["logreg", "randomforest", "extratrees"]:
            log(f"[celltype] {clf_name}__{mode}")
            oof = np.full(len(y), -1, dtype=int)
            for tr, te in cv.split(X, y):
                _, _, X_tr, X_te = make_transform(mode, X[tr], X[te])
                clf = clf_factory(clf_name)
                clf.fit(X_tr, y[tr])
                oof[te] = clf.predict(X_te)
            acc = accuracy_score(y, oof)
            bal = balanced_accuracy_score(y, oof)
            macro = f1_score(y, oof, average="macro")
            scaler, pca, X_full, _ = make_transform(mode, X)
            final = clf_factory(clf_name)
            final.fit(X_full, y)
            dump(
                {
                    "classifier": final,
                    "scaler": scaler,
                    "pca": pca,
                    "label_encoder": encoder,
                    "classifier_name": clf_name,
                    "mode": mode,
                    "embedding_dim": X.shape[1],
                    "pca_dim": PCA_DIM if pca is not None else None,
                    "n_classes": len(encoder.classes_),
                    "trained_on": "M1_patchseq_mouse_geneformer",
                },
                out_dir / f"{clf_name}__{mode}.joblib",
            )
            pd.DataFrame(
                {
                    "cell_id": common,
                    "true_label": encoder.inverse_transform(y),
                    "pred_label": encoder.inverse_transform(oof),
                }
            ).to_csv(out_dir / f"{clf_name}__{mode}_oof_predictions.csv", index=False)
            rows.append(
                {
                    "classifier": clf_name,
                    "mode": mode,
                    "n_cells": len(y),
                    "n_classes": len(encoder.classes_),
                    "accuracy": acc,
                    "balanced_accuracy": bal,
                    "macro_f1": macro,
                }
            )
    pd.DataFrame(rows).to_csv(out_dir / "summary.csv", index=False)


def maybe_log(y: np.ndarray, feature: str) -> tuple[np.ndarray, bool]:
    if feature in LOG_FEATURES and np.all(y >= 0):
        return np.log1p(y), True
    return y, False


def train_ephys(embs: pd.DataFrame, labels: pd.Series, ephys: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    common = embs.index.intersection(ephys.index)
    embs = embs.loc[common]
    ephys = ephys.loc[common]
    labels = labels.reindex(common)
    cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    rows = []
    for mode in ["embs", "pcs"]:
        for feature in ephys.columns:
            valid = ~ephys[feature].isna().values
            if valid.sum() < 50:
                continue
            log(f"[ephys] {feature}__{mode}")
            X = embs.values[valid]
            y_raw = ephys[feature].values[valid].astype(float)
            y, was_logged = maybe_log(y_raw, feature)
            strat = labels.iloc[np.where(valid)[0]].fillna("UNK").values
            oof_train_space = np.full(len(y), np.nan)
            fold_params = []
            for tr, te in cv.split(X, strat):
                _, _, X_tr, X_te = make_transform(mode, X[tr], X[te])
                model = ElasticNetCV(
                    l1_ratio=L1_RATIOS,
                    n_alphas=N_ALPHAS,
                    cv=INNER_CV,
                    max_iter=MAX_ITER,
                    n_jobs=-1,
                    random_state=SEED,
                    precompute=False,
                )
                model.fit(X_tr, y[tr])
                oof_train_space[te] = model.predict(X_te)
                fold_params.append({"alpha": float(model.alpha_), "l1_ratio": float(model.l1_ratio_)})
            pred = np.expm1(oof_train_space) if was_logged else oof_train_space
            r, p = pearsonr(y_raw, pred)
            mae = mean_absolute_error(y_raw, pred)
            r2 = r2_score(y_raw, pred)
            scaler, pca, X_full, _ = make_transform(mode, X)
            final = ElasticNetCV(
                l1_ratio=L1_RATIOS,
                n_alphas=N_ALPHAS,
                cv=INNER_CV,
                max_iter=MAX_ITER,
                n_jobs=-1,
                random_state=SEED,
                precompute=False,
            )
            final.fit(X_full, y)
            safe = sanitize(feature)
            dump(
                {
                    "model": final,
                    "scaler": scaler,
                    "pca": pca,
                    "was_logged": was_logged,
                    "feature": feature,
                    "mode": mode,
                    "pca_dim": PCA_DIM if pca is not None else None,
                    "embedding_dim": X.shape[1],
                    "trained_on": "M1_patchseq_mouse_geneformer",
                },
                out_dir / f"{safe}__{mode}.joblib",
            )
            pd.DataFrame({"y_true": y_raw, "y_pred": pred}).to_csv(out_dir / f"{safe}__{mode}_oof_predictions.csv", index=False)
            rows.append(
                {
                    "feature": feature,
                    "mode": mode,
                    "n_cells": int(valid.sum()),
                    "was_logged": was_logged,
                    "pearson_r": float(r),
                    "pearson_p": float(p),
                    "mae": float(mae),
                    "r2": float(r2),
                    "final_alpha": float(final.alpha_),
                    "final_l1_ratio": float(final.l1_ratio_),
                    "fold_params": fold_params,
                }
            )
    summary = pd.DataFrame([{k: v for k, v in row.items() if k != "fold_params"} for row in rows])
    summary.to_csv(out_dir / "summary.csv", index=False)
    summary.pivot(index="feature", columns="mode", values="pearson_r").to_csv(out_dir / "pearson_r_by_mode.csv")
    (out_dir / "fold_params.json").write_text(json.dumps(rows, indent=2))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", type=Path, default=DEFAULT_EMB)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--min-cells", type=int, default=8)
    parser.add_argument("--skip-celltype", action="store_true")
    parser.add_argument("--skip-ephys", action="store_true")
    args = parser.parse_args()
    seed()
    embs = load_embeddings(args.embeddings)
    labels = load_labels(args.min_cells)
    ephys = load_ephys()
    if not args.skip_celltype:
        train_celltype(embs, labels, args.out / "celltype_classifier")
    if not args.skip_ephys:
        train_ephys(embs, labels, ephys, args.out / "ephys_regressors")
    log(f"[done] outputs under {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
