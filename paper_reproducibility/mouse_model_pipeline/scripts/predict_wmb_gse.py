#!/usr/bin/env python3
"""Apply saved mouse PREPS cell-type classifiers and ephys regressors to WMB/GSE.

Examples are documented in ../README.md. The script expects bridge embedding files
with keys `embeddings` and `cell_ids` (the hidden+logits+probs concatenation used
by the saved models).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from joblib import load


PIPELINE_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path("/mnt/data/serinharmanci/preps_revision_052726")
MODELS_DIR = PIPELINE_DIR / "models"
EMBEDDINGS_DIR = PIPELINE_DIR / "embeddings"

DATASET_FILES = {
    "wmb": ["wmb_isocortex_1_embeddings.npz", "wmb_isocortex_2_embeddings.npz"],
    "gse": ["gse185862_ssv4_embeddings.npz"],
}

WMB_META = PROJECT_ROOT / "MOUSE_PATCHSEQ/Mouse_test_scell_data/metadata/WMB-10X_cell_metadata_with_cluster_annotation.tsv"
GSE_META = PROJECT_ROOT / "MOUSE_PATCHSEQ/Mouse_test_scell_data/metadata/GSE185862_metadata_ssv4.csv"

EPHYS_FEATURES = [
    "Input resistance (MOhm)",
    "Latency (ms)",
    "AP amplitude (mV)",
    "Rheobase (pA)",
    "Sag ratio",
    "Membrane time constant (ms)",
    "AP threshold (mV)",
    "Upstroke-to-downstroke ratio",
    "ISI adaptation index",
    "AP width (ms)",
]


def load_embeddings(path: Path) -> pd.DataFrame:
    arr = np.load(path, allow_pickle=True)
    df = pd.DataFrame(np.asarray(arr["embeddings"], dtype=np.float64), index=arr["cell_ids"])
    df.columns = [f"emb_{i:03d}" for i in range(df.shape[1])]
    df.index.name = "cell_id"
    return df


def iter_embedding_files(dataset: str, embedding_dir: Path) -> Iterable[Path]:
    for name in DATASET_FILES[dataset]:
        path = embedding_dir / name
        if not path.exists():
            raise FileNotFoundError(f"Missing embedding file: {path}")
        yield path


def apply_celltype(embs: pd.DataFrame, bundle: dict) -> pd.DataFrame:
    if bundle["embedding_dim"] != embs.shape[1]:
        raise ValueError(f"Cell-type model expects {bundle['embedding_dim']} dims, got {embs.shape[1]}")
    x_scaled = bundle["scaler"].transform(embs.values)
    if bundle.get("pca") is not None:
        x_scaled = bundle["pca"].transform(x_scaled)
    clf = bundle["classifier"]
    le = bundle["label_encoder"]
    pred_idx = clf.predict(x_scaled)
    pred = le.inverse_transform(pred_idx)
    out = pd.DataFrame({"y_pred": pred, "y_pred_top_prob": np.nan}, index=embs.index)
    if hasattr(clf, "predict_proba"):
        proba = clf.predict_proba(x_scaled)
        aligned = np.zeros((len(embs), len(le.classes_)), dtype=float)
        for local_i, class_id in enumerate(clf.classes_):
            aligned[:, int(class_id)] = proba[:, local_i]
        out["y_pred_top_prob"] = aligned.max(axis=1)
        for i, cls in enumerate(le.classes_):
            out[f"prob_{cls}"] = aligned[:, i]
    return out


def apply_ephys(embs: pd.DataFrame, models_dir: Path, mode: str) -> pd.DataFrame:
    cols: dict[str, np.ndarray] = {}
    for model_path in sorted(models_dir.glob(f"*__{mode}.joblib")):
        bundle = load(model_path)
        if bundle["embedding_dim"] != embs.shape[1]:
            raise ValueError(
                f"{model_path.name} expects {bundle['embedding_dim']} dims, got {embs.shape[1]}"
            )
        x_scaled = bundle["scaler"].transform(embs.values)
        if bundle.get("pca") is not None:
            x_scaled = bundle["pca"].transform(x_scaled)
        y = bundle["model"].predict(x_scaled)
        cols[bundle["feature"]] = np.expm1(y) if bundle["was_logged"] else y
    return pd.DataFrame(cols, index=embs.index)


def join_metadata(df: pd.DataFrame, dataset: str) -> pd.DataFrame:
    if dataset == "wmb":
        meta = pd.read_csv(WMB_META, sep="\t", low_memory=False)
        meta = meta[["cell_label", "class", "subclass", "supertype", "cluster"]].rename(
            columns={"cell_label": "cell_id"}
        )
    else:
        meta = pd.read_csv(GSE_META, low_memory=False)
        keep = ["sample_name", "class_label", "subclass_label", "cluster_label", "region_label", "platform_label"]
        meta = meta[[c for c in keep if c in meta.columns]].rename(columns={"sample_name": "cell_id"})
    meta = meta.drop_duplicates("cell_id")
    return df.reset_index().merge(meta, on="cell_id", how="left").set_index("cell_id")


def summarize_ephys(df: pd.DataFrame, dataset: str, out_dir: Path, mode: str) -> None:
    label_cols = ["class", "subclass"] if dataset == "wmb" else ["subclass_label", "class_label"]
    ephys_cols = [c for c in EPHYS_FEATURES if c in df.columns]
    for label_col in label_cols:
        if label_col not in df.columns:
            continue
        summary = (
            df.dropna(subset=[label_col])
            .groupby(label_col)[ephys_cols]
            .agg(["mean", "median", "count"])
        )
        summary.to_csv(out_dir / f"{dataset}_summary_by_{label_col}_{mode}.csv")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=["wmb", "gse", "both"], default="both")
    parser.add_argument("--classifier", default="logreg")
    parser.add_argument("--mode", default="embs", choices=["embs", "pcs"])
    parser.add_argument("--out", type=Path, default=PIPELINE_DIR / "example_outputs/predictions")
    args = parser.parse_args()

    celltype_model = MODELS_DIR / "celltype_classifier" / f"{args.classifier}__{args.mode}.joblib"
    ephys_models = MODELS_DIR / "ephys_regressors"
    embedding_dir = EMBEDDINGS_DIR
    if not celltype_model.exists():
        raise FileNotFoundError(celltype_model)
    celltype_bundle = load(celltype_model)

    datasets = ["wmb", "gse"] if args.dataset == "both" else [args.dataset]
    for dataset in datasets:
        out_dir = args.out / dataset
        out_dir.mkdir(parents=True, exist_ok=True)
        ct_pieces = []
        eph_pieces = []
        for emb_file in iter_embedding_files(dataset, embedding_dir):
            print(f"[{dataset}] {emb_file.name}")
            embs = load_embeddings(emb_file)
            ct = join_metadata(apply_celltype(embs, celltype_bundle), dataset)
            eph = join_metadata(apply_ephys(embs, ephys_models, args.mode), dataset)
            ct["source_file"] = emb_file.stem
            eph["source_file"] = emb_file.stem
            ct.to_csv(out_dir / f"{emb_file.stem}__{args.classifier}__{args.mode}_celltype_predictions.csv")
            eph.to_csv(out_dir / f"{emb_file.stem}_predicted_ephys_{args.mode}.csv")
            ct_pieces.append(ct)
            eph_pieces.append(eph)
        ct_all = pd.concat(ct_pieces)
        eph_all = pd.concat(eph_pieces)
        ct_all.to_csv(out_dir / f"{dataset}_combined__{args.classifier}__{args.mode}_celltype_predictions.csv")
        eph_all.to_csv(out_dir / f"{dataset}_combined_predicted_ephys_{args.mode}.csv")
        summarize_ephys(eph_all, dataset, out_dir, args.mode)
        print(f"[done] {dataset}: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
