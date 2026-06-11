import argparse
parser = argparse.ArgumentParser(description='Attention-Based DEEPS Cell Type Scoring.')
parser.add_argument('tumor', help='Input the tumor for retrieval of cell markers (e.g., DIPG, glioma).')
parser.add_argument('test_name', help='Input the directory name of the tokenized dataset (e.g., mouse, glioma).')
parser.add_argument('-s', '--species', choices=['human', 'mouse'], default='human', help='Input -s human or -s mouse to designate species (default human).')
parser.add_argument('-t', '--tissue', default='brain', help='Input the tissue for retrieval of cell markers (default brain).')
parser.add_argument('-g', '--gpu_name', choices=list(map(str, range(1000))), default='0', help='Input the idle GPU on which to run the code (e.g., 0, 1, 2).')

args = parser.parse_args()
tumor = args.tumor
test_name = args.test_name
species = args.species
tissue = args.tissue
gpu_name = args.gpu_name

import os
os.environ['CUDA_VISIBLE_DEVICES'] = gpu_name
os.environ['NCCL_DEBUG'] = 'INFO'

import random
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import spearmanr, pearsonr
import torch
import re
import csv
import time
import requests
import xml.etree.ElementTree as ET
import openai
from collections import defaultdict, Counter
import glob
import shutil
import pickle
import mygene
from joblib import Parallel, delayed
from tqdm import tqdm
from datasets import load_from_disk
from geneformer import Classifier, DataCollatorForCellClassification
from transformers import BertForSequenceClassification


os.environ['PYTHONHASHSEED'] = '0'
random.seed(0)
np.random.seed(0)
torch.manual_seed(0)
torch.set_grad_enabled(False)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


# Set your OpenAI API key
client = openai.OpenAI(api_key="")


def extract_markers_with_chatgpt(title, abstract, cell_type_name):
    prompt = (
        f"Given the following abstract about {tumor}:\n\n"
        f"Title: {title}\n"
        f"Abstract: {abstract}\n\n"
        f"Extract gene markers for cell type {cell_type_name} mentioned in the context of {tumor} or {tissue}. "
        f"Return results in a table with three columns: Cell Type and {species} Gene Symbol Marker Genes and Title"
    )
    try:
        response = client.chat.completions.create(
            model="gpt-4.1",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.3
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Error: {e}"


def fill_remake(markers):
    remake = []
    lines = markers.strip().split("\n")
    for line in lines:
        if line.startswith("|") and len(line.strip()) > 1 and'Cell Type'not in line and "----" not in line:
            parts = [part.strip() for part in line.split("|") if part.strip()]
            if len(parts) == 3:
                cell_type, marker, title = parts
                marker_genes = [m.strip() for m in marker.split(",") if m.strip()]
                for gene in marker_genes:
                    remake.append((cell_type, gene, title))
    return remake


def fetch_pubmed_abstracts(cell_type_name, max_results=500, output_file=f"{tumor}_cell_markers.txt"):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    query = f"{cell_type_name} AND marker AND cell type AND ({tumor} OR {tissue}) AND 2021:2024[PDAT]"
    esearch_url = f"{base_url}esearch.fcgi?db=pubmed&term={query}&retmax={max_results}&retmode=xml" # gives 500 results per cell 
    esearch_response = requests.get(esearch_url)
    esearch_data = ET.fromstring(esearch_response.text)
    pmids = [id_tag.text for id_tag in esearch_data.findall(".//Id")]
    if os.path.exists(output_file):
        mode = 'a'
    else:
        mode = 'w'

    with open(output_file,mode, encoding="utf-8", newline="") as file:
        csv_writer = csv.writer(file)
        if mode == 'w':
            csv_writer.writerow(["CellType", "MarkerGene", "Title"]) # write header once, file is fresh

        for start in range(0, len(pmids), 10):
            batch_pmids = pmids[start:start+10]
            efetch_url = f"{base_url}efetch.fcgi?db=pubmed&id={','.join(batch_pmids)}&retmode=xml&rettype=abstract"
            efetch_response = requests.get(efetch_url)
            try:
                efetch_data = ET.fromstring(efetch_response.text)
            except ET.ParseError as e:
                print(f"XML ParseError for PMIDs {batch_pmids}: {e}")
                with open("malformed_pubmed_response.xml", "w", encoding="utf-8") as f:
                    f.write(efetch_response.text)
                continue
            for article in efetch_data.findall(".//PubmedArticle"):
                title = article.findtext(".//ArticleTitle", default="No Title").strip()
                abstract = article.findtext(".//AbstractText", default="No Abstract").strip()
                pmid = article.findtext(".//PMID", default="No PMID")
                markers = extract_markers_with_chatgpt(title, abstract, cell_type_name)
                if markers.startswith("Error"):
                    print(markers)
                    continue

                tuples = fill_remake(markers)
                for cell_type, gene, title in tuples:
                    csv_writer.writerow([cell_type, gene, title])

            print(f"Processed abstracts {start + 1} to {start + len(batch_pmids)}")
            time.sleep(5) # avoid rate limits

    print(f"Results saved to {output_file}")


def cell_fetch_name(file_path=f"{tumor}_cell_types.txt"):
    if os.path.exists(file_path):
        with open(file_path, "r") as file:
            cell_types = [line.strip() for line in file if line.strip()]
            
            for cell_type_name in cell_types: # runs each cell name through the function
                fetch_pubmed_abstracts(cell_type_name)


# Run it
cell_fetch_name()


def get_cell_type_markers(cell_type, species="human", tissue="brain"):
    prompt = (
        f"List well-established gene markers for {cell_type} "
        f"in the {tissue} of {species}. Format the response as a table "
        f"with '{species} Gene Symbol' and 'Description'."
    )

    response = client.chat.completions.create(
        model="gpt-4.1",
        messages=[{"role": "user", "content": prompt}],
        temperature=0.2,
    )
    return response.choices[0].message.content


file_path = f"{tumor}_cell_types.txt"
cleaned_lines = []
cellmarker_file = "cellmarker.csv"

with open(cellmarker_file, 'w', newline='', encoding='utf-8') as new_file:
    csv_writer = csv.writer(new_file)
    csv_writer.writerow(["celltype", "geneName"])

    if os.path.exists(file_path):
        with open(file_path, "r") as file:
            cell_types = [line.strip() for line in file if line.strip()]

    for cell_type_input in cell_types:
        output = get_cell_type_markers(cell_type_input, species=species, tissue=tissue)

        lines = output.strip().split("\n")
        for line in lines:
            if line.startswith("|") and len(line.strip()) > 1 and'Gene Symbol'not in line and "----" not in line:
                parts = [part.strip() for part in line.split("|") if part.strip()]
                if len(parts) == 2:
                    gene, desc = parts
                    cleaned_lines.append((gene, desc)) 
                    print(f"{gene}\n")
                    csv_writer.writerow([cell_type_input, gene]) 
                    '''
                    jont +=1
                    if jont >= 2:
                        add = len(cell_type_input)*' '
                        csv_writer.writerow([f"|{add}| {gene}|"])
                        if jont == 10:
                                add = len(cell_type_input + gene)*'-'
                    else:
                            csv_writer.writerow([f"|{cell_type_input}| {gene}|"])''' # some code to make it more readable not needed 
                    # do not add not in right form


input_file = f"{tumor}_cell_markers.txt"
cellmarker_file = "cellmarker.csv"
combined_output_file = "combined_marker_reference.csv"

marker_refs = defaultdict(set)
with open(input_file, newline='', encoding='utf-8') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        cell = re.sub(r'["“”]', '', row["CellType"]).strip().lower().title()  
        gene = re.split(r'\s*\(', row["MarkerGene"])[0].strip().rstrip(')')
        gene = re.sub(r'[^A-Za-z0-9\-]+', '', gene) 
        if gene.lower().startswith("none") or gene.lower().startswith("but") or gene.lower().startswith("not") or gene.lower() in {""} or len(gene) > 15:
            continue  # skip non-gene entries
        title = row["Title"].strip()
        marker_refs[(cell, gene)].add(title)

cellmarker_set = set()
with open(cellmarker_file, newline='', encoding='utf-8') as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        cell = re.sub(r'["“”]', '', row["celltype"]).strip().lower() 
        gene = re.split(r'\s*\(', row["geneName"])[0].strip().rstrip(')')
        gene = re.sub(r'[^A-Za-z0-9\-]+', '', gene)  
        if gene.lower().startswith("none") or gene.lower().startswith("but") or  gene.lower().startswith("not") or gene.lower() in {""} or len(gene) > 15:
            continue  # skip non-gene entries
        cellmarker_set.add((cell, gene))

with open(combined_output_file, "w", newline='', encoding='utf-8') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(["CellType", "MarkerGene", "NumTitles", "InCellMarker"])

    for (cell, gene), titles in sorted(marker_refs.items()):
        if not cell or cell[0] in "*(-":
            continue
        key = (cell.lower(), gene) 
        in_cellmarker = "Yes" if key in cellmarker_set else "No"
        num_titles = len(titles)
        if num_titles > 5 or in_cellmarker == "Yes":
            writer.writerow([cell, gene, num_titles, in_cellmarker])


# ----------------------------
# Paths and data
tokenized_input_path = f'{test_name}.dataset'
output_directory = f'{test_name}_attention_ranked_genes/'
os.makedirs(output_directory, exist_ok=True)

# Make a working copy of the tokenized dataset
dst_tokenized = f'{output_directory}tokenized_copy.dataset'
if os.path.exists(dst_tokenized):
    shutil.rmtree(dst_tokenized)
shutil.copytree(tokenized_input_path, dst_tokenized)

# Reference models
ref_num_tups = [
    ('aldinger_2000perCellType', 21), ('allen_2000perCellType', 20),
    ('bhaduri_3000perCellType', 10), ('bhaduri_d2_4000perCellType', 10),
    ('codex_1000perCellType', 16), ('devbrain_3000perCellType', 10),
    ('dirks_primary_gbm_combined_2000perCellType', 13),
    ('primary_gbm_2000perCellType', 8),
    ('recurrent_gbm_1000perCellType', 14),
    ('TissueImmune_2000perCellType', 45)
]

# ----------------------------
# Load token dictionary (symbol -> token_id); build reverse map (token_id -> symbol)
with open("/mnt/data/serinharmanci/preps_gene_emb/token_dictionary_gc30M.pkl", "rb") as f:
    token_dict = pickle.load(f)
id_to_gene = {v: k for k, v in token_dict.items()}  # token_id -> gene_symbol

# Optional: Ensembl -> gene symbol (only used if your tokens are Ensembl IDs)
mg = mygene.MyGeneInfo()
def map_ensembl_to_symbol(ensembl_ids):
    if len(ensembl_ids) == 0:
        return {}
    results, batch_size = [], 1000
    for i in range(0, len(ensembl_ids), batch_size):
        batch = ensembl_ids[i:i+batch_size]
        results.extend(mg.querymany(batch, scopes='ensembl.gene', fields='symbol', species=species))
    return {item['query']: item.get('symbol', item['query']) for item in results}

# ----------------------------
# Load tokenized HF dataset once
dataset = load_from_disk(dst_tokenized)
# Add dummy labels if not present
if 'label' not in dataset.column_names:
    dataset = dataset.add_column("label", [0] * len(dataset))
# Expect columns: 'input_ids', 'attention_mask', maybe 'individual', 'length', ...
has_individual = 'individual' in dataset.column_names

# Collator for padding to batch (consistent with your pipeline)
collator = DataCollatorForCellClassification(token_dictionary=token_dict)
# Batching params
batch_size = 12

def batched(iterable, n):
    for i in range(0, len(iterable), n):
        yield range(i, min(i + n, len(iterable)))

# ----------------------------
# Main Loop Over Fine-Tuned Models
for ref_name, num_classes in tqdm(ref_num_tups, desc="References"):
    model_path = f'{ref_name}/finetune/240605_geneformer_CellClassifier_0_L2048_B12_LR5e-05_LSlinear_WU500_E10_Oadamw_F0/'
    # Load classifier; force attentions on
    model = BertForSequenceClassification.from_pretrained(
    model_path,
    num_labels=num_classes,
    output_attentions=True,
    output_hidden_states=False).to(device)
    model.eval()
    
    # ensure we get attentions from forward
    if hasattr(model, "config"):
        model.config.output_attentions = True
    
    all_cells = []
    # Process in batches
    for idx_range in tqdm(batched(range(len(dataset)), batch_size), total=(len(dataset)+batch_size-1)//batch_size, desc=f"{ref_name}"):
        # Build a batch list of dicts expected by collator
        # Retain original cell IDs
        raw_examples = [dataset[i] for i in idx_range]
        if has_individual:
            indivs = [ex['individual'] for ex in raw_examples]
        else:
            indivs = [f"{ref_name}_cell_{i}" for i in idx_range]

        # Keep only model-relevant keys
        examples = [{k: v for k, v in ex.items() if k in ["input_ids", "attention_mask", "label"]} for ex in raw_examples]
        batch = collator(examples)

        input_ids = torch.tensor(batch["input_ids"], dtype=torch.long, device=device)
        attention_mask = torch.tensor(batch["attention_mask"], dtype=torch.long, device=device)

        with torch.no_grad():
            outputs = model(input_ids=input_ids, attention_mask=attention_mask, output_attentions=True)

        # outputs.attentions: tuple(len = n_layers) of tensors [B, n_heads, L, L]
        attn_last = outputs.attentions[-1]                     # [B, H, L, L]
        attn_avg  = attn_last.mean(dim=1)                      # [B, L, L] (avg heads)

        # For each item in batch, extract [CELL](pos 0) -> gene tokens (1..len_i-1)
        # Use attention_mask to get true sequence length per sample
        seq_lens = attention_mask.sum(dim=1).tolist()          # list of ints

        # We also need raw input_ids to map gene tokens
        # NOTE: collator may pad; we must slice valid region per sample
        input_ids_np = input_ids.detach().cpu().numpy()

        # Carry optional identifiers
        if has_individual:
            indivs = [dataset[i]['individual'] for i in idx_range]
        else:
            indivs = [f"{ref_name}_cell_{i}" for i in idx_range]

        for bi, (cell_id, L) in enumerate(zip(indivs, seq_lens)):   
            if L <= 1:
                continue           
            attn_mat = attn_avg[bi, :L, :L]
            cell_to_gene = attn_mat[0, 1:]
            gene_token_ids = input_ids_np[bi, 1:L].tolist()
            gene_syms = [id_to_gene.get(tid, f"UNK_{tid}") for tid in gene_token_ids]
            sorted_idx = torch.argsort(cell_to_gene, descending=True).detach().cpu().numpy()
            ranked_genes = [gene_syms[i] for i in sorted_idx]
            ranked_scores = cell_to_gene[sorted_idx].detach().cpu().numpy()

            df_cell = pd.DataFrame({
                "cell_id": cell_id,
                "gene": ranked_genes,
                "attention": ranked_scores
            })
            all_cells.append(df_cell)
           
    if all_cells:
        df_all = pd.concat(all_cells, ignore_index=True)
        out_path = os.path.join(output_directory, f"{ref_name}_attn_rank_ALLCELLS.csv")
        df_all.to_csv(out_path, index=False)
        print(f"Saved: {out_path}")


mg = mygene.MyGeneInfo()

# ----------------------------
# Inputs
embedding_folder = f"/mnt/data/serinharmanci/preps_gene_emb/{test_name}_attention_ranked_genes"
pattern = os.path.join(embedding_folder, "*_attn_rank_ALLCELLS.csv")
#pattern = os.path.join(embedding_folder, "recurrent_gbm_1000perCellType_attn_rank_ALLCELLS.csv")
rank_files = glob.glob(pattern)

marker_file = "/mnt/data/serinharmanci/preps_gene_emb/combined_marker_reference.csv"
output_dir = "/mnt/data/serinharmanci/preps_gene_emb/output_predictions_attention"
os.makedirs(output_dir, exist_ok=True)

# ----------------------------
# Load marker gene sets
marker_df = pd.read_csv(marker_file)
assert {"CellType", "MarkerGene", "NumTitles"}.issubset(marker_df.columns)

marker_dict = {
    cell_type: dict(zip(group["MarkerGene"], group["NumTitles"]))
    for cell_type, group in marker_df.groupby("CellType")
}
universe_gene_set = set(marker_df["MarkerGene"].unique())

# ----------------------------
# Gene symbol conversion
def ensembl_to_symbol(ensembl_ids):
    result = mg.querymany(
        ensembl_ids,
        scopes="ensembl.gene",
        fields="symbol",
        species=species,
        as_dataframe=True
    )
    return result["symbol"].to_dict()

# ----------------------------
# ssGSEA scoring function
def ssgsea_score(gene_ranks: dict, marker_weight_dict: dict, p: float = 0.0):
    gene_list = sorted(gene_ranks.items(), key=lambda x: x[1])
    ranked_genes = [g for g, _ in gene_list]
    N = len(ranked_genes)

    marker_genes = [g for g in marker_weight_dict if g in gene_ranks]
    if len(marker_genes) < 2:
        return 0.0

    Nh = len(marker_genes)
    raw_weights = np.array([marker_weight_dict[g] for g in marker_genes], dtype=np.float64)
    normalized_weights = {g: w / raw_weights.sum() for g, w in zip(marker_genes, raw_weights)}

    running_score = 0
    score_profile = []

    for i, gene in enumerate(ranked_genes):
        hit = np.power(i + 1, p) * normalized_weights.get(gene, 0.0)
        miss = 0.0 if gene in normalized_weights else 1.0 / (N - Nh)
        running_score += hit - miss
        score_profile.append(running_score)

    return max(score_profile)

# ----------------------------
# Function to process a single cell
def process_cell(cell_id, group, marker_dict, full_gene_set, rank_file_id):
    group = group.sort_values("attention", ascending=False)
    ranked_genes = group["gene_symbol"].tolist()
    default_rank = len(ranked_genes)

    gene_ranks = defaultdict(lambda: default_rank)
    for i, gene in enumerate(ranked_genes):
        gene_ranks[gene] = i

    cell_scores = {"cell_id": cell_id, "rank_file_id": rank_file_id}
    score_values = {}

    for celltype, gene_weights in marker_dict.items():
        score = ssgsea_score(gene_ranks, gene_weights)
        cell_scores[celltype] = score
        score_values[celltype] = score

    if score_values:
        top_celltype = max(score_values, key=score_values.get)
        cell_scores["max_score"] = score_values[top_celltype]
    else:
        cell_scores["max_score"] = 0.0

    return cell_id, cell_scores, cell_scores.get("max_score"), top_celltype if score_values else None

# ----------------------------
# Score all cells in parallel
cell_top_predictions = defaultdict(list)
all_scores = []

for file in rank_files:
    rank_file_id = os.path.basename(file).replace("_attn_rank_ALLCELLS.csv", "")
    df = pd.read_csv(file)

    ensembl_to_symbol_map = ensembl_to_symbol(df["gene"].unique().tolist())
    df["gene_symbol"] = df["gene"].map(ensembl_to_symbol_map)
    df = df.dropna(subset=["gene_symbol"])

    assert {"cell_id", "gene_symbol", "attention"}.issubset(df.columns)
    observed_genes = set(df["gene_symbol"].unique())
    full_gene_set = universe_gene_set | observed_genes

    results = Parallel(n_jobs=32, backend="loky")(
        delayed(process_cell)(cell_id, group, marker_dict, full_gene_set, rank_file_id)
        for cell_id, group in tqdm(df.groupby("cell_id"), desc=rank_file_id)
    )

    for cell_id, cell_scores, max_score, top_celltype in results:
        all_scores.append(cell_scores)
        if top_celltype:
            cell_top_predictions[cell_id].append(top_celltype)

# ----------------------------
# Build output DataFrames
df_scores = pd.DataFrame(all_scores)

score_columns = df_scores.drop(columns=["cell_id", "rank_file_id", "max_score"], errors="ignore").select_dtypes(include="number").columns

max_scores_df = df_scores[["cell_id", "rank_file_id", "max_score"]].copy()
max_scores_df["top_cell_type"] = df_scores[score_columns].idxmax(axis=1)

max_scores_df.to_csv(
    os.path.join(output_dir, "attention_per_file_max_scores_and_types.csv"),
    index=False
)
print("Saved per-rank-file max scores and top cell types.")

# ----------------------------
# Majority vote across files
majority_vote = {
    cell_id: Counter(labels).most_common(1)[0][0]
    for cell_id, labels in cell_top_predictions.items()
}
df_scores["voted_cell_type"] = df_scores["cell_id"].map(majority_vote)

df_scores.to_csv(
    os.path.join(output_dir, "attention_weighted_ssgsea_scores.csv"),
    index=False
)
print("Saved full scores with majority-voted labels.")

majority_vote_df = pd.DataFrame([
    {"cell_id": cid, "voted_cell_type": ctype}
    for cid, ctype in majority_vote.items()
])
majority_vote_df.to_csv(
    os.path.join(output_dir, "attention_majority_voted_celltypes.csv"),
    index=False
)
print("Saved majority vote assignments to attention_majority_voted_celltypes.csv")


# ----------------------------
# Feature plot
scores_csv = os.path.join(output_dir, "attention_weighted_ssgsea_scores.csv")
adata_path = "adata.h5ad"
reduction_name = "X_umap_geneformer2"
out_root = "FeaturePlots_by_rank"
os.makedirs(out_root, exist_ok=True)

# Load data
scores = pd.read_csv(scores_csv)
adata = sc.read_h5ad(adata_path)

# Keep only cells present in AnnData
scores = scores[scores["cell_id"].isin(adata.obs_names)]

# Identify score columns
exclude_cols = ["cell_id", "rank_file_id", "max_score", "voted_cell_type"]
score_cols = [c for c in scores.columns if c not in exclude_cols]

# Average across rank files
mean_df = (
    scores.groupby("cell_id")[score_cols]
    .mean()
    .reset_index()
)

# Replace NaN with None
mean_df = mean_df.replace({np.nan: None})

# Align order with AnnData
mean_df = mean_df.set_index("cell_id").loc[adata.obs_names]
assert all(mean_df.index == adata.obs_names)

# Assign cell type based on threshold
# Assign cell type with highest scoring > threshold
threshold = 0.5
def assign_cell_type(row):
    top_type = row[score_cols].idxmax()
    top_score = row[top_type]
    return top_type if top_score > threshold else "Unassigned"

# Assign all cell types > threshold
threshold = 0.3
def assign_cell_type(row):
    row_scores = row[score_cols]
    if row_scores.max() < threshold:
        return "Unassigned"
    top_types = row_scores[row_scores >= threshold].sort_values(ascending=False).index.tolist()
    return ";".join(top_types)

mean_df["cell_type"] = mean_df.apply(assign_cell_type, axis=1)
mean_df.to_csv(os.path.join(output_dir, "attention_weighted_ssgsea_scores_with_celltype.csv"), index=False)

# Add to metadata
for col in score_cols:
    adata.obs[f"avg__{col}"] = mean_df[col]

# Plot each averaged feature
for col in score_cols:
    feature = f"avg__{col}"
    if not adata.obs[feature].isna().all():
        sc.pl.umap(
            adata,
            color=feature,
            color_map="viridis",
            show=False,
            title=f"{col} (mean across ranks)",
            save=f"_{col}_FeaturePlot_mean.png"
        )


# ----------------------------
# DEEPS vs Seurat AddModuleScore
deeps_path = os.path.join(output_dir, "attention_weighted_ssgsea_scores.csv")
seurat_path = "seurat_metadata_with_module_scores.csv"
umap_path = "umap_embeddings.csv"

deeps_df = pd.read_csv(deeps_path)
seurat_df = pd.read_csv(seurat_path)
umap_df = pd.read_csv(umap_path)

assert "cell_id" in deeps_df.columns and "cell_id" in seurat_df.columns
merged = pd.merge(deeps_df, seurat_df, on="cell_id", suffixes=("_deeps", "_seurat"))
merged = merged.merge(umap_df, on="cell_id", how="left")

# Normalize scores
deeps_cols = [c for c in merged.columns if c.endswith("_deeps")]
seurat_cols = [c for c in merged.columns if c.endswith("_seurat")]

scaler = MinMaxScaler()
merged[deeps_cols] = scaler.fit_transform(merged[deeps_cols])
merged[seurat_cols] = scaler.fit_transform(merged[seurat_cols])

# Correlation
results = []
for cell_type in [c.replace("_deeps", "") for c in deeps_cols]:
    d_col = f"{cell_type}_deeps"
    s_col = f"{cell_type}_seurat"
    if d_col in merged.columns and s_col in merged.columns:
        valid = merged[[d_col, s_col]].dropna()
        if len(valid) > 10:
            pearson_r, _ = pearsonr(valid[d_col], valid[s_col])
            spearman_r, _ = spearmanr(valid[d_col], valid[s_col])
            results.append({
                "CellType": cell_type,
                "Pearson_r": pearson_r,
                "Spearman_r": spearman_r,
                "N_cells": len(valid)
            })

corr_df = pd.DataFrame(results)
corr_df.to_csv("comparison_correlation_summary.csv", index=False)

# Scatter plot
sns.set(style="whitegrid")
for cell_type in [c.replace("_deeps", "") for c in deeps_cols]:
    d_col, s_col = f"{cell_type}_deeps", f"{cell_type}_seurat"
    if d_col in merged.columns and s_col in merged.columns:
        plt.figure(figsize=(5,5))
        sns.scatterplot(data=merged, x=s_col, y=d_col, s=10, alpha=0.5)
        sns.regplot(data=merged, x=s_col, y=d_col, scatter=False, color="red", ci=None)
        plt.title(f"{cell_type}: DEEPS vs Seurat (Spearman = {corr_df.loc[corr_df.CellType==cell_type, 'Spearman_r'].values[0]:.2f})")
        plt.xlabel("Seurat AddModuleScore")
        plt.ylabel("DEEPS Attention Score")
        plt.tight_layout()
        plt.savefig(f"scatter_DEEPS_vs_Seurat_{cell_type}.png", dpi=300)
        plt.close()

# UMAP
for cell_type in [c.replace("_deeps", "") for c in deeps_cols]:
    for method, col_suffix in [("DEEPS", "_deeps"), ("Seurat", "_seurat")]:
        score_col = f"{cell_type}{col_suffix}"
        if score_col in merged.columns:
            plt.figure(figsize=(5,5))
            plt.scatter(
                merged["UMAP_1"], merged["UMAP_2"],
                c=merged[score_col], cmap="viridis", s=6
            )
            plt.title(f"{method} {cell_type} score")
            plt.axis("off")
            plt.tight_layout()
            plt.savefig(f"UMAP_{method}_{cell_type}.png", dpi=300)
            plt.close()

# barplot of correlation
plt.figure(figsize=(8, 4))
sns.barplot(data=corr_df.melt(id_vars="CellType", value_vars=["Pearson_r","Spearman_r"]),
            x="CellType", y="value", hue="variable")
plt.xticks(rotation=45, ha="right")
plt.ylabel("Correlation Coefficient")
plt.title("DEEPS vs Seurat Correlation per Cell Type")
plt.tight_layout()
plt.savefig("DEEPS_vs_Seurat_CorrelationSummary.png", dpi=300)
plt.close()
