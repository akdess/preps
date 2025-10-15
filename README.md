# PREPS 
The microenvironment of glioma is heterogeneous, including tumor cells, neurons, and immune cells, making it difficult to develop an effective treatment. Our previous study also demonstrated neuronal behaviors of glioma tumor cells, especially firing an action potential, highlighting the importance of characterizing the electrophysiological properties of each single cell. The electrophysiology data is achieved through Patch-sequencing experiments. However, the available data size is limited due to the experimental difficulty. Here, we introduce **PREPS** (**Pr**edicting **E**lectrophysiological **P**roperties of **S**ingle-cell RNA-seq), a machine-learning-based ***computational framework*** that employs the state-of-the-art **GPT** (**G**enerative **P**re-trained **T**ransformer) models to predict electrophysiological features of glioma samples by single-cell RNA-sequencing. 
  
## Methodology
In the first step of PREPS, leveraging the foundational GPT model, Geneformer, which has captured the complexity within human gene networks based on a broad range of healthy tissues, we **fine-tuned** the model into a series of brain-specific cell type classifiers using the transcriptomes of various developing brain and glioma datasets. Besides clustering and annotating glioma cells, we extracted and concatenated **embeddings** from the intermediate layers of these classifiers to represent the comprehensive transcriptomic features of each cell. Next, we built a group of predictive Elastic Nets (i.e., PREPS models) that **map** the electrophysiological features of glioma cells to their embeddings, with models optimized through a systematic grid search of all parameter combinations. Finally, we applied PREPS models to **predict** electrophysiological features of a larger amount of glioma data, where conducting many Patch-seq experiments is time-consuming and labor-intensive.
  
We also developed a single-cell gene set enrichment-like method (`scoring.py`) to assign cell types using gene **attention scores** derived from our fine-tuned transformer models. For each cell, we averaged multi-head attention weights from the final transformer layer and ranked genes based on the [CLS] token’s attention vector. Gene identifiers were converted to symbols, producing ranked gene lists per cell. To define marker sets, we automatically extracted and weighted marker genes for each cell type using PubMed abstracts (2021-2024) and GPT-4.1, prioritizing genes frequently cited or included in canonical brain cell markers. Using these weighted marker lists, we calculated enrichment scores per cell via a modified ssGSEA approach, assigning each cell to the highest scoring type. Final cell type labels were determined by consensus across multiple ranked gene inputs, and both enrichment scores and final annotations were exported.

### Fine-tuning
#### finetune.py
- This script fine-tunes the fundamental GPT model loaded from the directory `./Geneformer/` for a more specific context using a single reference dataset `./[ref_name]/[ref_name].dataset`.
- The fine-tuned model will be saved in the folder `./[ref_name]/finetune`.

#### Usage
`$ python finetune.py [ref_name] --gpu_name [gpu_name]`

#### Examples
`$ python finetune.py aldinger_2000perCellType`
  
`$ python finetune.py bhaduri_3000perCellType -g 2`

#### Notes
- The reference dataset should have been ***tokenized*** using `tokenize.py` and saved as `./[ref_name]/[ref_name].dataset`. See **Application - (2) Tokenization** for how `tokenize.py` works.
- Run `$ nvidia-smi` to select an idle `[gpu_name]` with low Memory-Usage and GPU-Utility, default `0`.
  
## Application
With the GPT models fine-tuned and the predictive PREPS models trained, it is easy to predict the electrophysiological features of a new scRNA-seq dataset (either human or mouse). Users can choose to run either the single script with the whole workflow integrated or separate scripts for flexible adjustment. Starting from an input `[seuratObj].rda` or `adata.h5ad`, the workflow consists of **(1) R Data conversion**, **(2) Tokenization**, **(3) Annotation**, and **(4) Electrophysiological feature/celltype prediction**. Below, we demonstrate how PREPS works with a mouse scRNA-seq dataset.

### Single-script whole workflow
#### preps.py
This script takes `adata.h5ad` from the directory `./[test_name]/` as input and generates predictions as output, with each separate step integrated into a whole workflow, including (2) Tokenization, (3) Annotation, and (4) Electrophysiological feature/celltype prediction. The tokenized data can be found in the directory `./[test_name]/[test_name].dataset/`. The GPT-based cell embeddings and reference cell type annotations are saved in the directory `./[test_name]_preds/`. The predicted electrophysiological features or cell types are saved in the directory `./[test_name]_[models]/`. 

#### Usage
`$ python preps.py [test_name] --species [species] --gpu_name [gpu_name] --models [models]`

#### Examples
`$ python preps.py mouse -s mouse -g 0 -m patchseq`
  
`$ python preps.py glioma -s human -g 1 -m celltype`

#### Notes
- `--species`: data source of `human` or `mouse`, default `human`.
- `--gpu_name`: on which GPU to run the code, `0`-`999`, default `0`.
- `--models`: predictive `patchseq` or `celltype` models to use, default `patchseq`.
- The parameter settings are the same as in separate scripts. More details are provided below.
  
### Separate scripts: (1) R Data conversion
If the scRNA-seq dataset `adata.h5ad` is available, skip this step and proceed to **(2) Tokenization**. Otherwise, suppose `[seuratObj].rda` is in the directory `./mouse/`. In `R`, we convert `seuratObj` into `meta.tsv`, `matrix.mtx`, `genes.tsv`, and `barcodes.tsv`, saving them in the same directory.
```
library(Matrix)
library(Seurat)

load("mouse/malcolm_rao_cx3cr1_mouse_032525_seuratObj_ann.rda")
write.table(seuratObj@meta.data, file = "mouse/meta.tsv", 
            sep = "\t", row.names = T, col.names = T, quote = F)
writeMM(seuratObj@assays$RNA@layers$counts, file = "mouse/matrix.mtx")
write.table(rownames(seuratObj), file = "mouse/genes.tsv", 
            sep = "\t", row.names = F, col.names = F, quote = F)
write.table(colnames(seuratObj), file = "mouse/barcodes.tsv", 
            sep = "\t", row.names = F, col.names = F, quote = F)
```
#### Notes
- In `meta.tsv`, the `colname` of cell IDs (i.e., barcodes) should be `CellID`.
- In `matrix.mtx`, ***raw read counts*** should be saved instead of processed or scaled data.

### (2) Tokenization
#### tokenize.py
- This script loads the scRNA-seq data `adata.h5ad` or the equivalent set {`meta.tsv`, `matrix.mtx`, `genes.tsv`, `barcodes.tsv`} from the directory `./[test_name]/`, converts them into an intermediate `[test_name].loom`, and tokenizes `[test_name].loom`, saving the results in a new folder `./[test_name]/[test_name].dataset/`.
- Human (the default `[species]`) or mouse gene symbols will be mapped to human Ensembl IDs through the `GProfiler` online search.
  
#### Usage
`$ python tokenize.py [test_name] --species [species]`
  
#### Examples
`$ python tokenize.py mouse -s mouse`
  
`$ python tokenize.py glioma -s human`
  
#### Notes
- `adata.h5ad` or `matrix.mtx` should contain ***raw read counts***.
- Keep all genes and ***do not filter***.

### (3) Annotation
#### annotate.py
- This script loads the tokenized data `[test_name].dataset` from the directory `./[test_name]/`, extracts cell embeddings, and annotates cell types using fine-tuned GPT models, saving the results in a new folder `./[test_name]_preds/`.
- Loading `[test_name].dataset` generates many temporary files within the folder. This script creates and works with `./[test_name]_preds/tokenized_copy.dataset` to keep `[test_name].dataset` clean for future use, similar to `finetune.py`.

#### Usage
`$ python annotate.py [test_name] --gpu_name [gpu_name]`

#### Examples
`$ python annotate.py mouse`
  
`$ python annotate.py glioma -g 1`

#### Notes
- Each fine-tuned GPT model's folder should have been saved in the ***current*** directory (e.g., `./aldinger_2000perCellType`, `./bhaduri_3000perCellType`).
- `./[test_name]_preds/tokenized_copy.dataset` can be deleted afterwards.

### (4) Electrophysiological feature/celltype prediction
#### patchseq_predict.py
- This script loads the pre-trained PREPS `[models]` (`patchseq` or `celltype`, default `patchseq`) to predict the electrophysiological features or cell types of the dataset `[test_name]` based on its cell embeddings loaded from the directory `./[test_name]_preds/`.
- The predicted features or cell types are saved in the directory `./[test_name]_[models]/`.

#### Usage
`$ python patchseq_predict.py [test_name] --models [models]`

#### Examples
`$ python patchseq_predict.py mouse -m patchseq`
  
`$ python patchseq_predict.py glioma -m celltype`

#### Notes
- The PREPS models with parameters grid-searched have been saved in the directory `./combined_patchseq_all_preds/`. ***Do not change*** the folder or file names that contain keys to identify the optimal model for each feature or cell type prediction.

### Attention-based DEEPS celltype scoring
#### scoring.py
- This script serves as the attention-based cell type scoring framework that combines transformer-derived gene attention values with literature-curated marker gene sets.
- It integrates multi-head attention patterns from fine-tuned Geneformer models and GPT-extracted marker weights from PubMed (queries of `[tumor]` OR `[tissue]`).
- Each cell in the testing dataset `[test_name]` is assigned an enrichment-based score across candidate cell types.

#### Usage
`$ python scoring.py [tumor] [test_name] --species [species] --tissue [tissue] --gpu_name [gpu_name]`

#### Examples
`$ python scoring.py DIPG mouse -s mouse -t brain -g 0`
  
`$ python scoring.py glioma glioma -s human -t brain -g 1`

#### Notes
- The testing dataset `[test_name]` should have been tokenized before running `scoring.py`.
