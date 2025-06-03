# PREPS
The microenvironment of glioma is heterogeneous, including tumor cells, neurons, and immune cells, making it difficult to develop an effective treatment. Our previous study also demonstrated neuronal behaviors of glioma tumor cells, especially firing an action potential, highlighting the importance of characterizing the electrophysiological properties of each single cell. The electrophysiology data is achieved through Patch-sequencing experiments. However, the available data size is limited due to the experimental difficulty. Here, we introduce **PREPS** (**Pr**edicting **E**lectrophysiological **P**roperties of **S**ingle-cell RNA-seq), a machine-learning-based ***computational framework*** that employs the state-of-the-art **GPT** (**G**enerative **P**re-trained **T**ransformer) models to predict electrophysiological features of glioma samples by single-cell RNA-sequencing. 
  
## Methodology
In the first step of PREPS, leveraging the foundational GPT model, Geneformer, which has captured the complexity within human gene networks based on a broad range of healthy tissues, we **fine-tuned** the model into a series of brain-specific cell type classifiers using the transcriptomes of various developing brain and glioma datasets. Besides clustering and annotating glioma cells, we extracted and concatenated **embeddings** from the intermediate layers of these classifiers to represent the comprehensive transcriptomic features of each cell. Next, we built a group of predictive Elastic Nets (i.e., PREPS models) that **map** the electrophysiological features of glioma cells to their embeddings, with models optimized through a systematic grid search of all parameter combinations. Finally, we applied PREPS models to **predict** electrophysiological features of a larger amount of glioma data, where conducting many Patch-seq experiments is time-consuming and labor-intensive. 

### Fine-tuning
#### finetune.py
- This script fine-tunes the fundamental GPT model loaded from the directory `./Geneformer/` for a more specific context using a single reference dataset `./[name]/[name].dataset`.
- The fine-tuned model will be saved in the folder `./[name]/finetune`.

#### Usage
`$ python finetune.py [gpu] [name]`

#### Examples
`$ python finetune.py 0 aldinger_2000perCellType`
  
`$ python finetune.py 1 bhaduri_3000perCellType`

#### Notes
- Run `$ nvidia-smi` first to select an idle `[gpu]` with low Memory-Usage and GPU-Utility.
- The reference dataset should have been ***tokenized*** using `tokenize.py` and saved as `./[name]/[name].dataset`. See **Application - (2) Tokenization** for how `tokenize.py` works.
  
## Application
With the GPT models fine-tuned and the predictive PREPS models trained, it is easy to predict the electrophysiological features of a new scRNA-seq dataset (either human or mouse). Starting from an input `[seuratObj].rda`, the workflow consists of **(1) Data conversion**, **(2) Tokenization**, **(3) Annotation**, and **(4) Electrophysiological feature prediction**. Below, we demonstrate how PREPS works with a mouse scRNA-seq dataset.
  
### (1) Data conversion
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
- This script loads the scRNA-seq dataset `adata.h5ad` or {`meta.tsv`, `matrix.mtx`, `genes.tsv`, `barcodes.tsv`} from the directory `./[name]/`, converts them into an intermediate `[name].loom`, and tokenizes `[name].loom`, saving the results in a new folder `./[name]/[name].dataset/`.
- Human (default `[species]`) or mouse gene symbols will be mapped to human Ensembl IDs through the `GProfiler` online search.
  
#### Usage
`$ python tokenize.py [name] --species [species]`
  
#### Examples
`$ python tokenize.py mouse -s mouse`
  
`$ python tokenize.py glioma -s human`
  
#### Notes
- `adata.h5ad` or `matrix.mtx` should contain ***raw read counts***.
- Keep all genes and ***do not filter***.

### (3) Annotation
#### annotate.py
- This script loads the tokenized folder `[name].dataset` from the current directory, extracts cell embeddings, and annotates cell types using fine-tuned GPT models, saving the results in a new folder `./[name]_preds/`.
- Loading `[name].dataset` generates many temporary files within the folder. This script creates and works with `./[name]_preds/tokenized_copy.dataset` to keep `[name].dataset` clean for future use, similar to `finetune.py`.

#### Usage
`$ python annotate.py [gpu] [name]`

#### Examples
`$ python annotate.py 0 mouse`
  
`$ python annotate.py 1 glioma`

#### Notes
- Each fine-tuned GPT model's folder should have been saved in the ***current*** directory (e.g., `./aldinger_2000perCellType`, `./bhaduri_3000perCellType`).
- `./[name]_preds/tokenized_copy.dataset` can be deleted afterwards.

### (4) Electrophysiological feature prediction
#### patchseq_predict.py

#### Usage
`$ python patchseq_predict.py [name] --models [models]`

#### Examples
`$ python patchseq_predict.py mouse -m patchseq`
  
`$ python patchseq_predict.py glioma -m allen`

#### Notes
- 



## Installation

