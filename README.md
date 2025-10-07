# PREPS 
The microenvironment of glioma is heterogeneous, including tumor cells, neurons, and immune cells, making it difficult to develop an effective treatment. Our previous study also demonstrated neuronal behaviors of glioma tumor cells, especially firing an action potential, highlighting the importance of characterizing the electrophysiological properties of each single cell. The electrophysiology data is achieved through Patch-sequencing experiments. However, the available data size is limited due to the experimental difficulty. Here, we introduce **PREPS** (**Pr**edicting **E**lectrophysiological **P**roperties of **S**ingle-cell RNA-seq), a machine-learning-based ***computational framework*** that employs the state-of-the-art **GPT** (**G**enerative **P**re-trained **T**ransformer) models to predict electrophysiological features of glioma samples by single-cell RNA-sequencing. 
  
## Methodology
In the first step of PREPS, leveraging the foundational GPT model, Geneformer, which has captured the complexity within human gene networks based on a broad range of healthy tissues, we **fine-tuned** the model into a series of brain-specific cell type classifiers using the transcriptomes of various developing brain and glioma datasets. Besides clustering and annotating glioma cells, we extracted and concatenated **embeddings** from the intermediate layers of these classifiers to represent the comprehensive transcriptomic features of each cell. Next, we built a group of predictive Elastic Nets (i.e., PREPS models) that **map** the electrophysiological features of glioma cells to their embeddings, with models optimized through a systematic grid search of all parameter combinations. Finally, we applied PREPS models to **predict** electrophysiological features of a larger amount of glioma data, where conducting many Patch-seq experiments is time-consuming and labor-intensive. 

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
This script takes `adata.h5ad` from the directory `./[test_name]/` as input and generates annotations/predictions as output, with each separate step integrated into a whole workflow, including (2) Tokenization, (3) Annotation, and (4) Electrophysiological feature/celltype prediction.

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
- This script loads the tokenized folder `[test_name].dataset` from the directory `./[test_name]/`, extracts cell embeddings, and annotates cell types using fine-tuned GPT models, saving the results in a new folder `./[test_name]_preds/`.
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
