# PREPS
The microenvironment of glioma is heterogeneous, including tumor cells, neurons, and immune cells, making it difficult to develop an effective treatment. Our previous study also demonstrated neuronal behaviors of glioma tumor cells, especially firing an action potential, highlighting the importance of characterizing the electrophysiological properties of each single cell. The electrophysiology data is achieved through Patch-sequencing experiments. However, the available data size is limited due to the experimental difficulty. Here, we introduce **PREPS** (**Pr**edicting **E**lectrophysiological **P**roperties of **S**ingle-cell RNA-seq), a machine-learning-based ***computational framework*** that employs the state-of-the-art **GPT** (**G**enerative **P**re-trained **T**ransformer) models to predict electrophysiological features of glioma samples by single-cell RNA-sequencing. 
  
## Methodology
In the first step of PREPS, leveraging the foundational GPT model, Geneformer, which has captured the complexity within human gene networks based on a broad range of healthy tissues, we **fine-tuned** the model into a series of brain-specific cell type classifiers using the transcriptomes of various developing brain and glioma datasets. Besides clustering and annotating glioma cells, we extracted and concatenated **embeddings** from the intermediate layers of these classifiers to represent the comprehensive transcriptomic features of each cell. Next, we built a group of predictive Elastic Nets (i.e., PREPS models) that **map** the electrophysiological features of glioma cells to their embeddings, with models optimized through a systematic grid search of all parameter combinations. Finally, we applied PREPS models to **predict** electrophysiological features of a larger amount of glioma data, where conducting many Patch-seq experiments is time-consuming and labor-intensive. 
  
## Application
The workflow of PREPS starts from the tokenization of the gene expression matrices of all datasets. 
  
### Tokenization
PREPS loads [file_name].h5ad file

### Fine-tuning


### Annotation


### PREPS model training


### Electrophysiological feature prediction

## Installation

