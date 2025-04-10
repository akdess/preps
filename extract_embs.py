import os
import shutil
from geneformer import EmbExtractor
from tqdm import tqdm


test_name = 'spatial_raw_data_input'
output_directory = f'{test_name}_embs/'
os.mkdir(output_directory)
shutil.copytree(f'{test_name}.dataset', f'{output_directory}tokenized_copy.dataset')

ref_num_tups = [('aldinger_2000perCellType', 21), 
                ('allen_2000perCellType', 20), 
                ('bhaduri_3000perCellType', 10), 
                ('bhaduri_d2_4000perCellType', 10), 
                ('codex_1000perCellType', 16), 
                ('devbrain_3000perCellType', 10), 
                ('dirks_primary_gbm_combined_2000perCellType', 13), 
                ('primary_gbm_2000perCellType', 8), 
                ('recurrent_gbm_1000perCellType', 14), 
                ('TissueImmune_2000perCellType', 45)]

for ref_name, num_classes in tqdm(ref_num_tups):
    for emb_layer in [0, -1]:
        output_prefix = f'embs_by_{ref_name}_num_classes_{num_classes}_emb_layer_{emb_layer}'

        # initiate EmbExtractor
        embex = EmbExtractor(model_type='CellClassifier',
                            num_classes=num_classes,
                            emb_mode='cell', 
                            filter_data=None, 
                            max_ncells=1000*1000,
                            emb_layer=emb_layer,
                            emb_label=['individual', 'group'],
                            labels_to_plot=['group'],
                            forward_batch_size=100,
                            nproc=1)

        # extracts embedding from input data
        embs = embex.extract_embs(model_directory=f'{ref_name}/finetune/240605_geneformer_CellClassifier_0_L2048_B12_LR5e-05_LSlinear_WU500_E10_Oadamw_F0/', 
                                input_data_file=f'{output_directory}tokenized_copy.dataset',
                                output_directory=output_directory,
                                output_prefix=output_prefix)

        # plot UMAP of cell embeddings
        embex.plot_embs(embs=embs, 
                        plot_style='umap',
                        output_directory=output_directory,  
                        output_prefix=output_prefix, 
                        max_ncells_to_plot=1000*1000)
