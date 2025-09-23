import argparse
parser = argparse.ArgumentParser(description='Annotate one dataset using fine-tuned models. Output cell embeddings (_preds.csv) and cell-type scores (_scores.csv).')
parser.add_argument('test_name', help='Input the test dataset to be annotated (e.g., mouse, glioma).')
parser.add_argument('-g', '--gpu_name', choices=list(map(str, range(1000))), default='0', help='Input the idle GPU on which to run the code (e.g., 0, 1, 2).')
args = parser.parse_args()
test_name = args.test_name
gpu_name = args.gpu_name

import os
os.environ['CUDA_VISIBLE_DEVICES'] = gpu_name
os.environ['NCCL_DEBUG'] = 'INFO'

from datasets import load_from_disk
from sklearn.metrics import accuracy_score, f1_score
from transformers import BertForSequenceClassification
from transformers import Trainer
from transformers.training_args import TrainingArguments

from geneformer import DataCollatorForCellClassification
from scipy.special import softmax
import shutil
import numpy as np
import pandas as pd
from tqdm import tqdm


ann_output_directory = f'{test_name}_preds/'
os.mkdir(ann_output_directory)
shutil.copytree(f'{test_name}/{test_name}.dataset', ann_output_directory + 'tokenized_copy.dataset')

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
    output_prefix = f'preds_by_{ref_name}_num_classes_{num_classes}'
    model_directory=f'{ref_name}/finetune/240605_geneformer_CellClassifier_0_L2048_B12_LR5e-05_LSlinear_WU500_E10_Oadamw_F0/'


    # load data, labels, model
    tokenized_dataset = load_from_disk(ann_output_directory + 'tokenized_copy.dataset')
    # tokenized_dataset = tokenized_dataset.filter(lambda x: x['isTumor'] == 1, num_proc=1)
    labels = [0] * tokenized_dataset.num_rows
    tokenized_dataset = tokenized_dataset.add_column('label', labels)

    df_target_names = pd.read_excel(f'{ref_name}/finetune/target_names.xlsx', header=None)


    def compute_metrics(pred):
        labels = pred.label_ids
        preds = pred.predictions.argmax(-1)
        # calculate accuracy and macro f1 using sklearn's function
        acc = accuracy_score(labels, preds)
        macro_f1 = f1_score(labels, preds, average='macro')
        return {
        'accuracy': acc,
        'macro_f1': macro_f1
        }


    # set model parameters
    # max input size
    max_input_size = 2 ** 11  # 2048
    # number gpus
    num_gpus = 1
    # batch size for training and eval
    geneformer_batch_size = 12

    model = BertForSequenceClassification.from_pretrained(model_directory, 
                                                        num_labels=num_classes,
                                                        output_attentions=False,
                                                        output_hidden_states=False).to('cuda')

    # predict
    training_args = {
        "do_train": False,
        "do_eval": False,
        "evaluation_strategy": "epoch",
        "group_by_length": True,
        "length_column_name": "length",
        "disable_tqdm": False,
        "per_device_train_batch_size": geneformer_batch_size,
        "per_device_eval_batch_size": geneformer_batch_size,
        "output_dir": ann_output_directory
    }
    training_args_init = TrainingArguments(**training_args)
    trainer = Trainer(
        model=model,
        args=training_args_init,
        data_collator=DataCollatorForCellClassification(),
        train_dataset=tokenized_dataset,
        eval_dataset=tokenized_dataset
    )
    predictions = trainer.predict(tokenized_dataset)

    df_preds = pd.DataFrame(predictions.predictions, index=tokenized_dataset['individual'])
    df_preds.index.name = 'individual'
    df_preds.to_csv(ann_output_directory + f'{output_prefix}_preds.csv', sep=',')

    y_predict_prob = softmax(predictions.predictions, axis=1)
    y_predict_int = np.argmax(y_predict_prob, axis=1)
    y_predict_class = [df_target_names.iat[i, 0] for i in y_predict_int]
    y_predict_score = np.amax(y_predict_prob, axis=1)

    # save results
    col_ann = f'{ref_name}_ann'
    col_score = f'{ref_name}_score'

    df_predict_prob = pd.DataFrame(y_predict_prob, index=tokenized_dataset['individual'])
    df_predict_prob.index.name = 'individual'
    df_predict_prob.columns = df_target_names.iloc[:, 0].tolist()
    df_predict_prob[col_ann] = y_predict_class
    df_predict_prob[col_score] = y_predict_score
    df_predict_prob.to_csv(ann_output_directory + f'{output_prefix}_scores.csv', sep=',')
