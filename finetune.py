import argparse
parser = argparse.ArgumentParser(description='Fine-tune a pre-trained Geneformer model for a more specific context using a single tokenized dataset.')
parser.add_argument('ref_name', help='Input the reference dataset on which to fine-tune the model (e.g., aldinger_2000perCellType, bhaduri_3000perCellType).')
parser.add_argument('-g', '--gpu_name', choices=list(map(str, range(1000))), default='0', help='Input the idle GPU on which to run the code (e.g., 0, 1, 2).')
args = parser.parse_args()
ref_name = args.ref_name
gpu_name = args.gpu_name

import os
os.environ['CUDA_VISIBLE_DEVICES'] = gpu_name
os.environ['NCCL_DEBUG'] = 'INFO'

# imports
from collections import Counter
import datetime
import pickle
import subprocess
import seaborn as sns; sns.set()
from datasets import load_from_disk
from sklearn.metrics import accuracy_score, f1_score
from transformers import BertForSequenceClassification
from transformers import Trainer
from transformers.training_args import TrainingArguments

from geneformer import DataCollatorForCellClassification
import shutil
import pandas as pd


# load training dataset 
finetune_output_directory = f'{ref_name}/finetune/'
os.mkdir(finetune_output_directory)

# It is important to always create and work with a copy, as many temporary files will be created and mess the dataset
# When fine-tuning is done, the copy can be deleted while the original dataset keeps clean 
shutil.copytree(f'{ref_name}/{ref_name}.dataset', finetune_output_directory + 'tokenized_copy.dataset')
train_dataset = load_from_disk(finetune_output_directory + 'tokenized_copy.dataset')

dataset_list = []
evalset_list = []
organ_list = []
target_dict_list = []

for organ in Counter(train_dataset["isTumor"]).keys():
    assert organ in [0, 1]
    
    # Only those cells labeled with "isTumor = 0" are to be used for model finetuning 
    if organ == 1:
        continue
    
    organ_ids = [organ]
    organ_list += [organ]

    print(organ)

    num_proc = 1
    
    # filter datasets for given organ
    def if_organ(example):
        return example["isTumor"] in organ_ids
    trainset_organ = train_dataset.filter(if_organ, num_proc=num_proc)
    
    # per scDeepsort published method, drop cell types representing <0.5% of cells
    # celltype_counter = Counter(trainset_organ["group"])
    # total_cells = sum(celltype_counter.values())
    # cells_to_keep = [k for k,v in celltype_counter.items() if v>(0.005*total_cells)]
    # def if_not_rare_celltype(example):
    #     return example["group"] in cells_to_keep
    # trainset_organ_subset = trainset_organ.filter(if_not_rare_celltype, num_proc=num_proc)
    trainset_organ_subset = trainset_organ

    # shuffle datasets and rename columns
    trainset_organ_shuffled = trainset_organ_subset.shuffle(seed=42)
    trainset_organ_shuffled = trainset_organ_shuffled.rename_column("group", "label")
    trainset_organ_shuffled = trainset_organ_shuffled.remove_columns("isTumor")
    
    # change labels to one vs other
    # focus = ['interneuron']
    # def classes_to_focus(example):
    #     example["label"] = example["label"] if example["label"] in focus else 'other'
    #     return example
    # trainset_organ_shuffled = trainset_organ_shuffled.map(classes_to_focus, num_proc=num_proc)
    
    # create dictionary of cell types : label ids
    target_names = list(Counter(trainset_organ_shuffled["label"]).keys())
    target_name_id_dict = dict(zip(target_names,[i for i in range(len(target_names))]))
    target_dict_list += [target_name_id_dict]

    # save label id: cell type
    df_target_names = pd.DataFrame({'target_names': target_names})
    df_target_names.to_excel(f'{finetune_output_directory}target_names.xlsx', index=False, header=False)
    
    # change labels to numerical ids
    def classes_to_ids(example):
        example["label"] = target_name_id_dict[example["label"]]
        return example
    labeled_trainset = trainset_organ_shuffled.map(classes_to_ids, num_proc=num_proc)
    
    # create 80/20 train/eval splits
    labeled_train_split = labeled_trainset.select([i for i in range(0,round(len(labeled_trainset)*0.8))])
    labeled_eval_split = labeled_trainset.select([i for i in range(round(len(labeled_trainset)*0.8),len(labeled_trainset))])
    
    # filter dataset for cell types in corresponding training set
    trained_labels = list(Counter(labeled_train_split["label"]).keys())
    def if_trained_label(example):
        return example["label"] in trained_labels
    labeled_eval_split_subset = labeled_eval_split.filter(if_trained_label, num_proc=num_proc)

    dataset_list += [labeled_train_split]
    evalset_list += [labeled_eval_split_subset]

trainset_dict = dict(zip(organ_list,dataset_list))
traintargetdict_dict = dict(zip(organ_list,target_dict_list))

evalset_dict = dict(zip(organ_list,evalset_list))


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

# set training hyperparameters
# max learning rate
max_lr = 5e-5
# how many pretrained layers to freeze
freeze_layers = 0
# number gpus
num_gpus = 1
# number cpu cores
num_proc = 1
# batch size for training and eval
geneformer_batch_size = 12
# learning schedule
lr_schedule_fn = "linear"
# warmup steps
warmup_steps = 500
# number of epochs
epochs = 10
# optimizer
optimizer = "adamw"

for organ in organ_list:
    print(organ)
    organ_trainset = trainset_dict[organ]
    organ_evalset = evalset_dict[organ]
    organ_label_dict = traintargetdict_dict[organ]
    print('---')
    print(len(organ_label_dict.keys()))
    print('---')
    assert True
    
    # set logging steps
    logging_steps = round(len(organ_trainset)/geneformer_batch_size/10)
    
    # reload pretrained model
    model = BertForSequenceClassification.from_pretrained("Geneformer/", 
                                                      num_labels=len(organ_label_dict.keys()), 
                                                      ignore_mismatched_sizes=True, 
                                                      output_attentions = False,
                                                      output_hidden_states = False).to("cuda")
    
    # define output directory path
    current_date = datetime.datetime.now()
    datestamp = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}"
    output_dir = f"{finetune_output_directory}{datestamp}_geneformer_CellClassifier_{organ}_L{max_input_size}_B{geneformer_batch_size}_LR{max_lr}_LS{lr_schedule_fn}_WU{warmup_steps}_E{epochs}_O{optimizer}_F{freeze_layers}/"
    
    # ensure not overwriting previously saved model
    saved_model_test = os.path.join(output_dir, f"pytorch_model.bin")
    if os.path.isfile(saved_model_test) == True:
        raise Exception("Model already saved to this directory.")

    # make output directory
    subprocess.call(f'mkdir {output_dir}', shell=True)
    
    # set training arguments
    training_args = {
        "learning_rate": max_lr,
        "do_train": True,
        "do_eval": True,
        "evaluation_strategy": "epoch",
        "save_strategy": "epoch",
        "logging_steps": logging_steps,
        "group_by_length": True,
        "length_column_name": "length",
        "disable_tqdm": False,
        "lr_scheduler_type": lr_schedule_fn,
        "warmup_steps": warmup_steps,
        "weight_decay": 0.001,
        "per_device_train_batch_size": geneformer_batch_size,
        "per_device_eval_batch_size": geneformer_batch_size,
        "num_train_epochs": epochs,
        "load_best_model_at_end": True,
        "output_dir": output_dir,
    }
    
    training_args_init = TrainingArguments(**training_args)

    # create the trainer
    trainer = Trainer(
        model=model,
        args=training_args_init,
        data_collator=DataCollatorForCellClassification(),
        train_dataset=organ_trainset,
        eval_dataset=organ_evalset,
        compute_metrics=compute_metrics
    )
    # train the cell type classifier
    trainer.train()
    predictions = trainer.predict(organ_evalset)
    with open(f"{output_dir}predictions.pickle", "wb") as fp:
        pickle.dump(predictions, fp)
    trainer.save_metrics("eval",predictions.metrics)
    trainer.save_model(output_dir)

    print('---')
    print(ref_name)
    print(len(organ_label_dict.keys()))
    print('---')

# shutil.rmtree(f'{output_directory}tokenized_copy.dataset')
