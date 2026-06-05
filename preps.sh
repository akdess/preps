#!/bin/bash
set -e # Stops execution immediately if any python script fails

# Command Line Arguments
WORKFLOW=$1
DATASET=$2
SPECIES=${3:-human} # Defaults to human if not specified

if [ -z "$WORKFLOW" ] || [ -z "$DATASET" ]; then
    echo "Usage: ./preps.sh [finetune|train|apply] [dataset_name] [human|mouse]"
    echo "Example: ./preps.sh apply glioma2 human"
    exit 1
fi

echo "========================================"
echo "🚀 Starting PREPS Workflow: $WORKFLOW"
echo "📂 Dataset: $DATASET | 🧬 Species: $SPECIES"
echo "========================================"

if [ "$WORKFLOW" == "finetune" ]; then
    python tokenize_data.py "$DATASET" -s "$SPECIES"
    python finetune.py "$DATASET"

elif [ "$WORKFLOW" == "train" ]; then
    python tokenize_data.py "$DATASET" -s "$SPECIES"
    python annotate.py "$DATASET"
    python patchseq_glm.py
    python select_best_models.py

elif [ "$WORKFLOW" == "apply" ]; then
    python tokenize_data.py "$DATASET" -s "$SPECIES"
    python annotate.py "$DATASET"
    python patchseq_predict.py "$DATASET" -m patchseq

else
    echo "❌ Invalid workflow. Choose finetune, train, or apply."
    exit 1
fi

echo "✅ Workflow $WORKFLOW completed successfully!"