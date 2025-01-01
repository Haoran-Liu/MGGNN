# MGGNN
Marker Gene-Guided Graph Neural Networks for Enhanced Spatial Transcriptomics Clustering

# Setting Up

To get started, create a conda environment for MGGNN by following the instructions in conda/conda.sh.

# Spot Identification

The spots_identification folder contains code for identifying spots using marker genes. Users can also provide their own pre-defined spots as input for the model.

# Preparing for MGGNN

Before using MGGNN, update the file paths in MGGNN.py to match your local directory structure.

# Running MGGNN

To run the model, use the following command:

python MGGNN.py \
    --SAMPLE 'control' \
    --N_CLUSTERS 7 \
    --N_EPOCHS 100 \
    --W1 10 \
    --W2 1 \
    --EARLY_STOPPING 0.9

Parameters:
	•	--SAMPLE: Name of the sample to process (e.g., 'control').
	•	--N_CLUSTERS: Number of clusters to identify.
	•	--N_EPOCHS: Number of training epochs.
	•	--W1: Weight for reconstruction loss.
	•	--W2: Weight for contrastive loss.
	•	--EARLY_STOPPING: Early stopping threshold.
