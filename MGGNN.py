# mamba activate MGGNN
# cd /ssd/haoran/MGGNN/code
# .libPaths(): /ssd/haoran/modules/mambaforge/envs/MGGNN/lib/R/library
# R_HOME: /ssd/haoran/modules/mambaforge/envs/MGGNN/lib/R

import os
import torch
import pandas as pd
import scanpy as sc
from sklearn import metrics
import multiprocessing as mp

import numpy as np
import math
import GraphST_modified

import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--SAMPLE', type = str)
parser.add_argument('--N_CLUSTERS', type = int)
parser.add_argument('--N_EPOCHS', type = int)
parser.add_argument('--W1', type = int)
parser.add_argument('--W2', type = int)
parser.add_argument('--EARLY_STOPPING', type = float)
args = parser.parse_args()

SAMPLE = args.SAMPLE
N_CLUSTERS = args.N_CLUSTERS
N_EPOCHS = args.N_EPOCHS
W1 = args.W1
W2 = args.W2
EARLY_STOPPING = args.EARLY_STOPPING

# SAMPLE = 'control'
# N_CLUSTERS = 7
# N_EPOCHS = 100
# W1 = 10
# W2 = 1
# EARLY_STOPPING = 0.9



# Run device, by default, the package is implemented on 'cpu'. We recommend using GPU.
# device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
device = torch.device(f'cuda:{torch.cuda.current_device()}' if torch.cuda.is_available() else 'cpu')

# the location of R, which is necessary for mclust algorithm. Please replace the path below with local R installation path
os.environ['R_HOME'] = '/ssd/haoran/modules/mambaforge/envs/MGGNN/lib/R' # R RHOME

# read data
file_fold = '/ssd/haoran/MGGNN/data/Heme/' + str(SAMPLE) #please replace 'file_fold' with the download path
adata = sc.read_visium(file_fold, count_file = 'filtered_feature_bc_matrix.h5', load_images = True)
adata.var_names_make_unique()

# define model
model = GraphST_modified.GraphST(adata, device = device)

# train model
adata = model.train()

# set radius to specify the number of neighbors considered during refinement
radius = 50

tool = 'mclust' # mclust, leiden, and louvain

# clustering
from utils import clustering

if tool == 'mclust':
    clustering(adata, N_CLUSTERS, radius=radius, method=tool, refinement=True) # For DLPFC dataset, we use optional refinement step.
elif tool in ['leiden', 'louvain']:
    clustering(adata, N_CLUSTERS, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

# add ground_truth
df_meta = pd.read_csv('/ssd/haoran/MGGNN/data/Heme/' + SAMPLE + '/' + SAMPLE + '_anatomy.csv', index_col=0)
df_meta['anatomy'][df_meta['anatomy'] == 'hypothalamuis'] = 'hypothalamus'
df_meta_layer = df_meta['anatomy']
adata.obs['ground_truth'] = df_meta_layer.values

df_f = adata.obs[['array_row', 'array_col', 'domain', 'ground_truth']]
# df_f.to_csv("/ssd/haoran/MGGNN/results/MGGNN_" + SAMPLE + ".csv")

# filter out NA nodes
idx_notnull = ~pd.isnull(adata.obs['ground_truth'])
adata = adata[idx_notnull]

# calculate metric ARI
ARI = metrics.adjusted_rand_score(adata.obs['domain'], adata.obs['ground_truth'])
adata.uns['ARI'] = ARI



#####
##### fine_tuning
#####

x = adata.obsm['feat'].copy()
y = adata.obs['ground_truth'].copy()
# idx_train = np.random.choice(x.shape[0], math.floor(x.shape[0]*PERCENTAGE), replace=False)
sample_train = pd.read_csv('/ssd/haoran/MGGNN/data/Heme_selected_spots/' + SAMPLE + '_threshold.csv', index_col = 0)
sample_train = sample_train.dropna()

tmp_list1 = adata.obs.index.tolist()
tmp_list2 = sample_train.index.tolist()
idx_train = [tmp_list1.index(i) for i in tmp_list2]
y_train = sample_train.rank_prediction

# Y: layer to num
layer = np.sort(y.unique())
layer_idx = {name: id for id, name in enumerate(layer)}
f = lambda value: layer_idx[value]
y_num = np.array([f(value) for value in y])
y_train_num = np.array([f(value) for value in y_train])



# self.model.train()
# radius = 50
from spatialentropy import altieri_entropy

if len(sample_train) / len(df_f) > 0.2:
    df_f = adata.obs[['array_row', 'array_col', 'domain', 'ground_truth']]
else:
    pred_fine_tuning, hiden_feat, emb, epochs = model.fine_tuning(x, y_num, idx_notnull, idx_train, y_train_num, N_CLUSTERS, N_EPOCHS, W1, W2, EARLY_STOPPING)
    pred_fine_tuning = pred_fine_tuning.argmax(1).detach().cpu().numpy()
    hiden_feat = hiden_feat.detach().cpu().numpy()
    emb = emb.detach().cpu().numpy()
    adata.obsm['emb'] = hiden_feat
    clustering(adata, N_CLUSTERS, radius=radius, method=tool, refinement=True)
    df_f = adata.obs[['array_row', 'array_col', 'domain', 'ground_truth']]


#
df_f.to_csv("/ssd/haoran/MGGNN/results/MGGNN_" + SAMPLE + ".csv")

#
print('\n')
print('#' * 50)
print('SAMPLE:\t\t\t\t\t', SAMPLE)
print('EARLY_STOPPING:\t\t\t\t', EARLY_STOPPING)
print('MAX_EPOCHS:\t\t\t\t', N_EPOCHS)
print('EPOCHS:\t\t\t\t\t', epochs)
