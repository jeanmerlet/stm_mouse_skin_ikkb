import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.decomposition
import sklearn.manifold
import scprep
import time
import os
import re
import glob

np.random.seed(888)
random_state = 888

sample_ids_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/mesenchymal_count_mtxs/dermis_ids.txt'
matrix_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/prrx1_matrices/combined_all_prrx1_pos.tsv'
batch_matrix_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/prrx1_matrices/batch_combined_all_prrx1_pos.tsv'
clusters_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/scripts/human_ad/tsne_clustering/reynolds/clusters_20dim.tsv'

mtx = scprep.io.load_tsv(matrix_path, cell_axis='column')

# load in the batch corrected matrix for visualization
batch_mtx = scprep.io.load_tsv(batch_matrix_path, cell_axis='column')

# t-sne stuff
start = time.time()
pca_operator = sklearn.decomposition.PCA(n_components=100, random_state=random_state)
tsne_operator = sklearn.manifold.TSNE(n_components=2, random_state=random_state)
Y_tsne = tsne_operator.fit_transform(pca_operator.fit_transform(batch_mtx.values))
end = time.time()
print('Embedded t-SNE in {:.2f} seconds.'.format(end-start))

# plot by healthy vs. AD
cmap = ['lightcoral', 'dodgerblue']

labels = []
for bc in mtx.index.values:
    acc_id = re.search('_(.*)$', bc).groups()[0]
    if '809' in acc_id:
        labels.append('AD')
    else:
        labels.append('H')

fig, ax = plt.subplots(figsize=(8, 8))
scprep.plot.scatter2d(Y_tsne, label_prefix='t-SNE', title='PRRX1+ Mesenchymal Cells by Condition', legend_anchor=(1, 1),
                      c=labels, ticks=False, cmap=cmap, ax=ax, s=8, fontsize=14, discrete=True)
plt.tight_layout()
out_path = f'/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/plots/new/reynolds_prrx1-pos_H_vs_AD.pdf'
plt.savefig(out_path)
plt.clf()

# plot individual genes
cmap = ['lightgrey', 'red']
genes = ['ACTA2', 'MCAM', 'DCN', 'PDGFRA', 'COL1A2', 'COL3A1', 'CCL26', 'CCL11', 'CEBPB']

for gene in genes:
    labels = mtx.loc[:, mtx.columns.str.contains(gene)]
    if labels.shape[1] > 1:
        labels = labels.iloc[:, 0].values
    fig, ax = plt.subplots(figsize=(10, 8))
    scprep.plot.scatter2d(Y_tsne, label_prefix='t-SNE', title=f'{gene} in PRRX1+ Mesenchymal Cells',
                          c=labels, ticks=False, cmap=cmap, ax=ax, s=8, fontsize=14)
    plt.tight_layout()
    out_path = f'/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/plots/new/reynolds_prrx1-pos_H_vs_AD_{gene}.pdf'
    plt.savefig(out_path)
    plt.clf()

# for clusters
clusters = pd.read_csv(clusters_path, sep='\t', index_col=0, header=None)
labels = []
for bc in mtx.index.values:
    labels.append(clusters.loc[bc].values[0])

cmap = None

fig, ax = plt.subplots(figsize=(12, 8))
scprep.plot.scatter2d(Y_tsne, label_prefix='t-SNE', title='kNN Clusters', legend_anchor=(1, 1),
                      c=labels, ticks=False, cmap=cmap, ax=ax, s=10, fontsize=14, discrete=True)
plt.tight_layout()
out_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/plots/new/reynolds_prrx1-pos_H_vs_AD_clusters.pdf'
plt.savefig(out_path)
plt.clf()
