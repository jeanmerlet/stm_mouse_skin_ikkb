import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.decomposition
import sklearn.manifold
import scprep
import time
import os
import re

np.random.seed(888)
random_state = 888

mtx_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/combined_he-preprocessed_human-ad_mtx.tsv'
clusters_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/scripts/human_ad/tsne_clustering/he/clusters_20dim.tsv'

mtx = pd.read_csv(mtx_path, sep='\t', index_col=0, header=0)
mtx = mtx.fillna(0)
mtx = mtx.transpose()

# loading in clusters
clusters = pd.read_csv(clusters_path, sep='\t', index_col=0, header=None)
labels = []
for bc in mtx.index.values:
    labels.append(clusters.loc[bc].values[0])
        
# t-sne stuff
start = time.time()
pca_operator = sklearn.decomposition.PCA(n_components=100, random_state=random_state)
tsne_operator = sklearn.manifold.TSNE(n_components=2, random_state=random_state)
#Y_tsne = tsne_operator.fit_transform(pca_operator.fit_transform(mtx.values))
Y_tsne = tsne_operator.fit_transform(pca_operator.fit_transform(prrx1_mtx.values))
end = time.time()
print('Embedded t-SNE in {:.2f} seconds.'.format(end-start))

# plot clusters onto t-SNE
cmap = None

fig, ax = plt.subplots(figsize=(12, 8))
scprep.plot.scatter2d(Y_tsne, label_prefix='t-SNE', title='KNN Clusters', legend_anchor=(1, 1),
                      c=labels, ticks=False, cmap=cmap, ax=ax, s=10, fontsize=14, discrete=True)
plt.tight_layout()
out_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/plots/clusters_20dim.png'
plt.savefig(out_path, dpi=300)

# case control

bc_to_condition_map = {'SRR11396159': 'LS',
                       'SRR11396160': 'LS',
                       'SRR11396161': 'NL',
                       'SRR11396162': 'H',
                       'SRR11396163': 'LS',
                       'SRR11396164': 'H',
                       'SRR11396165': 'LS',
                       'SRR11396166': 'H',
                       'SRR11396167': 'H',
                       'SRR11396168': 'H',
                       'SRR11396169': 'NL',
                       'SRR11396170': 'H',
                       'SRR11396171': 'H',
                       'SRR11396172': 'NL',
                       'SRR11396173': 'NL',
                       'SRR11396174': 'NL',
                       'SRR11396175': 'H'}

bc_to_sample_map = {'S1': 'SRR11396159',
                   'S2': 'SRR11396160',
                   'S3': 'SRR11396161',
                   'S4': 'SRR11396162',
                   'S5': 'SRR11396163',
                   'S6': 'SRR11396164',
                   'S7': 'SRR11396165',
                   'S8': 'SRR11396166',
                   'S9': 'SRR11396167',
                   'S10': 'SRR11396168',
                   'S11': 'SRR11396169',
                   'S12': 'SRR11396170',
                   'S13': 'SRR11396171',
                   'S14': 'SRR11396172',
                   'S15': 'SRR11396173',
                   'S16': 'SRR11396174',
                   'S17': 'SRR11396175'}

labels = []
for bc in prrx1_mtx.index.values:
    bc_id = re.search('^(S\d+)_', bc).groups()[0]
    labels.append(bc_to_condition_map[bc_to_sample_map[bc_id]])

cmap = ['dodgerblue', 'lightcoral']

fig, ax = plt.subplots(figsize=(8, 8))
scprep.plot.scatter2d(Y_tsne, label_prefix='t-SNE', title='PRRX1+ Mesenchymal Cells by Condition', legend_anchor=(1, 1),
                      c=labels, ticks=False, cmap=cmap, ax=ax, s=8, fontsize=14, discrete=True)
plt.tight_layout()
out_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/plots/he_prrx1_pos_condition.pdf'
plt.savefig(out_path)
#plt.savefig(out_path, dpi=300)

# mesenchymal cells are composed of clusters:
# 3, 4, 5, 8
# where 3 is pericytes and 4, 5, 8 are fibroblursts

# subsetting to prrx1+ cells
prrx1_mtx_idx = mtx.loc[:, 'PRRX1'] > 0
prrx1_mtx = mtx.loc[prrx1_mtx_idx, :]

# plot individual genes from a list onto t-SNE
cmap = ['dodgerblue', 'lightcoral']
cmap = ['mediumblue', 'white', 'red']
genes = ['ACTA2', 'MCAM', 'CCL11', 'CCL26']
genes = ['PRRX1', 'DCN', 'PDGFRA']
genes = ['CCL26', 'CCL11', 'CEBPB']

for gene in genes:
    bcs = mtx.loc[:, gene]
    bcs = bcs[bcs > 0].index.values
    #labels = []
    labels = prrx1_mtx.loc[:, gene].values
    for bc in mtx.index.values:
        break
        if bc in bcs:
            labels.append(f'{gene}_pos')
        else:
            labels.append(f'{gene}_neg')
    fig, ax = plt.subplots(figsize=(10, 8))
    scprep.plot.scatter2d(Y_tsne, label_prefix='t-SNE', title=f'{gene} in PRRX1+ Mesenchymal Cells',
                          c=labels, ticks=False, cmap=cmap, ax=ax, s=8, fontsize=14)
    plt.tight_layout()
    out_path = f'/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/plots/he_prrx1_pos_{gene}.pdf'
    #plt.savefig(out_path, dpi=300)
    plt.savefig(out_path)
    plt.clf()
