import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.decomposition
import sklearn.manifold
import phate
import umap
import scprep
import time

np.random.seed(888)
random_state = 888

data_dir = "/Users/6j9/projects/mouse/data/"
mtx_tsv_path = '/Users/6j9/projects/mouse/data/combined_mtx.tsv'

cell_types_path = '/Users/6j9/projects/mouse/data/top8_v2.tsv'
cell_types = pd.read_csv(cell_types_path, sep='\t', header=None, index_col=0)
fibro_idx = (cell_types == 'Fibroblasts').values
fibro_bcs = cell_types.loc[fibro_idx, :].index.values

matrix = pd.read_csv(mtx_tsv_path, sep='\t', header=0, index_col=0)
matrix = matrix.transpose()
matrix = scprep.filter.filter_rare_genes(matrix, min_cells=10)
matrix = scprep.normalize.library_size_normalize(matrix)
matrix = scprep.transform.sqrt(matrix)

fibro_mtx = matrix.loc[fibro_bcs, :]

labels = ['WT'] * 6599 + ['cKO'] * 4485
fibro_labels = ['WT'] * 3609 + ['cKO'] * 1527

start = time.time()
pca_operator = sklearn.decomposition.PCA(n_components=2, random_state=random_state)
Y_pca = pca_operator.fit_transform(np.array(matrix))
end = time.time()
print("Embedded PCA in {:.2f} seconds.".format(end-start))

start = time.time()
pca_operator = sklearn.decomposition.PCA(n_components=100, random_state=random_state)
tsne_operator = sklearn.manifold.TSNE(n_components=2, random_state=random_state)
#Y_tsne = tsne_operator.fit_transform(pca_operator.fit_transform(np.array(matrix)))
Y_tsne = tsne_operator.fit_transform(pca_operator.fit_transform(fibro_mtx.values))
end = time.time()
print("Embedded t-SNE in {:.2f} seconds.".format(end-start))

start = time.time()
reducer = umap.UMAP()
#UMAP = reducer.fit_transform(matrix.values)
UMAP = reducer.fit_transform(fibro_mtx.values)
end = time.time()
print("Embedded UMAP in {:.2f} seconds.".format(end-start))

#start = time.time()
#phate_operator = phate.PHATE(n_jobs=-2)
#Y_phate = phate_operator.fit_transform(matrix)
#end = time.time()
#print("Embedded PHATE in {:.2f} seconds.".format(end-start))

#fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(16, 5))
plt.clf()
# for cre+ / cre-
cmap = ['dodgerblue', 'lightcoral']
# for gene expression
cmap = ['mediumblue', 'white', 'darkred']
gene = 'CCL11'
labels = fibro_mtx.loc[:, gene].values
# to check median:
# np.median(np.array(labels)[np.array(labels) > 0])
labels = [label if label < 6 else 6 for label in labels]

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(18, 6))
# plot PCA
#scprep.plot.scatter2d(Y_pca, label_prefix="PC", title="PCA",
#                      c=labels, ticks=False, cmap=cmap, ax=ax1)

# plot UMAP
scprep.plot.scatter2d(UMAP, label_prefix="UMAP", title="UMAP",
                      c=labels, ticks=False, cmap=cmap, ax=ax1, s=6)
# plot tSNE
scprep.plot.scatter2d(Y_tsne, label_prefix="t-SNE", title="t-SNE", legend=False,
                      c=labels, ticks=False, cmap=cmap, ax=ax2, s=6)
fig.legend()
# plot PHATE
#scprep.plot.scatter2d(Y_phate, label_prefix="PHATE", title="PHATE", legend=False,
#                      c=labels, ticks=False, cmap='Spectral', ax=ax3)

plt.tight_layout()
plt.show()


# need to do prrx1 high...again



