import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.decomposition
import sklearn.manifold
import phate
import scprep
import time

data_dir = "/Users/6j9/projects/mouse/data/"
mcl_file = "/Users/6j9/projects/mouse/phate/mcl/mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.25.txt"

experiments = ['FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848']
m = []
for experiment in experiments:
  matrix_dir = data_dir + experiment + '/filtered_feature_bc_matrix'
  mat = scprep.io.load_10X(matrix_dir, sparse=True, gene_labels='both')
  m.append(mat)

experiment_names = ['86846', '86847', '86848', '86849', '92848']
matrix, labels = scprep.utils.combine_batches(m, experiment_names, append_to_cell_names=True)
del m

matrix = scprep.filter.filter_rare_genes(matrix, min_cells=10)
#matrix = scprep.normalize.library_size_normalize(matrix)
#matrix = scprep.transform.sqrt(matrix)

start = time.time()
pca_operator = sklearn.decomposition.PCA(n_components=2)
Y_pca = pca_operator.fit_transform(np.array(matrix))
end = time.time()
print("Embedded PCA in {:.2f} seconds.".format(end-start))

start = time.time()
pca_operator = sklearn.decomposition.PCA(n_components=100)
tsne_operator = sklearn.manifold.TSNE(n_components=2)
Y_tsne = tsne_operator.fit_transform(pca_operator.fit_transform(np.array(matrix)))
end = time.time()
print("Embedded t-SNE in {:.2f} seconds.".format(end-start))

start = time.time()
phate_operator = phate.PHATE(n_jobs=-2)
Y_phate = phate_operator.fit_transform(matrix)
end = time.time()
print("Embedded PHATE in {:.2f} seconds.".format(end-start))

cluster_count = 0 
mcl_cluster_ids = {}
with open(mcl_file, 'rt') as in_file:
  for line in in_file:
    cluster_count += 1
    line_list = line.strip().split('\t')
    for item in line_list:
      mcl_cluster_ids[item] = cluster_count

# save each cluster as a triplot projection
base_filepath = '/Users/6j9/projects/mouse/phate/out/projections/cutoff0.04_inf1.25/triplot_proj_0.04cutoff_1.25inf_clust'
#for i in range(cluster_count):
for i in range(1):
  cluster_labels = labels.copy(deep = True)
  for barcode, cluster_number in mcl_cluster_ids.items():
    if cluster_number == i + 1:
      cluster_labels[barcode] = cluster_number
  cluster_labels[cluster_labels != i + 1] = 0

  cluster_labels = cluster_labels.astype(int).to_numpy()

  fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(16, 5)) 
  #plot PCA
  title = 'PCA - cluster ' + str(i + 1)
  scprep.plot.scatter2d(Y_pca[np.nonzero(cluster_labels == 0)[0], :], label_prefix="PC", title=title, s=10,
                        ticks=False, marker='o', edgecolors='black', facecolors='none', ax=ax1, legend=False)
  scprep.plot.scatter2d(Y_pca[np.nonzero(cluster_labels == i+1)[0], :], label_prefix="PC", title=title, s=10,
                        ticks=False, marker='o', ax=ax1)
  #plot tSNE
  title = 't-SNE - cluster ' + str(i + 1)
  scprep.plot.scatter2d(Y_tsne[np.nonzero(cluster_labels == 0)[0], :], label_prefix="t-SNE", title=title, s=10,
                        ticks=False, marker='o', edgecolors='black', facecolors='none', ax=ax2, legend=False)
  scprep.plot.scatter2d(Y_tsne[np.nonzero(cluster_labels == i+1)[0], :], label_prefix="t-SNE", title=title, s=10,
                        ticks=False, marker='o', ax=ax2)
  #plot PHATE
  title = 'PHATE - cluster ' + str(i + 1)
  scprep.plot.scatter2d(Y_phate[np.nonzero(cluster_labels == 0)[0], :], label_prefix="PHATE", title=title, s=10,
                        ticks=False, marker='o', edgecolors='black', facecolors='none', ax=ax3, legend=False)
  scprep.plot.scatter2d(Y_phate[np.nonzero(cluster_labels == i+1)[0], :], label_prefix="PHATE", title=title, s=10,
                        ticks=False, marker='o', ax=ax3)
  plt.tight_layout()
  filename = base_filepath + str(i + 1) + '.png'
  plt.show()
  #plt.savefig(filename)
  #plt.close()
  print(f"cluster {str(i + 1)} saved...")
