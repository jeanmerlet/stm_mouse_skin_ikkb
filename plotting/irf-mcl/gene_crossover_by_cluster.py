import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scprep
import re

prefix = "/Users/6j9/projects/mouse"
data_dir = prefix + "/data/"
inf = "1.30"
mcl_path = prefix + "/irf-mcl/mcl_clusters/mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_" + inf + ".txt"

experiments = ['FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848']
m = []
for experiment in experiments:
  print(experiment)
  matrix_dir = data_dir + experiment + '/filtered_feature_bc_matrix'
  mat = scprep.io.load_10X(matrix_dir, sparse=True, gene_labels='both')
  m.append(mat)

experiment_names = ['86846', '86847', '86848', '86849', '92848']
#matrix, labels = scprep.utils.combine_batches(m, experiment_names, append_to_cell_names=True)
matrix, labels = scprep.utils.combine_batches(m, experiment_names, append_to_cell_names=True)
print("batches combined...")
del m

gene_symbols = np.full(len(matrix.columns), "", dtype='U30')
for i, gene in enumerate(matrix.columns.values):
  gene_symbols[i] =  (re.search("(.+) \(", gene).group(1)).upper()

matrix.columns = gene_symbols
#matrix_ln = scprep.normalize.library_size_normalize(matrix)
genes = ["PRRX1", "IKBKB"]
matrix = matrix.loc[:, genes]
#matrix_ln = scprep.normalize.library_size_normalize(matrix.T).T

# 1.30
cluster_nums = [1, 2, 4, 7, 8, 11, 13, 15, 20, 23, 25, 31, 33, 34, 39, 38, 44, 45, 47]

# 1.20
#cluster_nums = [2, 49, 45, 32, 52, 3, 60, 41, 26, 12, 10, 6, 14, 63, 15, 56, 68, 7, 33, 44, 9, 66, 59, 17, 53, 62, 67, 25, 5, 72, 27, 70, 65, 73]

cluster_count = 0
cluster_labels = labels.copy()
keep = []
with open(mcl_path, 'rt') as in_file:
  for line in in_file:
    cluster_count += 1
    if cluster_count in cluster_nums:
      line_list = line.strip().split('\t')
      for item in line_list:
        cluster_labels[item] = cluster_count
        keep.append(item)

cluster_labels = cluster_labels[cluster_labels.index.isin(keep)]
idx = np.in1d(matrix.index, cluster_labels.index)
matrix = matrix[idx]

cre_pos_idx = matrix.index.str.match(r".*_\d{2}84[67]")
cre_pos = matrix[cre_pos_idx]

for cluster in cluster_nums:
  idx = cluster_labels == cluster
  cluster_mat = matrix[idx]['PRRX1']
  cluster_mat = cluster_mat.sparse.to_dense()
  median = cluster_mat.median()
  print(f"cluster {cluster}: {median}")
  if median >= 1:
    print("MEOW")
  cluster_mat = cluster_mat.iloc[np.nonzero(cluster_mat)[0]]
  median = cluster_mat.median()
  print(f"nnz cluster {cluster}: {median}")
  print("")

for gene in genes:
  break
  gene_col = matrix[gene]
  nnz_cols = gene_col[np.nonzero(gene_col)[0]]
  nnz_cols = np.sort(nnz_cols)
  cutoff_idx = round(0.9 * len(nnz_cols))

gene = 'PRRX1'
figs, ax = plt.subplots(1, 1, figsize=(10, 5))
title = 'PRRX1 vs. IKBKB in Cre+ Fibroblasts'
ax.set_title(title)
ax.set_xlabel(gene + " UMI count")
ax.set_ylabel("Number of cells")
data = cre_pos
data['IKBKB_count'] = np.NaN
umi_counts = np.unique(data[gene])
nnz_cre_pos = cre_pos.iloc[np.nonzero(cre_pos[gene])[0], :]
nnz_ikbkb_count = np.count_nonzero(nnz_cre_pos['IKBKB'])
for n in umi_counts:
  print(np.count_nonzero(nnz_ikbkb_count['IKBKB']))
  continue
  if n == 0:
    continue
  else:
    ikbkb_idx = data[gene] == n
    if n == 1:
      del_ikbkb_count = nnz_ikbkb_count - np.count_nonzero(data[data[gene] == n]['IKBKB'])
    else:
      del_ikbkb_count = del_ikbkb_count - np.count_nonzero(data[data[gene] == n]['IKBKB'])
    print(n)
    print(del_ikbkb_count)
    print("")
  ikbkb_count = del_ikbkb_count
  k = 0
  for bc in data.index:
    if k >= ikbkb_count:
      break
    if ikbkb_idx[bc]:
      data.loc[bc, 'IKBKB_count'] = n
      k += 1

x_tick_labels = np.arange(0, max(umi_counts) + 1)
ax.set_xticks(x_tick_labels)
bins = np.arange(0, max(umi_counts) + 1) - 0.5
ax.hist(data['PRRX1'], bins=bins, label='Cells with n UMIs of PRRX1', alpha=0.5)
ax.hist(data['IKBKB_count'], bins=bins, label='Cumulative IKBKB-deleted Cells', cumulative=True, alpha=0.5)

ax.legend()
plt.show()
