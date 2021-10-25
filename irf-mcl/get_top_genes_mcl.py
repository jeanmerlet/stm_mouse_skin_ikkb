import numpy as np
import pandas as pd
import scprep
import re

prefix = "/Users/6j9/projects/mouse"
data_dir = prefix + "/data/"

cutoff, inf = "0.04", "1.20"
mcl_suffix = "mouseUPenn_iRF-LOOP_normalized_" + cutoff + "_mcl_clusters_" + inf + ".txt"
mcl_file = prefix + "/irf-mcl/mcl_clusters/" + suffix
pearson_suffix = "irf-mcl_" + cutoff + "_cut_" + inf + "_inf_" + "cluster_pearson.tsv"
out_path = prefix + "/irf-mcl/cluster_pearson/" + suffix

experiments = ['FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848']
m = []
for experiment in experiments:
  print(experiment)
  matrix_dir = data_dir + experiment + '/filtered_feature_bc_matrix'
  mat = scprep.io.load_10X(matrix_dir, sparse=True, gene_labels='both')
  m.append(mat)

experiment_names = ['86846', '86847', '86848', '86849', '92848']
matrix, labels = scprep.utils.combine_batches(m, experiment_names, append_to_cell_names=True)
del m

matrix = scprep.filter.filter_rare_genes(matrix, min_cells=1)
matrix = scprep.normalize.library_size_normalize(matrix)
mito_genes = np.array([g.startswith('mt-') for g in matrix.columns])
mito_exp = matrix.loc[:, mito_genes].mean(axis=1)
matrix = scprep.filter.filter_values(matrix, values=mito_exp, percentile=95, keep_cells='below')
labels = labels[matrix.index]
matrix = scprep.transform.sqrt(matrix)

cluster_count = 0 
mcl_cluster_ids = {}
with open(mcl_file, 'rt') as in_file:
  for line in in_file:
    cluster_count += 1
    line_list = line.strip().split('\t')
    for item in line_list:
      if item in labels:
        labels[item] = cluster_count

num_genes = np.shape(matrix)[1]
mean_exp_vectors = np.zeros((num_genes, cluster_count))
ratio_gene_exp = np.zeros((num_genes, cluster_count))
for i in range(cluster_count):
  print(i)
  cluster_matrix = matrix.loc[labels.index[labels == i+1]]
  vector = np.mean(cluster_matrix, axis=0)
  ratios = np.count_nonzero(cluster_matrix, axis=0)
  mean_exp_vectors[:, i] = vector
  ratio_gene_exp[:, i] = ratios
pearson = np.corrcoef(mean_exp_vectors, rowvar=False)

useful_genes_idx = np.array([not g.startswith(('mt-', 'Rpl', 'Rps', 'Gm42418')) for g in matrix.columns])
useful_genes = matrix.columns[useful_genes_idx]
useful_exp_vectors = mean_exp_vectors[useful_genes_idx, :]
useful_ratios = ratio_gene_exp[useful_genes_idx, :]

def get_top_n_genes(cluster, n):
  idx = np.argsort(useful_exp_vectors[:, cluster - 1])
  n_idx = idx[-n:]
  return [useful_genes[n_idx], useful_exp_vectors[:, cluster - 1][n_idx]]

def get_top_n_ratios(cluster, n):
  idx = np.argsort(useful_ratios[:, cluster - 1])
  n_idx = idx[-n:]
  return [useful_genes[n_idx], useful_ratios[:, cluster - 1][n_idx]]

def test(cluster, n):
  top_genes = get_top_n_genes(cluster, n)
  a = top_genes[0]
  b = get_top_n_ratios(cluster, n)[0]
  idx = np.isin(a, b, invert=True)
  return [a[idx], top_genes[1][idx]]
  #return a[np.isin(a, b, invert=True)]

#get_top_n_genes(1, 20)

#gene_names = np.full(num_genes, "", dtype='U18')
#for i, gene in enumerate(matrix.columns):
#  gene_names[i] = gene[-19:-1]
#np.column_stack((gene_names, pearson))

np.savetxt(out_path, pearson, fmt="%s", delimiter="\t")
