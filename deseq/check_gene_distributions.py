import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import scprep
import os
import re
from scipy import stats

inf = "1.20"
prefix = "/Users/6j9/projects/mouse"
mcl_suffix = "mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_" + inf + ".txt"
mcl_file = os.path.join(prefix, "irf-mcl/mcl_clusters", mcl_suffix)

experiments = ['FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848']
m = []
for experiment in experiments:
  print(experiment)
  matrix_dir = os.path.join(prefix, "data", experiment, "filtered_feature_bc_matrix")
  mat = scprep.io.load_10X(matrix_dir, sparse=True, gene_labels='both')
  m.append(mat)

experiment_names = ['86846', '86847', '86848', '86849', '92848']
matrix, labels = scprep.utils.combine_batches(m, experiment_names, append_to_cell_names=True)
del m

gene_symbols = []
for gene in matrix.columns:
  gene_symbols.append(re.search("(.+) \(.+", gene).group(1).upper())

matrix.columns = gene_symbols
# remove zero / low-count genes
matrix = matrix.loc[:, matrix.sum(axis=0) > 10]

# 1.20 fibro clusters
cre_pos = [5, 7, 9, 15, 17, 25, 27]
cre_neg = [2, 3, 6, 10, 14, 32, 33]
clusters = cre_neg

# 1.25 fibro clusters
cre_pos = [7, 8, 11, 17, 22, 23, 28, 35, 38, 44, 51]
cre_neg = [2, 3, 4, 9, 10, 25, 26, 27, 41, 45, 46, 52, 53, 54]
#clusters = cre_neg

# 1.30 fibro clusters
cre_pos = [1, 11, 23, 33, 38, 39, 44, 45]
cre_neg = [2, 4, 7, 8, 13, 15, 20, 25, 31, 34, 47]
#clusters = cre_neg

def plot_gene_distribution_by_cluster(matrix, clusters, gene):
  cluster_bcs = {}
  cluster_num = 0
  with open(mcl_file, 'rt') as in_file:
    for line in in_file:
      cluster_num += 1
      if cluster_num in clusters:
        cluster_bcs[cluster_num] = line.strip().split('\t')
  n_ax_rows, n_ax_cols = 2, 4
  figs, axes = plt.subplots(n_ax_rows, n_ax_cols)
  ax_x_idx, ax_y_idx = 0, 0
  for key in cluster_bcs.keys():
    data = matrix.loc[cluster_bcs[key], gene]
    #ax = sns.distplot(data, bins=np.arange(0, 10, 1), kde=False, ax=axes[ax_y_idx, ax_x_idx])
    ax = sns.distplot(data, bins=np.arange(0, 10, 1), kde=False, norm_hist=True, ax=axes[ax_y_idx, ax_x_idx])
    ax.set(title=f"cluster {key}")
    ax.set_xticks(np.arange(0, 10, 1))
    ax.set_yticks(np.arange(0.0, 1.1, 0.2))
    if ax_x_idx == n_ax_cols - 1:
      ax_y_idx += 1
      ax_x_idx = 0
    else:
      ax_x_idx += 1
  plt.tight_layout()
  plt.show()
  #return(cluster_bcs)

plot_gene_distribution_by_cluster(matrix, clusters, 'PRRX1')

#cluster_bcs = plot_gene_distribution_by_cluster(matrix, clusters, 'PRRX1')
#for key in cluster_bcs:
  #print(f"{key}: {len(cluster_bcs[key])}")

cluster_num = 0 
barcodes = []
with open(mcl_file, 'rt') as in_file:
  for line in in_file:
    cluster_num += 1
    if cluster_num in clusters:
      barcodes.append(line.strip().split('\t'))

barcodes = [barcode for cluster in barcodes for barcode in cluster]
cluster_matrix = matrix[np.in1d(matrix.index, barcodes)]

i = 0
for gene in cluster_matrix.columns:
  break
  ax = sns.distplot(cluster_matrix[gene], bins=160, kde=False)
  plt.show()
  i += 1
  k2, p = stats.normaltest(cluster_matrix[gene][cluster_matrix[gene] != 0])
  worst_p = 0.0
  if p > worst_p:
    worst_p = p
  if i % 1000 == 0:
    print(i)

#print(worst_p)

def subset_to_clusters(data, clusters_a, clusters_b):
#def subset_to_clusters(data, clusters):
  #num_clusters = len(clusters)
  #num_subcomparisons = int(num_clusters * (num_clusters - 1) / 2)
  #new_data = pd.DataFrame(data=np.zeros((np.shape(data)[0], num_subcomparisons)), index=data.index)
  num_subcomparisons = int(len(clusters_a) * len(clusters_b))
  new_data = pd.DataFrame(data=np.zeros((np.shape(data)[0], num_subcomparisons)), index=data.index)
  i = 0
  for column in data.columns:
    cluster_numbers = re.match("(\d+)_vs_(\d+)", column).groups()
    a, b = cluster_numbers[0], cluster_numbers[1]
    #if int(a) in clusters and int(b) in clusters:
    if (int(a) in clusters_a and int(b) in clusters_b) or (int(a) in clusters_b and int(b) in clusters_a):
      new_data = new_data.rename(columns={i: column})
      i += 1
      new_data[column] = data[column]
  return new_data

def plot_distribution(data, title=None, x_label=None):
  ax = sns.distplot(data, bins=160, kde=False)
  #k2, p = stats.normaltest(data)
  #print("p = {:g}".format(p))
  #ax = sns.distplot(data, kde=False)
  #ax = sns.distplot(data)
  #ax.set(xlabel=x_label, title=title)
  #ax.set_xticks(np.arange(-40, 40, 2))
  plt.show()

#print(cluster_matrix)
#plot_distribution(cluster_matrix['PRRX1'])

def get_column_from_deseq2_tables(table_paths, tablecol, fdr_cutoff=False):
  for i, path in enumerate(table_paths):
    print(i)
    table = pd.read_table(path, sep="\t", index_col=0)
    if i == 0:
      genes = table.index
      data = np.zeros((len(genes), len(table_paths)))
      matrix = pd.DataFrame(data, index=genes)
    head, tail = os.path.split(path)
    tail = re.search("(\d+_vs_\d+)", tail).group(1)
    colname = {i: tail}
    matrix = matrix.rename(columns=colname)
    if fdr_cutoff:
      na_idx = table['padj'].isna()
      table = table[~na_idx]
      idx = table['padj'] < fdr_cutoff
      lfc = table[idx][tablecol]
      matrix[colname[i]] = lfc[matrix.index].values
      matrix[colname[i]][~matrix.index.isin(lfc.index)] = 0 
    else:
      lfc = table[tablecol]
      matrix[colname[i]] = lfc[matrix.index].values
  return matrix.values.flatten()

#matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
#data = subset_to_clusters(matrix, cre_pos)
#data = subset_to_clusters(matrix, cre_neg, cre_pos)
#data = data.values.flatten()
#data = np.delete(data, np.where(data == 0))

#data = get_column_from_deseq2_tables(table_paths, 'stat')
#print(data)

#title = 'lfc distribution - nonzero values with p-adjusted < 0.01 - cre+ vs. cre- only'
title = 'wald statistic distribution - all'
x_label = 'wald test statistic'
#plot_distribution(data, title, x_label)
