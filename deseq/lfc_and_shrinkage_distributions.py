import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
import re
from scipy import stats

inf = "1.25"
prefix = "/Users/6j9/projects/mouse/deseq"
de_dir = os.path.join(prefix, "results", inf, "new_tables")
matrix_name = "upenn_mouse_lfc_pairwise_matrix_inf-" + inf + "_alpha-0.01_ikbkb-neg.tsv"
matrix_path = os.path.join(prefix, "matrix", matrix_name)

table_paths = []
for root, dirs, files in os.walk(de_dir):
  for table_path in files:
    if ".tsv" in table_path:
      table_paths.append(os.path.join(root, table_path))

cre_pos = [1, 38, 39, 23, 45, 44, 33, 11]
cre_neg = [4, 47, 34, 20, 15, 8, 7, 25, 2, 31, 13]

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

def plot_distribution(data, title, x_label):
  ax = sns.distplot(data, bins=160, kde=False)
  ax.set(xlabel=x_label, title=title)
  plt.show()

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

matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
data = matrix
#data = subset_to_clusters(matrix, cre_pos)
#data = subset_to_clusters(matrix, cre_neg, cre_pos)
data = data.values.flatten()
data = np.delete(data, np.where(data == 0))
title = 'lfc distribution - nonzero values with p-adjusted < 0.01'
x_label = 'lfc'

#data = get_column_from_deseq2_tables(table_paths, 'stat')
#title = 'wald statistic distribution - all'
#x_label = 'wald test statistic'

plot_distribution(data, title, x_label)
