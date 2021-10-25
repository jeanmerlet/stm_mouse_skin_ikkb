import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scprep
import time

data_dir = "/Users/6j9/projects/mouse/data/"
mcl_file = "/Users/6j9/projects/mouse/irf-mcl/mcl_clusters/mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.2.txt"

experiments = ['FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848']
m = []
for experiment in experiments:
  matrix_dir = data_dir + experiment + '/filtered_feature_bc_matrix'
  mat = scprep.io.load_10X(matrix_dir, sparse=True, gene_labels='both')
  m.append(mat)

experiment_names = ['86846', '86847', '86848', '86849', '92848']
matrix, labels = scprep.utils.combine_batches(m, experiment_names, append_to_cell_names=True)
del m

#matrix = scprep.filter.filter_rare_genes(matrix, min_cells=1)
#matrix = scprep.normalize.library_size_normalize(matrix)
#matrix = scprep.transform.sqrt(matrix)

cluster_count = 0 
mcl_cluster_ids = {}
with open(mcl_file, 'rt') as in_file:
  for line in in_file:
    cluster_count += 1
    line_list = line.strip().split('\t')
    for item in line_list:
      mcl_cluster_ids[item] = cluster_count

for barcode, cluster_number in mcl_cluster_ids.items():
  labels[barcode] = cluster_number

cluster_matrices = []
for i in range(cluster_count):
  cluster_matrix = matrix.loc[labels.index[labels == i+1]]
  cluster_matrix = scprep.filter.filter_rare_genes(cluster_matrix, min_cells=1)
  cluster_matrices.append(cluster_matrix)
