import numpy as np
import pandas as pd
import scprep
import os

prefix = "/Users/6j9/projects/mouse"
data_dir = os.path.join(prefix, "data")
out_path = os.path.join(prefix, "cytoscape/coloring/umi_counts.tsv")

experiments = ['FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848']
m = []
for experiment in experiments:
  print(experiment)
  matrix_dir = os.path.join(data_dir, experiment, "filtered_feature_bc_matrix")
  mat = scprep.io.load_10X(matrix_dir, sparse=True, gene_labels='both')
  m.append(mat)

experiment_names = ['86846', '86847', '86848', '86849', '92848']
matrix, labels = scprep.utils.combine_batches(m, experiment_names, append_to_cell_names=True)
print("batches combined...")
del m

umi_counts = matrix.sum(axis=1)
umi_counts.to_csv(out_path, sep="\t")
