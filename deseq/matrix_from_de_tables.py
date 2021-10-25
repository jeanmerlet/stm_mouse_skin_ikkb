import pandas as pd
import numpy as np
import os
import re

inf = "1.30"
prefix = "/Users/6j9/projects/mouse/deseq"
de_dir = os.path.join(prefix, "de_results", inf, "ikbkb-neg_tables")
matrix_name = "upenn_mouse_lfc_pairwise_matrix_inf-" + inf + "_padj-0.01_ikbkb-neg.tsv"
matrix_path = os.path.join(prefix, "matrix", matrix_name)

table_paths = []
for root, dirs, files in os.walk(de_dir):
  for table_path in files:
    if ".tsv" in table_path:
      table_paths.append(os.path.join(root, table_path))

p_cutoff = 0.01
lfcs = []

for i, path in enumerate(table_paths):
  table = pd.read_table(path, sep="\t", index_col=0)
  print(i)
  if i == 0:
    genes = table.index
    data = np.zeros((len(genes), len(table_paths)))
    matrix = pd.DataFrame(data, index=genes)
  head, tail = os.path.split(path)
  tail = re.search("(\d+_vs_\d+)", tail).group(1)
  colname = {i: tail}
  matrix = matrix.rename(columns=colname)
  table.loc[table['padj'].isna(), 'log2FoldChange'] = 0
  table.loc[table['padj'] > p_cutoff, 'log2FoldChange'] = 0
  matrix[colname[i]] = table['log2FoldChange']
  lfcs.append(np.count_nonzero(matrix[colname[i]] != 0))

print(f"max: {max(lfcs)}")
print(f"min: {min(lfcs)}")
print(f"mean: {np.mean(lfcs)}")
print(f"median: {np.median(lfcs)}")

print(matrix_path)
matrix = matrix.fillna(0)
matrix.to_csv(matrix_path, sep="\t")
