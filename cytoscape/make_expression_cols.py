import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scprep
import re

prefix = "/Users/6j9/projects/mouse"
data_dir = prefix + "/data/"

experiments = ['FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848']
m = []
for experiment in experiments:
  print(experiment)
  matrix_dir = data_dir + experiment + '/filtered_feature_bc_matrix'
  mat = scprep.io.load_10X(matrix_dir, sparse=True, gene_labels='both')
  m.append(mat)

experiment_names = ['86846', '86847', '86848', '86849', '92848']
matrix, labels = scprep.utils.combine_batches(m, experiment_names, append_to_cell_names=True)
print("batches combined...")
del m

matrix = scprep.normalize.library_size_normalize(matrix)

gene_symbols = np.full(len(matrix.columns), "", dtype='U30')
for i, gene in enumerate(matrix.columns.values):
  symbol = re.search("(.+) \(", gene).group(1)
  gene_symbols[i] = symbol.upper()

#genes = ["PRRX1", "IKBKB", "CCL11", "CCL19", "CXCL1", "CXCL12"]
#genes = ["PRRX1", "IKBKB"]
genes = ['CXCL12', 'CCL11', 'CCL7']

out_cols = pd.DataFrame(index=matrix.index, columns=genes)
for gene in genes:
  idx = np.where(gene_symbols == gene)[0][0]
  out_cols[gene] = matrix.iloc[:, idx]
  nnz_cols = out_cols[gene][np.nonzero(out_cols[gene])[0]]
  nnz_cols = np.sort(nnz_cols)
  cutoff_idx = round(0.9 * len(nnz_cols))
  print(f"{gene} 90% cutoff: {nnz_cols[cutoff_idx]}")

#figs, axs = plt.subplots(2, 3)
#plt.clf()
#figs, ax = plt.subplots(1, 1)
#axs = [ax for list in axs for ax in list]
for i, gene in enumerate(genes):
  break
  title = gene + " Nonzero Expression"
  title = gene + " Expression"
  ax.set_title(title)
  ax.set_xlabel("UMIs per cell / total UMIs per cell * 10^4")
  ax.set_ylabel("Number of cells")
  data = out_cols[gene]
  data = data[np.nonzero(data)[0]]
  x_tick_labels = np.arange(round(max(data)) + 1)
  ax.set_xticks(x_tick_labels)
  #nnz_data = data[np.nonzero(data)[0]]
  #ax.hist(nnz_data, bins=50)
  #ax.axvline(9.73, color='r')
  #ax.legend(labels=['0.90 cutoff'])
  ax.hist(data, bins=len(x_tick_labels))
  #ax.set_xticks([], [])

#plt.tight_layout()
#plt.show()



out_path = prefix + "/cytoscape/coloring/chemokine_cols.tsv"
out_cols.to_csv(out_path, sep="\t")
