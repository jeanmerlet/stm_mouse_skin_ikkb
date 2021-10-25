import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
import re
import sys
sys.setrecursionlimit(10**6)

inf = "1.30"
prefix = "/Users/6j9/projects/mouse/deseq"
matrix_name = "upenn_mouse_lfc_pairwise_matrix_inf-" + inf + "_padj-0.01_ikbkb-neg.tsv"
matrix_path = os.path.join(prefix, "matrix", matrix_name)
heatmap_name = "upenn_mouse_lfc_pairwise_heatmap_inf-" + inf + "_padj-0.01_ikbkb-neg.png"
heatmap_path = os.path.join(prefix, "plots/heatmaps", heatmap_name)

matrix = pd.read_table(matrix_path, sep="\t", index_col=0)
cluster = "1"
lfc_min = 4
cluster_idx = np.logical_or(matrix.columns.str.contains(f"^\d+_vs_{cluster}$"), matrix.columns.str.contains(f"^{cluster}_vs_\d+$"))
matrix = matrix.loc[:, cluster_idx]
gene_maxes = abs(matrix).max(axis=1)
matrix = matrix[gene_maxes > lfc_min]

clustermap = sns.clustermap(matrix, cmap='inferno', yticklabels=[], method='ward', figsize=(15, 9))
#clustermap = sns.clustermap(matrix, cmap='inferno', method='ward', figsize=(17, 9))
clustermap.ax_heatmap.set_xlabel(f"Pairwise LFC ({np.shape(matrix)[1]})")
clustermap.ax_heatmap.set_xticklabels(clustermap.ax_heatmap.get_xticklabels(), rotation=90)
clustermap.ax_heatmap.set_ylabel(f"Genes ({np.shape(matrix)[0]}) with max(abs(lfc)) > {lfc_min}")

plt.show()
