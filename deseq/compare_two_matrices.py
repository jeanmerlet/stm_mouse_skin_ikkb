from scipy.stats import mannwhitneyu
import pandas as pd
import numpy as np
import os

root = "/Users/6j9/projects/mouse/deseq/matrix"
path_a = os.path.join(root, "upenn_mouse_lfc_pairwise_matrix_inf-1.30_padj-0.01.tsv")
path_b = os.path.join(root, "upenn_mouse_lfc_pairwise_matrix_inf-1.30_padj-0.01_ikbkb-neg.tsv")

def compare_matrices(b, a):
  a, b = a.values.flatten(), b.values.flatten()
  print(mannwhitneyu(a, b, alternative='two-sided'))
  return
  u, p = mannwhitneyu(a, b)
  effect_size = u / (len(a) * len(b))
  print(f"common language effect size: {effect_size}, overlap p-statistic: {p}")

mat_a = pd.read_csv(path_a, sep="\t", index_col=0)
mat_b = pd.read_csv(path_b, sep="\t", index_col=0)

test = compare_matrices(mat_a, mat_b)
