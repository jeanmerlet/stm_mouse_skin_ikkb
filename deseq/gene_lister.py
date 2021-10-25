import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os
import re
import sys

sys.setrecursionlimit(10**6)

# make text in pdf export editable
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

inf = '1.30'
prefix = '/Users/6j9/projects/mouse/data/deseq2'
matrix_dir = os.path.join(prefix, "matrix")

# use for fibro (ikbkb deletion)
matrix_name = 'upenn_mouse_lfc_pairwise_matrix_inf-' + inf + '_padj-0.01_ikbkb-neg.tsv'
# use for immune (don't want to ikbkb delete for immune)
#matrix_name = "upenn_mouse_lfc_pairwise_matrix_inf-" + inf + "_padj-0.01.tsv"

matrix_path = os.path.join(prefix, "matrix", matrix_name)
out_dir = os.path.join(prefix, "subcomparisons")

if inf == "1.20":
  fibro_cre_pos = [5, 7, 9, 15, 17, 25, 27]
  fibro_cre_neg = [2, 3, 6, 10, 14, 32, 33]
  fibro_prrx1_cre_pos = [5, 7, 17, 25, 27]
  fibro_prrx1_cre_neg = [2, 3, 6, 10, 14, 32]
  immune_cre_pos = [1, 11, 18, 28]
  immune_cre_neg = [4, 20, 30]
  immune = [1, 4, 8, 11, 18, 20, 21, 23, 28, 30]
elif inf == "1.25":
  fibro_cre_pos = [7, 8, 11, 17, 22, 23, 28, 35, 38, 44, 51]
  fibro_cre_neg = [2, 3, 4, 9, 10, 25, 26, 27, 41, 45, 46, 52, 53, 54]
  fibro_prrx1_cre_pos = [7, 8, 11, 17, 22, 35, 38, 51]
  fibro_prrx1_cre_neg = [2, 3, 4, 10, 25, 27, 45, 46, 53, 54]
  fibro_prrx1_neg_cre_pos = [23, 28, 44]
  fibro_prrx1_neg_cre_neg = [9, 26, 41, 52]
  immune_cre_pos = [1, 5, 18, 20, 24]
  immune_cre_neg = [6, 12, 14, 16, 32, 33, 34, 39, 55]
elif inf == "1.30":
  fibro_cre_pos = [1, 11, 23, 33, 38, 39, 44, 45]
  fibro_cre_neg = [2, 4, 7, 8, 13, 15, 20, 25, 31, 34, 47]
  fibro_prrx1_cre_pos = [1, 11, 33, 38, 45]
  fibro_prrx1_cre_neg = [4, 7, 13, 15, 31, 34, 47]
  fibro_prrx1_neg_cre_pos = [23, 39, 44]
  fibro_prrx1_neg_cre_neg = [2, 8, 20, 25]
  immune_cre_neg_macro = [5, 16, 37, 41, 14, 17]
  immune_cre_pos_macro = [30]
  immune_cre_neg = [5, 16, 37, 41, 14, 17, 35, 18, 36, 40, 28]
  immune_cre_pos = [30, 19, 24, 3, 26, 6, 12, 9]

#references = [1]
#levels = fibro_cre_pos

#references = fibro_cre_neg
#levels = fibro_cre_pos

references = fibro_prrx1_cre_neg
levels = fibro_prrx1_cre_pos

matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)

def subset_comparisons(matrix, levels, references):
  num_comparisons = int(len(levels) * len(references))
  new_matrix = pd.DataFrame(data=np.zeros((np.shape(matrix)[0], num_comparisons)), index=matrix.index)
  i = 0
  for column in matrix.columns:
    cluster_numbers = re.match("(\d+)_vs_(\d+)", column).groups()
    lev, ref = cluster_numbers[0], cluster_numbers[1]
    if int(ref) in references and int(lev) in levels:
      new_matrix = new_matrix.rename(columns={i: column})
      new_matrix[column] = matrix[column]
      i += 1
    elif int(lev) in references and int(ref) in levels:
      colname = str(ref) + "_vs_" + str(lev)
      new_matrix = new_matrix.rename(columns={i: colname})
      new_matrix[colname] = -matrix[column]
      i += 1
  return new_matrix

def find_matching_genes(matrix, cutoff, drop_any_na=False, na_to_zero=True, drop_zero_rows=False, abs_sort=True, no_zeros=False):
  matrix = matrix[abs(matrix) >= cutoff]
  if drop_any_na: matrix = matrix[matrix.isna().sum(axis=1) == 0]
  if na_to_zero: matrix = matrix.fillna(0)
  if drop_zero_rows: matrix = matrix[matrix.sum(axis=1) != 0]
  if abs_sort: matrix['sums'] = abs(matrix).sum(axis=1)
  else: matrix['sums'] = matrix.sum(axis=1)
  matrix = matrix.sort_values(by=['sums'], ascending=False)
  matrix = matrix.drop('sums', axis=1)
  matrix = matrix[matrix.notna().sum(axis=1) > 0]
  if no_zeros:
    matrix = matrix[~(matrix == 0).any(axis=1)]
  return matrix

def combine_pairwise_comparisons(matrix, references, levels):
  first = True
  for ref in references:
    refs = [ref]
    levels = [lev for lev in levels if lev is not ref]
    comparisons = subset_comparisons(matrix, levels, refs)
    if first == True:
      combined = comparisons
      first = False
    else:
      combined = pd.concat([combined, comparisons], axis=1, sort=False)
  return(combined)


combined = combine_pairwise_comparisons(matrix, references, levels)
cutoff = 0.01

# save here #
matching = find_matching_genes(combined, cutoff, drop_zero_rows=True, abs_sort=False)
top_n = 500
#matching = pd.concat([matching.iloc[:top_n, :], matching.iloc[-top_n:, :]])
matching = matching.iloc[-top_n:, :]
matching = matching.iloc[::-1]
out_name = "cre-pos_vs_cluster-1_top-" + str(top_n) + "-upreg-genes_inf-" + inf + "_ikbkb-neg.tsv"
#out_path = os.path.join(out_dir, out_name)
#matching.to_csv(out_path, sep="\t")

# plot here #
#matching = find_matching_genes(combined, cutoff, drop_zero_rows=True, abs_sort=False, no_zeros=True)
matching = find_matching_genes(combined, cutoff, drop_zero_rows=True, abs_sort=False)

# genes starting with #
#idx_a = matching.index.str.startswith('CCL')
#idx_b = matching.index.str.startswith('CXCL')
#idx = [any(tup) for tup in zip(idx_a, idx_b)]
#print(idx)
#idx_c = matching.index.str.startswith('IL')
#idx = [any(tup) for tup in zip(idx_c)]
starts_list = ['CCL', 'CXCL']
all_idx = []
for prefix in starts_list:
  break
  idx = matching.index.str.startswith(prefix)
  all_idx = [any(tup) for tup in zip(all_idx, idx)]
#matching = matching.loc[all_idx]

# if doing a startswith heatmap use this block
#starts_list = ['COL', 'CXCL', 'IL', 'CD']
#starts_list = ['COL', 'CXCL', 'IL']
#starts_list = ['COL', 'CXCL']
#starts_list = ['COL', 'CCL']
#starts_list = ['CCL']

# if doing a set of specific genes use one of these blocks#
# fibro list
#gene_list = ['IL6', 'CCL11', 'CCL7', 'IL6RA', 'IL1R2', 'IL33', 'IL1R1', 'IL6ST', 'IL4RA', 'IL17RA', 'CXCL1', 'CXCL12', 'CCL19']
#gene_list = ['COL13A1', 'COL14A1']
#gene_list = ['DCN', 'GPX3', 'SPARC', 'PLAC8']
#gene_list = ['CD207', 'MFGE8', 'H2-M2', 'PXDC1', 'TBC1D4', 'LTC4S', 'GRASP', 'COL16A1', 'TNFAIP2']
#gene_list = ['CCL7', 'CCL11', 'CXCL1', 'CCL19']



# immune list
#gene_list = ['TGFB1', 'IL1B', 'IL6', 'IL18', 'CCL7', 'TSLP']
# test
# based on 1.30 immune ikbkb-neg results from above list
#gene_list = ['AREG', 'IL1B', 'CXCL1', 'CCL7']

# hypodermal fibro cluster for 1.30
#gene_list = ['CCL11']

#for gene in gene_list:
  #if gene in matching.index: print(f"{gene} in index")

# idx vs. top n #
#idx = np.in1d(matching.index, gene_list)
#matching = matching.loc[idx]
#matching = pd.concat([matching.iloc[:25, :], matching.iloc[-25:, :]])
#matching = matching.iloc[-25:, :]
#print(np.shape(matching))
#out_name = "cre-neg-vs-cre-pos_top-500" + "_inf-" + inf + ".tsv"
#out_path = os.path.join(out_dir, out_name)
#matching = matching.iloc[-500:, :]
#matching.to_csv(out_path, sep="\t")

#title = "Immune " + inf + " M2 vs. M1 " + str(single)
#title = 'Immune ' + inf +  ' cre- vs. cre+'
#title = "Cre- vs. cre+ at 1.25 inf: CCL- and CXCL- genes"

#title = "1.30 Inflation " + inf + " cre- vs. cluster 1"
def make_heatmap(data, hide_trees=False, title=None, cbar_pos=(0.02, 0.8, 0.05, 0.15), figsize=(15, 5),
                 sort_genes=True, sort_clusters=True):
  cmap = plt.get_cmap("RdBu_r")
  clustermap = sns.clustermap(data, method='ward', xticklabels=data.columns, yticklabels=data.index, figsize=figsize,
                              cmap=cmap, center=0, row_cluster=sort_genes, col_cluster=sort_clusters)
  hm = clustermap.ax_heatmap
  if hide_trees:
    clustermap.ax_row_dendrogram.set_visible(False)
    clustermap.ax_col_dendrogram.set_visible(False)
  if title:
    if hide_trees: hm.set_title(title)
    else: clustermap.fig.suptitle(title)
  hm.set_xlabel(f"LFC2 Comparisons ({np.shape(data)[1]})")
  hm.set_xticks(np.arange(np.shape(data)[1]) + 0.5)
  hm.set_xticklabels(hm.get_xticklabels(), rotation=90, size=8)
  hm.set_ylabel(f"Genes ({np.shape(data)[0]})")
  hm.set_yticks(np.arange(np.shape(data)[0]) + 0.5)
  hm.set_yticklabels(hm.get_yticklabels(), rotation=0)
  hm.set_xticks([], [])
  hm.figure.subplots_adjust(bottom=0.20, right=0.9)
  clustermap.cax.set_position(cbar_pos)
  #hm.figure.subplots_adjust(right=0.9)
  #hm.figure.subplots_adjust(bottom=0.18, right=0.9)
  plt.savefig('/Users/6j9/projects/mouse/plots/final/fig2/de_cre-neg-vs-cre-pos_fibro_all.pdf',
              transparent=True)
  plt.show()

### FIGURE 2 GENE LIST ###

#all
#idx_a = matching.index.str.startswith('CCL')
#idx_b = matching.index.str.startswith('CX')
#idx = [any(tup) for tup in zip(idx_a, idx_b)]
idx = matching.index.str.startswith('IL')

# manual
#gene_list = ['CCL7', 'CXCL12', 'CCL8', 'CCL19', 'CCL11', 'CXCL2', 'CXCL14', 'CCL2', 'CXCL10', 'CXCL1']
#gene_list = ['CCL7', 'CXCL12', 'CCL11']
#gene_list = ['CXLCL1', 'CXCL10', 'CCL7', 'CXCL12', 'CCL11', 'CCL19', 'CCL8']
#gene_list = ['CCL7', 'CCL11', 'CCL19', 'CXCL1', 'CXCL10', 'CXCL12', 'TSLP', 'IL1A', 'IL1B', 'IL1RA', 'IL1B', 'TNFA', 'IFNG',
#             'IL17A', 'IL18', 'IL1RA', 'IL4', 'IL5', 'IL6', 'IL10', 'IL11', 'EBI3', 'TGFB1', 'AREG', 'TSLP']

# manual - downregulated cytokines / chemokines
#gene_list = ['CXCL1', 'CXCL10']

#idx = np.in1d(matching.index, gene_list)

matching = matching.loc[idx]

title = 'DE of Chemokines between PRRX1+ Cre- and Cre+ Fibroblast Clusters'
figsize = (10, 4)
#print(matching)
make_heatmap(matching, title=title, hide_trees=True, cbar_pos=(0.12, .2, .04, .62), sort_genes=False, sort_clusters=False, figsize=figsize)
