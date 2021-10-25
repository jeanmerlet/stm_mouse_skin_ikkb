import pandas as pd
import numpy as np
import re

mcl_file = '/Users/6j9/projects/mouse/data/irf-mcl/mcl_clusters/mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.30.txt'
counts_out_path = '/Users/6j9/projects/mouse/data/irf-mcl/counts/inf-1.30_cre_counts.tsv'
ratio = 0.9
cat_out_path = '/Users/6j9/projects/mouse/data/irf-mcl/counts/inf-1.30_cre_categories_' + str(ratio) + '_ratio.tsv'

cluster_count = 0
mcl_clusters = {}
with open(mcl_file, 'rt') as in_file:
  for line in in_file:
    cluster_count += 1
    line_list = line.strip().split('\t')
    for item in line_list:
      mcl_clusters[cluster_count] = line_list

# create df with cre pos and neg counts
cre_neg = ['86846', '86847']
index = np.arange(1, len(mcl_clusters.keys()) + 1)
cre_counts = pd.DataFrame(data=np.zeros((len(index), 2), dtype=int), index=index, columns=['cre_neg', 'cre_pos'])
for cluster, bcs in mcl_clusters.items():
  total_cells = len(bcs)
  cre_counts.loc[cluster, 'cre_neg'] = sum(acc_id in s for acc_id in cre_neg for s in bcs)
  cre_counts.loc[cluster, 'cre_pos'] = total_cells - cre_counts.loc[cluster, 'cre_neg']

cre_counts.to_csv(counts_out_path, sep='\t', header=True, index=True)

# separate by ratio into cre-, cre+, mixed
cre_category = pd.DataFrame(data=np.full((len(index), 1), 'mixed'), index=index, columns=['category'])
for cluster in index:
  total_bcs = np.sum(cre_counts.loc[cluster, :])
  if cre_counts.loc[cluster, 'cre_pos'] / total_bcs >= ratio:
    cre_category.loc[cluster, 'category'] = 'cre_pos'
  elif cre_counts.loc[cluster, 'cre_neg'] / total_bcs >= ratio:
    cre_category.loc[cluster, 'category'] = 'cre_neg'

cre_category.to_csv(cat_out_path, sep='\t', header=True, index=True)
