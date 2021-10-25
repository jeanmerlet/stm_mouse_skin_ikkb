import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

ct_path = '/Users/6j9/projects/mouse/data/cell_typing/cellassign_results/jean-kang_final/jean-kang2_cts_all.tsv'
mcl_path = '/Users/6j9/projects/mouse/data/irf-mcl/mcl_clusters/mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.30.txt'

celltypes = pd.read_csv(ct_path, sep='\t', names=['cell', 'type'], index_col=None)
celltypes.index = celltypes['cell']

#print(celltypes.shape)
#print(len([bc for bc in celltypes['cell'] if bc[-5:] in ['86846', '86847']]))

immune = [3, 5, 6, 9, 12, 14, 16, 17, 18, 19, 24, 26, 28, 30, 35, 36, 37, 40, 41]
cluster_num = 0
cluster_cts = {}
with open(mcl_path) as clusters:
  for line in clusters:
    cluster_num += 1
    if cluster_num in immune:
      bcs = line.strip('\n').split('\t')
      cluster_cts[cluster_num] = celltypes.reindex(index=bcs)

# cre+ vs. cre- cts per cluster
cre_cluster_cts = {}
num_cells = []
for cluster, cts in cluster_cts.items():
  grouped = cts.groupby('type')['cell'].apply(list).to_frame()
  cre_cts = {}
  total = 0
  for ct in grouped.index:
    #print(len(grouped.loc[ct, 'cell']))
    if len(grouped.loc[ct, 'cell']) > 5:
      cre = []
      #cre.append(int(len([bc for bc in grouped.loc[ct, 'cell'] if bc[-5:] in ['86846', '86847']])*.836))
      #cre.append(int(len([bc for bc in grouped.loc[ct, 'cell'] if bc[-5:] not in ['86846', '86847']])*1.25))
      cre.append(int(len([bc for bc in grouped.loc[ct, 'cell'] if bc[-5:] in ['86846', '86847']])))
      cre.append(int(len([bc for bc in grouped.loc[ct, 'cell'] if bc[-5:] not in ['86846', '86847']])*1.49))
      cre_cts[ct] = cre
      total += np.sum(cre)
  cre_cluster_cts[cluster] = cre_cts
  num_cells.append(total)

i = 0
for cluster, cts in cre_cluster_cts.items():
  break
  print(cluster)
  total_cells = num_cells[i]
  fig, ax = plt.subplots(figsize=(8, 5))
  cre_pos = ax.barh(np.arange(len(cts))-0.25, [x[0]/total_cells for x in cts.values()], height=0.4)
  cre_neg = ax.barh(np.arange(len(cts))+0.25, [x[1]/total_cells for x in cts.values()], height=0.4)
  ax.set_yticks(range(len(cts)))
  ax.set_yticklabels(cts.keys())
  ax.legend((cre_pos[0], cre_neg[0]), ('cre-', 'cre+'))
  ax.set_xticks(np.arange(0, 1.1, 0.1))
  ax.set_title('Cluster ' + str(cluster))
  ax.set_xlabel('Proportion of Cell Type by Condition')
  ax.set_ylabel('Cell Type')
  plt.tight_layout()
  out_path = os.path.join('/Users/6j9/projects/mouse/plots/final/fig3/ct_dists/jean-kang1_markers', str(cluster) + '.pdf')
  plt.savefig(out_path)
  #plt.show()
  plt.clf()
  i += 1

# cre+ vs. cre- combined from above clusters
cre_total_cts = {}
total_cells = np.sum(num_cells)
for cluster, cts in cre_cluster_cts.items():
  for ct, count in cts.items():
    if ct not in cre_total_cts.keys():
      cre_total_cts[ct] = [0, 0]
    cre_total_cts[ct][0] += count[0]
    cre_total_cts[ct][1] += count[1]

fig, ax = plt.subplots(figsize=(8, 5))
cre_pos = ax.barh(np.arange(len(cre_total_cts))-0.25, [x[0]/total_cells for x in cre_total_cts.values()], height=0.4)
cre_neg = ax.barh(np.arange(len(cre_total_cts))+0.25, [x[1]/total_cells for x in cre_total_cts.values()], height=0.4)
ax.set_yticks(range(len(cre_total_cts)))
ax.set_yticklabels(cre_total_cts.keys())
ax.legend((cre_pos[0], cre_neg[0]), ('cre-', 'cre+'))
ax.set_xticks(np.arange(0, 1.1, 0.1))
ax.set_title('All Immune Clusters')
ax.set_xlabel('Proportion of Cell Type by Condition')
ax.set_ylabel('Cell Type')
plt.tight_layout()
out_path = os.path.join('/Users/6j9/projects/mouse/plots/final/fig3/ct_dists/jean-kang1_markers', 'all.pdf')
plt.savefig(out_path)
#plt.show()
