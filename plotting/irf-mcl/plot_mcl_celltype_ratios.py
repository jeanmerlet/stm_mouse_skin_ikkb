import csv
import numpy as np
import re
import matplotlib.pyplot as plt

barcodes_dir = "/Users/6j9/projects/mouse/data/"
cell_type_dir = "/Users/6j9/projects/mouse/cytoscape/coloring/"
mcl_dir = "/Users/6j9/projects/mouse/irf-mcl/mcl_clusters/"
mcl_file = mcl_dir + "mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.20.txt"
experiments = ["FGC2063_5_86846", "FGC2063_5_86847", "FGC2063_5_86848", "FGC2063_5_86849", "FGC2091_7_92848"]
cols = []

for i, experiment in enumerate(experiments):
  name = experiment[-5:]
  prefix = barcodes_dir + experiment + "/filtered_feature_bc_matrix/"
  barcodes_path = prefix + "barcodes.tsv"
  barcodes = [row[0] for row in csv.reader(open(barcodes_path, mode='rt'), delimiter="\t")]
  barcodes = [barcode + '_' + name for barcode in barcodes]
  bc_col = np.asarray(barcodes)
  bc_len = len(bc_col)

  cell_type_path = cell_type_dir + name + "_cell_types_kang2.tsv"
  cell_types = [row[0:2] for row in csv.reader(open(cell_type_path, mode='rt'), delimiter="\t")]
  ct_col = np.asarray(cell_types)

  cell_type_col = np.full(bc_len, "other", dtype='U35')
  idx = np.in1d(bc_col, ct_col[:,0])
  cell_type_col[idx] = ct_col[:,1]
  combined = np.vstack((bc_col, cell_type_col)).T
  cols.append(combined)
cols = np.concatenate(cols, axis=0)

mcl_clusters = []
cluster_count = 0 
with open(mcl_file, 'rt') as in_file:
  for line in in_file:
    cluster_count += 1
    mcl_clusters.append(line.strip().split('\t'))

cell_names = np.unique(cols[:, 1])
all_ratios = []
for i, cluster in enumerate(mcl_clusters):
  #if i == 50:
    #break
  current_ratios = []
  total_cells = len(cluster)
  for cell_name in cell_names:
    # get barcodes matching cell type
    current_bc = cols[np.nonzero(cols[:, 1] == cell_name)][:, 0]
    # use those barcodes to count that cell type in cluster
    bc_count = len(np.intersect1d(cluster, current_bc))
    ratio = bc_count / total_cells
    ratio_desc = "cluster " + str(i + 1) + ": " + str(ratio) + " " + cell_name + " cells"
    #print(ratio_desc)
    current_ratios.append(ratio)
  #print("")
  all_ratios.append(np.reshape(np.asarray(current_ratios), (1, len(current_ratios))))
all_ratios = np.concatenate(all_ratios, axis=0)

inflation = re.search('(?<=clusters_)\d.\d+', mcl_file).group(0)
title = 'Cumulative cell type ratios - ' + inflation + ' inf'
#colors = ['green', 'white', 'purple', 'white', 'salmon', 'yellow', 'white', 'teal', 'black', 'blue', 'white', 'pink', 'white', 'red', 'grey']
colors = ['blue', 'black', 'pink', 'grey', 'red']

fig, ax = plt.subplots(figsize=(8, 5))
ax.set_title(title)
ax.axvline(0.80, label='cutoff', color='green')
ax.hist(all_ratios, cluster_count, histtype='step', label=cell_names, cumulative=True, density=True, color=colors)
ax.set_xlim(0.0, 1.0)
ax.set_xticks(np.arange(0, 1.1, 0.1))
ax.set_ylim(1.0, 0.0)
ax.set_yticks(np.arange(0.0, 1.1, 0.1))
ax.set_xlabel("ratio")
ax.set_ylabel("likelihood of occurence")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper right')
ax.grid(True)

fig.tight_layout()
plt.show()
