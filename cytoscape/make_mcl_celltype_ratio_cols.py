import csv
import numpy as np
import re
import matplotlib.pyplot as plt

mcl_dir = "/Users/6j9/projects/mouse/irf-mcl/mcl_clusters/"
mcl_file = mcl_dir + "mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.30.txt"
barcodes_dir = "/Users/6j9/projects/mouse/data/"
cell_type_dir = "/Users/6j9/projects/mouse/cytoscape/coloring/"
#out_file = cell_type_dir + "kang2_celltype-filter_cytoscape_cols_inf-1.20.txt"
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

  cell_type_path = cell_type_dir + "raw/" + name + "_cell_types_kang2.tsv"
  cell_types = [row[0:2] for row in csv.reader(open(cell_type_path, mode='rt'), delimiter="\t")]
  ct_col = np.asarray(cell_types)

  cell_type_col = np.full(bc_len, "other", dtype='U35')
  idx = np.in1d(bc_col, ct_col[:,0])
  cell_type_col[idx] = ct_col[:,1]

  exp_col = np.full(bc_len, name)
  if i < 2:
    exp_type = "control"
  else:
    exp_type = "case"
  type_col = np.full(bc_len, exp_type)

  combined = np.vstack((bc_col, cell_type_col, type_col)).T
  cols.append(combined)
cols = np.concatenate(cols, axis=0)

mcl_clusters = []
cluster_count = 0 
with open(mcl_file, 'rt') as in_file:
  for line in in_file:
    cluster_count += 1
    mcl_clusters.append(line.strip().split('\t'))

cell_names = np.unique(cols[:, 1])
#cell_names = ["B cells", "Dendritic cells", "Erythroid-like and erythroid precur", "Gamma delta T cells", "Mast cells", "NK cells", "Nuocytes", "Plasma cells", "T cells", "T regulatory cells"]
#cell_names = ["Keratinocytes"]
cell_names = ['other']
all_ratios = []
all_descs = []
cluster_list = []
for i, cluster in enumerate(mcl_clusters):
  break
  current_ratios = []
  current_descs = []
  total_cells = len(cluster)
  total_ratio = 0
  total_count = 0
  for cell_name in cell_names:
    # get barcodes matching cell type
    current_bc = cols[np.nonzero(cols[:, 1] == cell_name)][:, 0]
    # use those barcodes to count that cell type in cluster
    bc_count = len(np.intersect1d(cluster, current_bc))
    ratio = bc_count / total_cells
    ratio_desc = "cluster " + str(i + 1) + ": " + str(ratio) + " " + cell_name + ", " + str(bc_count) + " cells"
    #print(ratio_desc)
    current_ratios.append(ratio)
    current_descs.append(ratio_desc)
    total_ratio += ratio
    total_count += bc_count
  #if total_ratio >= 0.50:
  #if total_ratio >= 0.50 and total_count >= 10:
  if total_ratio <= 0.80:
    cluster_list.append(i + 1)
    #print(ratio_desc)
    print("")
  #print("")
  all_descs.append(current_descs)
  all_ratios.append(np.reshape(np.asarray(current_ratios), (1, len(current_ratios))))
#all_ratios = np.concatenate(all_ratios, axis=0)
#print(cluster_list)
#print(np.unique(cluster_list))

new_cols = np.vstack((cols[:, 0], np.full(len(cols[:, 0]), "delete"))).T
for i, ratio in enumerate(all_ratios):
  break
  if ratio[np.where(cell_names == 'Prrx1+ fibroblasts')] + ratio[np.where(cell_names == 'hypodermal fibroblasts')]  >= 0.50:
    print('\n'.join(all_descs[i]))
    print("")
    #new_cols[:, 1][np.in1d(new_cols[:, 0], mcl_clusters[i])] = "keep"

#np.savetxt(out_file, new_cols, delimiter=",", fmt="%s")

def count_cells_per_mcl_cluster(mcl_clusters):
  for i, mcl_cluster in enumerate(mcl_clusters):
    num_cells = len(mcl_cluster)
    print(f"cluster {i+1}: {num_cells} cells - {mcl_cluster[0]}")

count_cells_per_mcl_cluster(mcl_clusters)
