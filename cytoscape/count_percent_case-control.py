import csv
import numpy as np
import re
import matplotlib.pyplot as plt

mcl_dir = "/Users/6j9/projects/mouse/irf-mcl/mcl_clusters/"
mcl_file = mcl_dir + "mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.20.txt"
barcodes_dir = "/Users/6j9/projects/mouse/data/"
cell_type_dir = "/Users/6j9/projects/mouse/cytoscape/coloring/"
out_file = cell_type_dir + "kang2_celltype-filter_cytoscape_cols_inf-1.20.txt"
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

  cell_type_path = cell_type_dir + "raw/" + name + "_cell_types2.tsv"
  #cell_type_path = cell_type_dir + "raw/" + name + "_cell_types_kang2.tsv"
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

cell_names = np.unique(cols[:, 1])
#print(cell_names)
#cell_names = ["NK cells"]
#cell_names = ["B cells", "Basophils", "Dendritic cells", "Eosinophils", "Gamma delta T cells", "Mast cells", "Monocytes", "Neutrophils", "NK cells", "Nuocytes", "Plasma cells", "Plasmacytoid dendritic cells", "T cells", "T regulatory cells", "Macrophages"]
for cell_name in cell_names:
  break
  cols[:, 1][np.nonzero(cols[:, 1] == cell_name), ] = "immune"
#cell_names = ["immune"]
#cell_names = ["Keratinocytes"]

total_count = np.shape(cols)[0]
case = np.nonzero(cols[:, 2] == "case")
case_count = np.shape(case)[1]
control = np.nonzero(cols[:, 2] == "control")
control_count = np.shape(control)[1]
case_mult = case_count / total_count
control_mult = control_count / total_count
#f = open("/Users/6j9/projects/mouse/cytoscape/cols1.txt", 'ab')
for cell_name in cell_names:
  break
  cells = np.nonzero(cols[:, 1] == cell_name)
  #np.savetxt(f, cols[:, 0][cells], fmt='%s')
  cells_count = np.shape(cells)[1]
  cells_control = cols[np.intersect1d(cells, control)]
  cells_control_count = np.shape(cells_control)[0]
  cells_case_count = cells_count - cells_control_count

  scaled_control = cells_control_count / control_mult
  scaled_case = cells_case_count / case_mult
  scaled_total = scaled_control + scaled_case
  if scaled_total > 0:
    print(f"{cell_name}: {cells_count}, cells in control: {cells_control_count}")
    print(f"% case: {scaled_case / scaled_total}")
    print(f"% control: {scaled_control / scaled_total}")
    print(f"% total: {cells_count / total_count}")
    print("")
#f.close()

def count_cells_per_celltype(cols, cell_names):
  for cell_name in cell_names:
    idx = np.nonzero(cols[:, 1] == cell_name)
    num_cells = np.shape(cols[idx])[0]
    print(f"{cell_name}: {num_cells}")

#count_cells_per_celltype(cols, cell_names)
