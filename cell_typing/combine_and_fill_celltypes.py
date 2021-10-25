import csv
import gzip
import scipy.io
import numpy as np
import re
import os

barcodes_dir = "/Users/6j9/projects/mouse/data/"
cell_type_dir = "/Users/6j9/projects/mouse/annotation/cellassign/celltypes/jean4"
experiments = ["FGC2063_5_86846", "FGC2063_5_86847", "FGC2063_5_86848", "FGC2063_5_86849", "FGC2091_7_92848"]
cols = []

for i, experiment in enumerate(experiments):
  name = experiment[-5:]
  prefix = os.path.join(barcodes_dir, experiment, "filtered_feature_bc_matrix")
  barcodes_path = os.path.join(prefix, "barcodes.tsv")
  barcodes = [row[0] for row in csv.reader(open(barcodes_path, mode='rt'), delimiter="\t")]
  barcodes = [barcode + '_' + name for barcode in barcodes]
  bc_col = np.asarray(barcodes)
  bc_len = len(bc_col)

  cell_type_path = os.path.join(cell_type_dir, name + "_cell_types_jean4.tsv")
  cell_types = [row[0:2] for row in csv.reader(open(cell_type_path, mode='rt'), delimiter="\t")]
  ct_col = np.asarray(cell_types)

  cell_type_col = np.full(bc_len, "no markers", dtype='U35')
  idx = np.in1d(bc_col, ct_col[:,0])
  cell_type_col[idx] = ct_col[:,1]
  combined = np.vstack((bc_col, cell_type_col)).T
  cols.append(combined)

cols = np.concatenate(cols, axis=0)
out_dir = "/Users/6j9/projects/mouse/cytoscape/coloring/"
filename = os.path.join(out_dir, "celltype_cytoscape_cols_jean4.tsv")
np.savetxt(filename, cols, fmt='%s', delimiter='\t')
