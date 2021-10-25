import csv
import gzip
import os
import scipy.io
import numpy as np

experiments = ["/FGC2063_5_86846", "/FGC2063_5_86847", "/FGC2063_5_86848", "/FGC2063_5_86849", "/FGC2091_7_92848"]
cols = []

for i, experiment in enumerate(experiments):
  name = experiment[-5:]
  print(name)
  prefix = "./data" + experiment + "/filtered_feature_bc_matrix"
  barcodes_path = os.path.join(prefix, "barcodes.tsv.gz")
  barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode='rt'), delimiter="\t")]
  barcodes = [barcode + '_' + name for barcode in barcodes]
  bc_col = np.asarray(barcodes)
  bc_len = len(bc_col)
  exp_col = np.full(bc_len, name)
  if i < 2:
    exp_type = "control"
  else:
    exp_type = "case"
  type_col = np.full(bc_len, exp_type)
  cols.append(np.vstack((bc_col, exp_col, type_col)).T)

print(np.shape(cols))
cols = np.concatenate(cols, axis=0)
print(np.shape(cols))

filename = "full_cytoscape_cols.txt"
np.savetxt(filename, cols, delimiter=",", fmt="%s")
