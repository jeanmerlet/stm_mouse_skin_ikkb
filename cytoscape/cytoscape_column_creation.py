import csv
import os
import scipy.io
import numpy as np

data_dir = "/Users/6j9/projects/mouse/data/"
mcl_file = "/Users/6j9/projects/mouse/phate/mcl/mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.25.txt"

experiments = ['FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848']
experiment_names = ['86846', '86847', '86848', '86849', '92848']

for experiment in experiments:
  working_dir = data_dir + experiment + "/filtered_feature_bc_matrix/"
  exp_name = experiment[-5:]
  print(exp_name)
  barcodes_path = working_dir + "barcodes.tsv"
  barcodes = [row[0] for row in csv.reader(open(barcodes_path, mode='rt'), delimiter="\t")]
  barcodes = [barcode + '_' + exp_name for barcode in barcodes]
  bc_col = np.asarray(barcodes)
  exp_col = np.full(bc_len, exp_name)
  type_col = np.full(bc_len, "control")
  cols = np.vstack((bc_col, exp_col, type_col)).T
  filename = exp_name + "_cytoscape_cols.txt"
  np.savetxt(prefix + '/' + filename, cols, delimiter=",", fmt="%s")
  
