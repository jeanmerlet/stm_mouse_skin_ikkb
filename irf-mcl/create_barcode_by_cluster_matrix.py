import numpy as np
import pandas as pd

data_dir = "/Users/6j9/projects/mouse/data/"
mcl_file = "/Users/6j9/projects/mouse/irf-mcl/mcl_clusters/mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.30.txt"
bc_cluster_out_filepath = "/Users/6j9/projects/mouse/irf-mcl/barcodes_by_cluster/mouse_barcodes_by_clusters_1.30inf.tsv"

experiments = ["FGC2063_5_86846", "FGC2063_5_86847", "FGC2063_5_86848", "FGC2063_5_86849", "FGC2091_7_92848"]
bc = []
for experiment in experiments:
  name = "_" + experiment[-5:]
  bc_path = data_dir + experiment + "/filtered_feature_bc_matrix/barcodes.tsv"
  barcodes = pd.read_csv(bc_path, sep="\t", header=None)
  barcodes = [s + name for s in np.squeeze(barcodes.values)]
  bc.append(barcodes)

bc = [barcode for bc_set in bc for barcode in bc_set]
num_clusters = 316
clusters = np.arange(1, num_clusters + 1)
mat = np.zeros((len(bc), num_clusters))
print(np.shape(mat))
out_mat = pd.DataFrame(data=mat, index=bc, columns=clusters)

cluster_count = 0 
with open(mcl_file, "rt") as in_file:
  for line in in_file:
    cluster_count += 1
    barcode_list = line.strip().split("\t")
    for barcode in barcode_list:
      out_mat.loc[barcode, cluster_count] = 1

out_mat = out_mat.astype(int)
out_mat.to_csv(bc_cluster_out_filepath, sep="\t")
