import csv
import gzip
import os
import scipy.io
import numpy as np

matrix_dir = "./data/FGC2091_7_92848/filtered_feature_bc_matrix"
m = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))

barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode='rt'), delimiter="\t")]
barcodes = [barcode + '_92848' for barcode in barcodes]

features_path = os.path.join(matrix_dir, "features.tsv.gz")
#features_path1 = os.path.join(matrix_dir1, "features.tsv.gz")
ensemble_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode='rt'), delimiter="\t")]
#ensemble_ids1 = [row[0] for row in csv.reader(gzip.open(features_path1, mode='rt'), delimiter="\t")]

#for i in range(len(ensemble_ids)):
  #if ensemble_ids[i] != ensemble_ids1[i]:
  #  print(ensemble_ids[i], ensemble_ids1[i])

ensemble_ids = np.asarray(ensemble_ids)
ensemble_ids = np.insert(ensemble_ids, 0, "")



m = m.todense()
print("densed")
m = np.vstack((barcodes, m))
print("vstacked")
m = np.column_stack((ensemble_ids, m))
print("hstacked")
print(m[0:10, 0:10])
print(np.shape(m))

np.savetxt("filtered_matrix_92848_with_geneids.csv.gz", m, delimiter=",", fmt="%s")
