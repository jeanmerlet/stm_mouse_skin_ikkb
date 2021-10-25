import csv
import gzip
import os
import scipy.io
import numpy as np

matrix_dir1 = "./data/filtered_feature_bc_matrix/FGC2091_7_92847"
matrix_dir2 = "./data/filtered_feature_bc_matrix/FGC2091_7_92848"

features_path1 = os.path.join(matrix_dir1, "features.tsv.gz")
ensemble_ids1 = [row[0] for row in csv.reader(gzip.open(features_path1, mode='rt'), delimiter="\t")]
ensemble_ids1 = np.asarray(ensemble_ids1)

features_path2 = os.path.join(matrix_dir2, "features.tsv.gz")
ensemble_ids2 = [row[0] for row in csv.reader(gzip.open(features_path2, mode='rt'), delimiter="\t")]
ensemble_ids2 = np.asarray(ensemble_ids2)

for i, (id1, id2) in enumerate(zip(ensemble_ids1, ensemble_ids2)):
  if id1 != id2:
    print(f"{id1} != {id2}")
    break
