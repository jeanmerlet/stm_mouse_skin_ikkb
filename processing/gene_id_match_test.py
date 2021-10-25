import csv
import gzip
import os
import scipy.io
import numpy as np

matrix_dir1 = "Users/6j9/projects/mouse/data/FGC2063_5_86846/filtered_feature_bc_matrix/"
matrix_dir1 = "Users/6j9/projects/mouse/data/FGC2063_5_86847/filtered_feature_bc_matrix/"

bc_path1 = os.path.join(matrix_dir1, "barcodes.tsv")
bc_path2 = os.path.join(matrix_dir2, "barcodes.tsv")
bc1 = [row[0] for row in csv.reader(open(bc_path1, mode='rt'), delimiter="\t")]
bc2 = [row[0] for row in csv.reader(open(bc_path2, mode='rt'), delimiter="\t")]
bc1 = np.asarray(bc1)
bc2 = np.asarray(bc2)

print(len(np.intersect1d(bc1, bc2)))
