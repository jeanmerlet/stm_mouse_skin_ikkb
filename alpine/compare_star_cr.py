import csv
import os
import scipy.io
import numpy as np

experiment = 'FGC2091_7_92848'
matrix_dir1 = './data/' + experiment + '/raw_feature_bc_matrix'
matrix_dir2 = './data/' + experiment + '/star_cr3/Solo.out/Gene/raw'
m1 = scipy.io.mmread(os.path.join(matrix_dir1, "matrix.mtx"))
m2 = scipy.io.mmread(os.path.join(matrix_dir2, "matrix.mtx"))

barcodes_path1 = os.path.join(matrix_dir1, "barcodes.tsv")
barcodes_path2 = os.path.join(matrix_dir1, "barcodes.tsv")
barcodes1 = [row[0] for row in csv.reader(open(barcodes_path1, mode='rt'), delimiter="\t")]
barcodes2 = [row[0] for row in csv.reader(open(barcodes_path2, mode='rt'), delimiter="\t")]
barcodes1 = np.asarray(barcodes1)
barcodes2 = np.asarray(barcodes2)

genes_path1 = os.path.join(matrix_dir1, "features.tsv")
genes_path2 = os.path.join(matrix_dir2, "features.tsv")
gene_ids1 = [row[0] for row in csv.reader(open(genes_path1, mode='rt'), delimiter="\t")]
gene_ids2 = [row[0] for row in csv.reader(open(genes_path2, mode='rt'), delimiter="\t")]
gene_ids1 = np.asarray(gene_ids1)
gene_ids2 = np.asarray(gene_ids2)

nnz_col_idx1 = np.nonzero(m1.getnnz(axis=0))
nnz_col_idx2 = np.nonzero(m2.getnnz(axis=0))
nnz_bc1 = barcodes1[nnz_col_idx1]
nnz_bc2 = barcodes2[nnz_col_idx2]

differences = np.isin(nnz_bc2, nnz_bc1)
print(len(nnz_bc1))
print(len(nnz_bc2))
count = np.count_nonzero(differences)
percent = count / len(nnz_bc2)
print(f"BCs in cellranger alignment: {len(nnz_bc1)}")
print(f"BCs in star alignment: {len(nnz_bc2)}")
print(f"% of star BCs in cellranger: {percent}")
print("")

def comment():
  nnz_cols1 = m1.getnnz(axis=0)
  nnz_cols2 = m2.getnnz(axis=0)
  pos_umi_counts1 = nnz_cols1[nnz_cols1 > 0]
  pos_umi_counts2 = nnz_cols2[nnz_cols2 > 0]
  top_hund1 = nnz_cols1[nnz_cols1 > 100]
  top_hund2 = nnz_cols2[nnz_cols2 > 100]
  print(len(pos_umi_counts1))
  print(len(top_hund1))
  print("")
  print(len(pos_umi_counts2))
  print(len(top_hund2))
  print("")
