import csv
import gzip
import os
import scipy.io
import numpy as np
import heapq
import re

experiment = 'FGC2091_7_92848'
matrix_dir = './data/' + experiment + '/raw_feature_bc_matrix'
m = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))
#matrix_dir = './data/' + experiment + '/star_correct/Solo.out/GeneFull/raw'
#m = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx"))
#matrix_dir = './data/mouse_ref_genome/test5/Solo.out/GeneFull/raw'
#m = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx"))

#barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
#barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode='rt'), delimiter="\t")]
#barcodes = np.asarray(barcodes)

#genes_path = os.path.join(matrix_dir, "features.tsv.gz")
#gene_ids = [row[0] for row in csv.reader(gzip.open(genes_path, mode='rt'), delimiter="\t")]
#gene_ids = np.asarray(gene_ids)

# looking at cell ranger umi cutoff for raw matrix
num_expected = 5000
nnz_cols = m.getnnz(axis=0)
pos_umi_count = nnz_cols[nnz_cols > 0]
print(len(pos_umi_count))
top_slice = np.sort(nnz_cols)[-num_expected:]
percentile_idx = np.round(0.99 * num_expected).astype('intp')
#print(percentile_idx)
umi_cutoff = top_slice[percentile_idx]
#print(umi_cutoff)
umi_cutoff = umi_cutoff // 10
#print(umi_cutoff)
print("")
print(len(nnz_cols[nnz_cols > umi_cutoff]))
print("")

# background RNA profile creation
ambient_cutoff = 100
ambient_bcs = pos_umi_count[pos_umi_count <= 100]
#print(len(ambient_bcs))


# max 10 genes expressed by one cell:
#largest_hund = heapq.nlargest(100, nnz_cols)

def cells_lesseq_genes(x):
  return len(nnz_cols[nnz_cols <= x])

def cells_greateq_genes(x):
  return len(nnz_cols[nnz_cols >= largest_hund[x] - x])
