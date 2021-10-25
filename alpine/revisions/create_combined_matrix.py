import scipy.io
import gzip
import csv
import pandas as pd
import numpy as np
import os

data_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data'
out_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/combined_mtx.tsv'

mtx_dirs = []
for r, d, f in os.walk(data_dir):
    for exp_dir in d:
        if 'FGC' in exp_dir:
            mtx_dirs.append(os.path.join(r, exp_dir, 'filtered_feature_bc_matrix'))

mtx_dirs.sort()

exp_names = ['86846', '86847', '86848', '86849', '92848']
for i, mtx_dir in enumerate(mtx_dirs):
    mtx = scipy.io.mmread(os.path.join(mtx_dir, 'matrix.mtx.gz'))
    mtx = mtx.todense()
    barcodes_path = os.path.join(mtx_dir, "barcodes.tsv.gz")
    barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode='rt'), delimiter="\t")]
    barcodes = [bc + '_' + exp_names[i] for bc in barcodes]
    features_path = os.path.join(mtx_dir, "features.tsv.gz")
    genes = [row[1].upper() for row in csv.reader(gzip.open(features_path, mode='rt'), delimiter="\t")]
    labeled_mtx = pd.DataFrame(data=mtx, index=genes, columns=barcodes)
    if i == 0:
        combined_mtx = labeled_mtx.copy()
    else:
        combined_mtx = pd.concat([combined_mtx, labeled_mtx], axis=1)

combined_mtx.to_csv(out_path, sep='\t', index=True, header=True)
