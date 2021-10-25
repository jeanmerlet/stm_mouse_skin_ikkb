import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.decomposition
import sklearn.manifold
import scprep
import time
import os
import re
import glob

np.random.seed(888)
random_state = 888

root_mtx_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/mesenchymal_count_mtxs/relevant_samples'
sample_ids_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/mesenchymal_count_mtxs/dermis_ids.txt'

# combining the matrices
sample_metadata = pd.read_csv(sample_ids_path, header=0, index_col=None, sep='\t')
sample_ids = sample_metadata.loc[:, 'sample_id'].values.astype(str).tolist()

sample_dirs = []
for sample_id in sample_ids:
    sample_dirs.append(os.path.join(root_mtx_dir, sample_id))

samples = []
for i, d in enumerate(sample_dirs):
    samples.append(scprep.io.load_10X(d, gene_labels='both'))

mtx, labels = scprep.utils.combine_batches(data=samples, batch_labels=sample_ids, append_to_cell_names=True)
del(samples)

# normalization stuff
mtx = scprep.filter.remove_rare_genes(mtx, cutoff=0, min_cells=5)
mtx = scprep.normalize.library_size_normalize(mtx)
mtx = scprep.transform.log(mtx)

# save matrix for use by R / Seurat for batch effect removal
mtx = mtx.transpose()
mtx = mtx.sparse.to_dense()
out_mtx_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/mesenchymal_count_mtxs/combined_human-ad_mtx.tsv'
mtx.to_csv(out_mtx_path, sep='\t', header=True, index=True)
