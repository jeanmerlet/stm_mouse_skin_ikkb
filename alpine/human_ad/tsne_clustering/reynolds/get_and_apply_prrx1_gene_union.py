import pandas as pd
import numpy as np
import os

genes_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/prrx1_matrices'
gene_paths = [os.path.join(genes_dir, path) for path in os.listdir(genes_dir) if 'genes' in path]

for i, gene_path in enumerate(gene_paths):
    if i == 0:
        genes = pd.read_csv(gene_path)
    else:
        genes = np.union1d(genes, pd.read_csv(gene_path))

print(len(genes))
print(genes)

mtx_paths = [os.path.join(genes_dir, path) for path in os.listdir(genes_dir) if 'unprocessed' in path]

for i, mtx_path in enumerate(mtx_paths):
    print(i)
    mtx = pd.read_csv(mtx_path, sep='\t', header=0, index_col=0)
    mtx = mtx.transpose()
    mtx = mtx.reindex(genes, fill_value=0)
    _, tail = os.path.split(mtx_path)
    out_path = os.path.join(genes_dir, 'gene_union_' + tail)
    mtx.to_csv(out_path, sep='\t', header=True, index=True)
