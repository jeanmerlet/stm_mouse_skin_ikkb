import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import re

mtx_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/combined_mtx.tsv'
mtx = pd.read_csv(mtx_path, sep='\t', index_col=0, header=0)

cell_types_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/figure2v2.tsv'
cell_types = pd.read_csv(cell_types_path, index_col=0, sep='\t', header=None)
# convert to series for str matching
cell_types = cell_types.iloc[:, 0]
idx1 = cell_types.str.contains('fibroblast')
idx2 = cell_types.str.contains('Dermal')
idx = np.logical_or(idx1, idx2)

fibroblasts = cell_types.index[idx].values

genes = ['CEBPB', 'PRRX1', 'CCL11']
genes = ['CXCL12', 'CCL7', 'CCL11']
genes = ['CXCL1', 'CXCL8']
genes = ['CCL11']
gene_mtx = mtx.loc[genes, fibroblasts]

gene_mtx.loc['co_exp', :] = (gene_mtx != 0).all(axis=0)
control_fibro_idx= [True if re.search(r'8684[6,7]', value) else False for value in gene_mtx.columns.values]
control_gene_mtx = gene_mtx.iloc[:, control_fibro_idx]
exp_gene_mtx = gene_mtx.iloc[:, np.logical_not(control_fibro_idx)]

cebpb_control_count = np.count_nonzero(control_gene_mtx.loc['CEBPB', :] > 0)
cebpb_exp_count = np.count_nonzero(exp_gene_mtx.loc['CEBPB', :] > 0)

tri_control_count = np.count_nonzero(control_gene_mtx.loc['co_exp', :])
tri_exp_count = np.count_nonzero(exp_gene_mtx.loc['co_exp', :])

gene = 'PRRX1'
gene_control_count = np.count_nonzero(control_gene_mtx.loc[gene, :] > 0)
gene_exp_count = np.count_nonzero(exp_gene_mtx.loc[gene, :] > 0)

# stuff by cell type for ccl11
ccl11 = mtx.loc['CCL11', :]
bcs = ccl11.index.values
uniq_cell_types = np.unique(cell_types)

ccl11_pos_counts = np.zeros(len(uniq_cell_types), dtype=np.int32)
for i, cell_type in enumerate(uniq_cell_types):
    bcs_idx = cell_types == cell_type
    ccl11_pos_counts[i] = np.count_nonzero(ccl11[bcs_idx])

x = np.arange(len(ccl11_pos_counts))
fig, ax = plt.subplots(figsize=(12, 8))
ax.bar(x, ccl11_pos_counts)
ax.set_xticks(x)
ax.set_xticklabels(uniq_cell_types, rotation=90)
ax.set_xlabel('Cell Type')
ax.set_ylabel('Number of Cells')
ax.set_title('CCL11 Expressing Cells')

fig_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/ccl11_bar_plot.pdf'
plt.tight_layout()
plt.savefig(fig_path, format='pdf')





