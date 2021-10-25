import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.decomposition
import sklearn.manifold
import scprep
import time
import os
import re

np.random.seed(888)
random_state = 888

mtx_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/combined_he-preprocessed_human-ad_mtx.tsv'

mtx = pd.read_csv(mtx_path, sep='\t', index_col=0, header=0)
mtx = mtx.fillna(0)
mtx = mtx.transpose()

# subsetting to prrx1+ cells
prrx1_mtx_idx = mtx.loc[:, 'PRRX1'] > 0
mtx = mtx.loc[prrx1_mtx_idx, :]

bc_to_condition_map = {'SRR11396159': 'LS',
                       'SRR11396160': 'LS',
                       'SRR11396161': 'NL',
                       'SRR11396162': 'H',
                       'SRR11396163': 'LS',
                       'SRR11396164': 'H',
                       'SRR11396165': 'LS',
                       'SRR11396166': 'H',
                       'SRR11396167': 'H',
                       'SRR11396168': 'H',
                       'SRR11396169': 'NL',
                       'SRR11396170': 'H',
                       'SRR11396171': 'H',
                       'SRR11396172': 'NL',
                       'SRR11396173': 'NL',
                       'SRR11396174': 'NL',
                       'SRR11396175': 'H'}

bc_to_sample_map = {'S1': 'SRR11396159',
                    'S2': 'SRR11396160',
                    'S3': 'SRR11396161',
                    'S4': 'SRR11396162',
                    'S5': 'SRR11396163',
                    'S6': 'SRR11396164',
                    'S7': 'SRR11396165',
                    'S8': 'SRR11396166',
                    'S9': 'SRR11396167',
                    'S10': 'SRR11396168',
                    'S11': 'SRR11396169',
                    'S12': 'SRR11396170',
                    'S13': 'SRR11396171',
                    'S14': 'SRR11396172',
                    'S15': 'SRR11396173',
                    'S16': 'SRR11396174',
                    'S17': 'SRR11396175'}

# do for each comparison
comparisons = [['H', 'LS'], ['NL', 'LS']]
for comparison in comparisons:
    valid_bcs = []
    for bc in mtx.index.values:
        bc_id = re.search('(.*)_', bc).groups()[0]
        if bc_to_condition_map[bc_to_sample_map[bc_id]] in comparison:
            valid_bcs.append(bc)
    comparison_mtx = mtx.loc[valid_bcs, :]
    # t-sne stuff
    start = time.time()
    pca_operator = sklearn.decomposition.PCA(n_components=100, random_state=random_state)
    tsne_operator = sklearn.manifold.TSNE(n_components=2, random_state=random_state)
    Y_tsne = tsne_operator.fit_transform(pca_operator.fit_transform(comparison_mtx.values))
    end = time.time()
    print('Embedded t-SNE in {:.2f} seconds.'.format(end-start))
    # case control t-SNE
    cmap = ['dodgerblue', 'lightcoral']
    labels = []
    for bc in comparison_mtx.index.values:
        bc_id = re.search('^(S\d+)_', bc).groups()[0]
        labels.append(bc_to_condition_map[bc_to_sample_map[bc_id]])
    fig, ax = plt.subplots(figsize=(8, 8))
    scprep.plot.scatter2d(Y_tsne, label_prefix='t-SNE', title='PRRX1+ Mesenchymal Cells by Condition', legend_anchor=(1, 1),
                          c=labels, ticks=False, cmap=cmap, ax=ax, s=8, fontsize=14, discrete=True)
    plt.tight_layout()
    out_path = f'/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/plots/condition_pairs/he_prrx1-pos_{comparison[0]}_vs_{comparison[1]}_condition.pdf'
    plt.savefig(out_path)
    plt.clf()
    # individual genes t-SNE
    cmap = ['lightgrey', 'red']
    genes = ['CCL26', 'CCL11', 'CEBPB']
    for gene in genes:
        labels = comparison_mtx.loc[:, gene].values
        gene_max = 2*np.median(np.array(labels)[np.array(labels) > 0])
        labels = [label if label < gene_max else gene_max for label in labels]
        fig, ax = plt.subplots(figsize=(10, 8))
        scprep.plot.scatter2d(Y_tsne, label_prefix='t-SNE', title=f'{gene} in PRRX1+ Mesenchymal Cells',
                              c=labels, ticks=False, cmap=cmap, ax=ax, s=8, fontsize=14)
        plt.tight_layout()
        out_path = f'/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/plots/condition_pairs/he_prrx1-pos_{comparison[0]}_vs_{comparison[1]}_{gene}.pdf'
        plt.savefig(out_path)
        plt.clf()
