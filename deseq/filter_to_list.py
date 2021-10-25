import pandas as pd
import numpy as np

filter_path = '/Users/6j9/projects/mouse/data/de_results/dana_list.txt'

table_path = '/Users/6j9/projects/mouse/data/de_results/cre-pos_v_cre-neg_fibroblasts_prrx1-high/table_prrx1-high_cre-pos_v_cre-neg.tsv'
out_path = '/Users/6j9/projects/mouse/data/de_results/cre-pos_v_cre-neg_fibroblasts_prrx1-high/filtered_table_prrx1-high_cre-pos_v_cre-neg.tsv'

genes = []
with open(filter_path, 'r') as in_file:
    for line in in_file:
        genes.append(line.strip().upper())

table = pd.read_csv(table_path, sep='\t', header=0, index_col=0)
list_idx = np.in1d(table.index.values, genes)

not_in_table = np.array(genes)[np.in1d(genes, table.index.values, invert=True)]
for gene in not_in_table:
    print(gene)

table.loc[list_idx, :].to_csv(out_path, sep='\t')
