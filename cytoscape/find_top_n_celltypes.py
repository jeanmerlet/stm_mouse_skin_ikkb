import pandas as pd
import numpy as np

path = '/Users/6j9/projects/mouse/cytoscape/coloring/final_celltype_cols/top8_v2.tsv'

celltypes = pd.read_csv(path, sep='\t', names=['cell', 'type'], index_col=None)
grouped = celltypes.groupby('type')['cell'].apply(list).apply(len)
print(grouped.sort_values())
