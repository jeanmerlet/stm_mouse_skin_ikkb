import pandas as pd
import numpy as np
import scprep
import re

celltypes_in_path = "/Users/6j9/projects/mouse/annotation/cellassign/celltypes/mod1/combined_cell_types_mod1.tsv"
celltypes_out_path = "/Users/6j9/projects/mouse/annotation/cellassign/celltypes/mod1/posthoc_combined_cell_types_mod1.tsv"
data_dir = "/Users/6j9/projects/mouse/data/"

experiments = ['FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848']
m = []
for experiment in experiments:
  print(experiment)
  matrix_dir = data_dir + experiment + '/filtered_feature_bc_matrix'
  m.append(scprep.io.load_10X(matrix_dir, gene_labels='both'))

experiment_names = ['86846', '86847', '86848', '86849', '92848']
matrix, labels = scprep.utils.combine_batches(m, experiment_names, append_to_cell_names=True)
del m

symb_ens = matrix.columns
symbols = symb_ens.tolist()
for i, notation in enumerate(symb_ens):
  symbols[i] = (re.match("(.+) \(", notation).groups()[0]).upper()

matrix.columns = symbols
print("batches combined...")

def check_expression(genes, exp_vector):
  for gene in genes:
    print("meow")
    #if


check_types = ["Dendritic cells", "Gamma delta T cells", "ILC1", "ILC2", "T cells"]
check_genes = [["CD3G", "CD19", "LY6G", "LY6C"], ["CD4", "CD8A"], ["CD3G"], ["CD3G"], ["TRDV4"]]
df = pd.read_csv(celltypes_in_path, sep="\t", header=0, index_col=0)

for bc in df.index:
  break
  cell_type = df.loc[bc]
  if cell_type in check_types:
    genes = check_genes[check_types.index(cell_type)]
    exp_vector = matrix.iloc[np.where(labels == bc), :]
    if check_expression(genes, exp_vector):
      print("meow")
      # change the cell type

#df.to_frame()
#df.to_csv(celltypes_out_path, sep="\t")
