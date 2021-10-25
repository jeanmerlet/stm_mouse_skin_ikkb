import numpy as np
import pandas as pd
import scprep
import scrublet as scr

data_dir = "/Users/6j9/projects/mouse/data/"
mcl_file = "/Users/6j9/projects/mouse/phate/mcl/mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.25.txt"

experiments = ['FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848']
#experiments = ['FGC2063_5_86846', 'FGC2063_5_86847']
doublet_thresholds = [0.27, 0.25, 0.24, 0.2, 0.2]
#doublet_thresholds = [0.27, 0.25]

m, all_predicted_doublets = [], []
for i, experiment in enumerate(experiments):
  matrix_dir = data_dir + experiment + '/filtered_feature_bc_matrix'
  mat = scprep.io.load_10X(matrix_dir, sparse=True, gene_labels='both')
  scrub = scr.Scrublet(mat, expected_doublet_rate=0.06)
  doublet_scores, predicted_doublets = scrub.scrub_doublets()
  scrub.call_doublets(threshold=doublet_thresholds[i])
  all_predicted_doublets.append(predicted_doublets)
  m.append(mat)

predicted_doublets = np.concatenate(all_predicted_doublets)
experiment_names = ['86846', '86847', '86848', '86849', '92848']
#experiment_names = ['86846', '86847']
matrix, labels = scprep.utils.combine_batches(m, experiment_names, append_to_cell_names=True)
del m

# remove unexpressed genes
print('rare genes')
matrix = scprep.filter.filter_rare_genes(matrix, min_cells=1)

# old upper library size cutoff method
# matrix_pp = scprep.filter.filter_library_size(matrix, cutoff=25000, keep_cells='below')
# doublet/multiplet detection via Scrublet
print('doublets')
doublets_barcodes = np.array(matrix[predicted_doublets == True].index)

# normalize
print('normalize')
matrix_ln = scprep.normalize.library_size_normalize(matrix)
del matrix

# remove high mitochondrial gene expression cells
print('mito')
mitochondrial_gene_list = np.array([gene.startswith('mt-') for gene in matrix_ln.columns])
mito_expression = matrix_ln.loc[:, mitochondrial_gene_list].mean(axis=1)
mito_barcodes = np.array(scprep.filter.filter_values(matrix_ln, values=mito_expression, percentile=95, keep_cells = 'below').index)

# square-root transform for euclidean distance to work
#matrix_pp = scprep.transform.sqrt(matrix_pp)

print(len(np.intersect1d(doublets_barcodes, mito_barcodes)))
int_barcodes = np.array(np.intersect1d(doublets_barcodes, mito_barcodes))

print('combining')
barcodes = np.array(matrix_ln.index)
cyto_tags = np.full(len(barcodes), "mito", dtype="U12")
cyto_tags[np.in1d(barcodes, mito_barcodes) == True] = "preprocessed"
cyto_tags[np.in1d(barcodes, doublets_barcodes) == True] = "doublets"
cyto_tags[np.in1d(barcodes, int_barcodes) == True] = "both"
out = np.vstack((barcodes, cyto_tags)).T

filepath = "/Users/6j9/projects/mouse/cytoscape/pp_cytoscape_cols.txt"
np.savetxt(filepath, out, delimiter=",", fmt="%s")
