from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import numpy as np
import pandas as pd
import scprep

prefix = "/Users/6j9/projects/mouse"
data_dir = prefix + "/data/"

cutoff, inf = "0.04", "1.20"
mcl_suffix = "mouseUPenn_iRF-LOOP_normalized_" + cutoff + "_mcl_clusters_" + inf + ".txt"
mcl_file = prefix + "/irf-mcl/mcl_clusters/" + mcl_suffix
pearson_out_path = prefix + "/irf-mcl/cluster_pearson_correlations/pearson_corr_matrix_" + inf + ".tsv"
hierarchical_groups_out_path = prefix + "/irf-mcl/cluster_pearson_correlations/hierarchical_groups_" + inf + ".txt"
cluster_col_out_path = prefix + "/cytoscape/coloring/cluster_mean_exp_order_cols_" + inf + ".tsv"

experiments = ['FGC2063_5_86846', 'FGC2063_5_86847', 'FGC2063_5_86848', 'FGC2063_5_86849', 'FGC2091_7_92848']
m = []
for experiment in experiments:
  print(experiment)
  matrix_dir = data_dir + experiment + '/filtered_feature_bc_matrix'
  mat = scprep.io.load_10X(matrix_dir, sparse=True, gene_labels='both')
  m.append(mat)

experiment_names = ['86846', '86847', '86848', '86849', '92848']
matrix, labels = scprep.utils.combine_batches(m, experiment_names, append_to_cell_names=True)
print("batches combined...")
del m

# assign clusters pre processing to use in cytoscape columns
cluster_count = 0 
cluster_labels = labels.copy()
keep = []
with open(mcl_file, 'rt') as in_file:
  for line in in_file:
    cluster_count += 1
    line_list = line.strip().split('\t')
    for item in line_list:
      cluster_labels[item] = cluster_count
      keep.append(item)

cluster_labels = cluster_labels[cluster_labels.index.isin(keep)]

matrix = scprep.filter.filter_rare_genes(matrix, min_cells=1)
matrix = scprep.normalize.library_size_normalize(matrix)
#mito_genes = np.array([g.startswith('mt-') for g in matrix.columns])
#mito_exp = matrix.loc[:, mito_genes].mean(axis=1)
#matrix = scprep.filter.filter_values(matrix, values=mito_exp, percentile=95, keep_cells='below')
matrix = scprep.transform.sqrt(matrix)
labels = labels[labels.index.isin(matrix.index)]
print("preprocessing done...")

# reassign clusters post processing to use in pearson corr
cluster_count = 0 
with open(mcl_file, 'rt') as in_file:
  for line in in_file:
    cluster_count += 1
    line_list = line.strip().split('\t')
    for item in line_list:
      if item in labels:
        labels[item] = cluster_count

print("cluster ids done...")

num_genes = np.shape(matrix)[1]
mean_exp_vectors = np.zeros((num_genes, cluster_count))
for i in range(cluster_count):
  print(i)
  cluster_matrix = matrix.loc[labels.index[labels == i+1]]
  vector = np.mean(cluster_matrix, axis=0)
  mean_exp_vectors[:, i] = vector

df = pd.DataFrame(data=mean_exp_vectors, index=matrix.columns, columns=np.arange(1, cluster_count + 1).tolist())
df.to_csv(mean_exp_out_path, sep="\t")

pearson = np.corrcoef(mean_exp_vectors, rowvar=False)
pearson = (pearson + pearson.T)/2
np.fill_diagonal(pearson, 1)
idx = np.arange(1, cluster_count + 1)
df = pd.DataFrame(data=pearson, index=idx, columns=idx)
df.to_csv(pearson_out_path, sep="\t")
print("pearson correlation done...")

def hierarchical_clustering(pearson=pearson):
  dissimilarity = 1 - pearson
  hierarchy = linkage(squareform(dissimilarity), method='ward')
  assignments = fcluster(hierarchy, 0.5, criterion='distance')
  num_groups = max(assignments)
  groups = [[] for i in range(num_groups)]
  for i in range(cluster_count):
    groups[assignments[i] - 1].append(i+1)
  groups = [np.asarray(group) for group in groups]
  groups = np.asarray(groups)
  groups = groups[np.argsort([len(group) for group in groups])[::-1]]
  return groups

cluster_groups = hierarchical_clustering()
with open(hierarchical_groups_out_path , "w") as out_file:
  out_file.write("\n".join(str(group.tolist()) for group in cluster_groups))

def sort_hierarchical_clusters(cluster_groups=cluster_groups, pearson=pearson):
  np.fill_diagonal(pearson, 0)
  new_cluster_groups = []
  for cluster_group in cluster_groups:
    cluster = cluster_group[0] - 1
    pearson[:, cluster] = 0
    new_cluster_group = [cluster]
    num_clusters = len(cluster_group)
    cluster_inds = cluster_group - 1
    while len(new_cluster_group) < num_clusters:
      col = pearson[cluster, :].copy()
      mask = np.ones(len(col), dtype=bool)
      mask[cluster_inds] = False
      col[mask] = 0
      next_cluster = np.where(col == max(col))[0][0]
      new_cluster_group.append(next_cluster)
      pearson[:, next_cluster] = 0
    new_cluster_groups.append(new_cluster_group)
  return(new_cluster_groups)

cluster_order = sort_hierarchical_clusters()
print(cluster_order)
cluster_order = list(np.concatenate(cluster_order))
print(cluster_order)
print("\nclustering done...")

def write_cyto_cluster_order_cols(cluster_order=cluster_order, cluster_labels=cluster_labels):
  cluster_labels = cluster_labels.to_frame()
  start_order_name = "start_order_" + inf
  end_order_name = "end_order_" + inf
  cluster_labels.columns = [start_order_name]
  cluster_labels[end_order_name] = np.full(np.shape(cluster_labels)[0], 0)
  for bc in cluster_labels.index:
    idx = cluster_labels.loc[bc, start_order_name] - 1
    cluster_labels.loc[bc, end_order_name] = str(cluster_order.index(idx) + 1)
  return cluster_labels

cluster_labels = write_cyto_cluster_order_cols(cluster_order, cluster_labels)
cluster_labels.to_csv(cluster_col_out_path, sep="\t")
print("cluster labels saved...")
