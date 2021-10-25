import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.decomposition
import scprep

#inf = "1.20"
#corr_path = "/Users/6j9/projects/mouse/irf-mcl/cluster_pearson_correlations/pearson_corr_matrix_" + inf + ".tsv"
#groups_path = "/Users/6j9/projects/mouse/irf-mcl/cluster_pearson_correlations/hierarchical_groups_" + inf + ".txt"
#out_fig_path = "/Users/6j9/projects/mouse/plots/out/irf_corr_clusters_pca/clusters_pca_" + inf + ".png"
out_fig_path = "/Users/6j9/projects/mouse/plots/out/irf_corr_clusters_pca/clusters_pca_all.png"

fig, axes = plt.subplots(1,3, figsize=(15, 5)) 

infs = ["1.20", "1.25", "1.30"]
for ax_idx, inf in enumerate(infs):
  corr_path = "/Users/6j9/projects/mouse/irf-mcl/cluster_pearson_correlations/pearson_corr_matrix_" + inf + ".tsv"
  groups_path = "/Users/6j9/projects/mouse/irf-mcl/cluster_pearson_correlations/hierarchical_groups_" + inf + ".txt"
  corr = pd.read_csv(corr_path, sep="\t", index_col=0)

  groups = []
  with open(groups_path) as in_file:
    for line in in_file:
      groups.append([int(x) for x in line.strip("\n").strip('][').split(', ')])

  cluster_labels = np.zeros(np.shape(corr)[0])
  for i, group in enumerate(groups):
    for cluster in group:
      cluster_labels[cluster - 1] = i + 1

  cluster_labels = [int(x) for x in cluster_labels]
  pca_operator = sklearn.decomposition.PCA(n_components=2)
  Y_pca = pca_operator.fit_transform(np.array(corr))

  ax = fig.axes[ax_idx]
  title = 'iRF Clusters Corr PCA - inf ' + inf
  scprep.plot.scatter2d(Y_pca, label_prefix="PC", title=title, ticks=False, ax=ax, c=cluster_labels)
  plt.tight_layout()

#plt.savefig(out_fig_path)
plt.show()
