import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import sklearn.decomposition
import sklearn.manifold
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform, pdist
import scprep
import phate

prefix = "/Users/6j9/projects/mouse"
data_dir = prefix + "/irf-mcl/cluster_pearson_correlations"
cutoff= "0.04"

#fig, axes = plt.subplots(3, 1, figsize=(15, 9)) 
#fig, axes = plt.subplots(1, 3, figsize=(15, 5)) 
fig, axes = plt.subplots(1, 1, figsize=(15, 6)) 

infs = ["1.20", "1.25", "1.30"]
cluster_counts = [74, 170, 316]

group_names = [["fibro", "immune", "other tissue", "adipo-myo", "single", "single", "single", "single"],
               ["fibro 1", "fibro 2", "immune 1", "other tissue 1", "immune 2", "immune 3", "other tissue 2",
                "immune 4", "dead", "mast", "adipo-myo", "single", "single"],
               ["fibro 1", "fibro 2", "fibro 3", "immune 1", "immune 2", "immune 3", "keratinocytes",
                "other tissue 1", "immune 4", "other tissue 2", "immune 5", "dead 1", "adipo-myo",
                "dead 2", "mast", "single", "single", "single"]]

#def subset_clusters_by_group_type(clusters, groups, group_type):
  #for group in groups:
    #if re.match("(group_type)", group):
      

for idx, inf in enumerate(infs):
  #if idx < 2:
    #continue
  print(inf)
  pearson_path = data_dir + "/pearson_corr_matrix_" + inf + ".tsv"
  groups_path = data_dir + "/hierarchical_groups_" + inf + ".txt"
  cluster_count = cluster_counts[idx]

  pearson = pd.read_csv(pearson_path, sep="\t", header=0, index_col=0)
  similarity = 1 - pearson

  groups = []
  with open(groups_path) as in_file:
    for line in in_file:
      groups.append([int(x) for x in line.strip("\n").strip('][').split(', ')])

  method = 'ward'
  threshold = 0.5

  #ax = fig.axes[idx]
  ax = fig.axes[0]
  cm1 = plt.get_cmap('tab10')
  cm2 = plt.get_cmap('tab20')
  colors = [cm1(i) for i in range(10)] + [cm2((j*2)+1) for j in range(10)]
  color_hex = [mpl.colors.rgb2hex(x) for x in colors]
  sch.set_link_color_palette(color_hex)
  ax.set_title(f"Mean Exp Dendrogram ({method}) - inf {inf}")
  hierarchy = sch.linkage(squareform(similarity), method)
  #c, coph_dists = sch.cophenet(hierarchy, pdist(similarity))
  #print(c)
  dendro = sch.dendrogram(hierarchy, ax=ax, orientation='top', color_threshold=threshold,
                          labels=similarity.index, above_threshold_color='black')

  ticks = ax.get_xticklabels()
  groups_idx = []
  for i in ticks:
    groups_idx.append(int(i.get_text()) - 1)

  cluster_names = group_names[idx]
  cluster_labels = np.full(cluster_count, "", dtype='U20')
  for i, group in enumerate(groups):
    for cluster in group:
      cluster_labels[groups_idx.index(cluster - 1)] = cluster_names[i]

  new_labels = ['threshold']
  for label in cluster_labels:
    if label not in new_labels and label != 'single':
      new_labels.append(label)

  cluster_labels = new_labels + ['single']

  ax.axhline(y=threshold, c='black', lw=1, linestyle='dashed')
  ax.legend(labels=cluster_labels, bbox_to_anchor=(1.04, 1), loc='upper left', ncol=1)
  sch.set_link_color_palette(None)
  break

plt.tight_layout()
plt.show()
