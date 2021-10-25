import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import sklearn.decomposition
import sklearn.manifold
import scipy.cluster.hierarchy as sch
import scprep
import phate

prefix = "/Users/6j9/projects/mouse"
data_dir = prefix + "/irf-mcl/cluster_pearson_correlations"
cutoff= "0.04"

#fig, axes = plt.subplots(3, 3, figsize=(15, 9)) 
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

infs = ["1.20", "1.25", "1.30"]
#group_names = [["fibro", "immune", "other tissue", "adipo-myo", "mast", "plasma", "neutro", "dead"],
#               ["fibro 1", "fibro 2", "immune 1", "other tissue 1", "immune 2", "immune 3", "other tissue 2",
#                "immune 4", "dead", "mast", "adipo-myo", "plasma", "neutrophils"],
#               ["fibro 1", "fibro 2", "fibro 3", "immune 1", "immune 2", "immune 3", "keratinocytes",
#                "other tissue 1", "immune 4", "other tissue 2", "immune 5", "dead", "adipo-myo",
#                "dead", "mast", "neutrophils", "fibro 4", "plasma"]]
group_names = [["fibro", "immune", "other tissue", "adipo-myo", "single", "single", "single", "single"],
               ["fibro 1", "fibro 2", "immune 1", "other tissue 1", "immune 2", "immune 3", "other tissue 2",
                "immune 4", "dead", "mast", "adipo-myo", "plasma", "neutrophils"],
               ["fibro 1", "fibro 2", "fibro 3", "immune 1", "immune 2", "immune 3", "keratinocytes",
                "other tissue 1", "immune 4", "other tissue 2", "immune 5", "dead", "adipo-myo",
                "dead", "mast", "dead", "dead", "dead"]]

dendro_groups = [[2, 3, 5, 6, 7, 9, 10, 12, 14, 15, 17, 25, 26, 27, 32, 33, 39, 41, 44, 45, 49, 52, 53, 56, 59, 60, 62, 63, 65, 66, 67, 68, 70, 72, 73], 
[1, 4, 8, 11, 18, 20, 21, 28, 30, 35, 37, 40, 42, 43, 46, 47, 54, 57, 61, 64, 69], 
[13, 16, 24, 29, 31, 34, 36, 50, 51, 55, 71, 74], 
[19, 38], 
[58],
[22],
[48],
[23]]

# 1.20
#colors = ['green', 'blue', 'red', 'orange', 'black', 'black', 'black', 'black']

# 1.30
infs = ['1.30']
colors = ['blue', 'green', 'orange', 'turquoise', 'grey', 'magenta', 'lawngreen', 'firebrick', 'orchid', 'darkorange', 'purple', 'yellow', 'brown', 'red', 'cornflowerblue', 'black', 'black', 'black']

for idx, inf in enumerate(infs):
  #if idx < 2:
    #continue
  mean_exp_vectors_path = data_dir + "/mean_exp_vectors_" + inf + ".tsv"
  groups_path = data_dir + "/hierarchical_groups_" + inf + ".txt"

  mean_exp_vectors = pd.read_csv(mean_exp_vectors_path, sep="\t", header=0, index_col=0)
  mean_exp_vectors = np.transpose(mean_exp_vectors)
  cluster_count = np.shape(mean_exp_vectors)[0]

  groups = []
  with open(groups_path) as in_file:
    for line in in_file:
      groups.append([int(x) for x in line.strip("\n").strip('][').split(', ')])

  cluster_names = group_names[idx]
  cluster_labels = np.full(cluster_count, "", dtype='U20')
  for i, group in enumerate(groups):
    for cluster in group:
      cluster_labels[cluster - 1] = cluster_names[i]

  #cluster_labels = [int(x) for x in cluster_labels]
  
  for i in range(3):
    #ax = fig.axes[(idx*3) + i]
    ax = fig.axes[i]
    #if i % 3 == 0:
    #if i == 0:
    if False:
      # PCA
      title = 'Mean Exp PCA - inf ' + inf
      pca_operator = sklearn.decomposition.PCA(n_components=2)
      Y_pca = pca_operator.fit_transform(np.array(mean_exp_vectors))
      scprep.plot.scatter2d(Y_pca, label_prefix="PC", title=title, ticks=False, ax=ax, c=cluster_labels, legend=False)
      for i, ((x, y),) in enumerate(zip(Y_pca)):
        break
        value = i + 1
        color_id = next(i for i, v in enumerate(dendro_groups) if value in v)
        ax.text(x, y, i+1, ha="center", va="center", color=colors[color_id])
      #scprep.plot.scatter2d(Y_pca, label_prefix="PC", title=title, ticks=False, ax=ax, c=cluster_labels, legend=False)
      print("PCA done...")
    #elif i % 3 == 1:
    #elif i == 1:
    elif True:
      # t-SNE
      title = 'Mean Exp t-SNE - inf ' + inf
      #pca_operator = sklearn.decomposition.PCA(n_components=100)
      tsne_operator = sklearn.manifold.TSNE(n_components=2, learning_rate=100)
      Y_tsne = tsne_operator.fit_transform(np.array(mean_exp_vectors))
      #Y_tsne = tsne_operator.fit_transform(pca_operator.fit_transform(np.array(mean_exp_vectors)))
      scprep.plot.scatter2d(Y_tsne, label_prefix="t-SNE", title=title, ticks=False, ax=ax, c=cluster_labels, legend=False)
      for i, ((x, y),) in enumerate(zip(Y_tsne)):
        break
        value = i + 1
        color_id = next(i for i, v in enumerate(dendro_groups) if value in v)
        ax.text(x, y, i+1, ha="center", va="center", color=colors[color_id])
      print("t-SNE done...")
    else:
      # PHATE
      title = 'Mean Exp PHATE - inf ' + inf
      phate_operator = phate.PHATE(n_jobs=-2)
      Y_phate = phate_operator.fit_transform(mean_exp_vectors)
      scprep.plot.scatter2d(Y_phate, label_prefix="PHATE", title=title, ticks=False, ax=ax, c=cluster_labels,
                            legend=True, legend_anchor=(1.04, 1))
      for i, ((x, y),) in enumerate(zip(Y_phate)):
        break
        value = i + 1
        color_id = next(i for i, v in enumerate(dendro_groups) if value in v)
        ax.text(x, y, i+1, ha="center", va="center", color=colors[color_id])
      print("PHATE done...")
  print(f"\nplots for {inf} inflation done...")
  break

plt.tight_layout()
plt.show()
