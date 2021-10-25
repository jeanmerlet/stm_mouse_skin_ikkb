import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import re

features_path = '/Users/6j9/projects/mouse/data/irf-mcl/importance_scores/top_features_1.30/top_genes_full_1.30.tsv'
out_dir = '/Users/6j9/projects/mouse/plots/final/importance_bargraphs'

#clusters = [1, 2, 4, 7, 8, 3, 5, 6, 9, 10]
clusters = np.arange(1, 317, 1)

features = pd.read_csv(features_path, sep='\t', index_col=None)
features = features.sort_values('cluster')
features.index = features['cluster']

num_genes = 5
for cluster in clusters:
  print(cluster)
  fig, ax = plt.subplots(1, 1, figsize=(5, 5))
  scores = features.loc[cluster, 'gene-score_1.30'].split(',')
  if len(scores) < num_genes:
    x = np.arange(len(scores))
  else:
    x = np.arange(num_genes)
  #scores = features.loc[cluster, 'gene-score_1.30'].split(',')[:5]
  scores = scores[:num_genes]
  mpl_scores = []
  genes = []
  for score in scores:
    genes.append(re.search('(.+): ', score).group(1))
    mpl_scores.append(float(re.search('(\d+\.\d+)', score).group(1)))
  ax.bar(x, mpl_scores)
  ax.set_ylim([0, 100])
  ax.set_xticks(x)
  ax.set_xticklabels(genes, rotation=90)
  ax.set_ylabel('Importance Score')
  ax.set_xlabel('Genes')
  ax.set_title('Cluster ' + str(cluster))
  plt.tight_layout()
  cluster_name = '0' * (3 - len(str(cluster))) + str(cluster)
  filename = '1.30_inf_top5-imp-scores_cluster' + cluster_name + '_bar.pdf'
  out_path = os.path.join(out_dir, filename)
  plt.savefig(out_path)
  plt.clf()
  #plt.show()
