import numpy as np

mcl_dir = "/Users/6j9/projects/mouse/phate/mcl/"
mcl1 = mcl_dir + "mouseUPenn_iRF-LOOP_normalized_mcl_clusters_1.2.txt"
mcl2 = mcl_dir + "mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.25.txt"
#mcl2 = mcl_dir + "mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.3.txt"
mcl_files = [mcl1, mcl2]

processed_mcls = []
for mcl_file in mcl_files:
  mcl_clusters = []
  cluster_count = 0 
  with open(mcl_file, 'rt') as in_file:
    for line in in_file:
      cluster_count += 1
      mcl_clusters.append(line.strip().split('\t'))
  processed_mcls.append(mcl_clusters)

mcl1_clusters, mcl2_clusters = processed_mcls[0], processed_mcls[1]

def get_best_n_composite_cluster_sizes(matching_clusters, n):
  best_n_composite_clusters = {}
  for i in range(n):
    cluster, num_matches = list(matching_clusters.keys()), list(matching_clusters.values())
    best_cluster = cluster[num_matches.index(max(num_matches))]
    best_max = matching_clusters[best_cluster]
    best_n_composite_clusters[best_cluster] = best_max
    matching_clusters.pop(best_cluster)
    if len(matching_clusters) == 0:
      break
  return(best_n_composite_clusters)

composite_clusters = []
num_clusters = 5
for i, cluster1 in enumerate(mcl1_clusters):
  matching_clusters = {}
  size = len(cluster1)
  for j, cluster2 in enumerate(mcl2_clusters):
    if len(cluster2) <= size:
      intersection = np.intersect1d(cluster1, cluster2).size
      if not intersection == 0:
        pct_matches = intersection / size
        name = str(i + 1) + '-' + str(j + 1)
        matching_clusters[name] = pct_matches
  if len(matching_clusters) > 0:
    composite_clusters.append(get_best_n_composite_cluster_sizes(matching_clusters, num_clusters))
  print(f"cluster {i+1} barcode comparison done...")

for cluster_set in composite_clusters:
  print(f"{list(cluster_set.keys())}: {sum(list(cluster_set.values()))} {list(cluster_set.values())}")
