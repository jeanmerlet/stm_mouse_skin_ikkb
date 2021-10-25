import numpy as np

mcl_dir = "/Users/6j9/projects/mouse/irf-mcl/mcl_clusters"
mcl1 = mcl_dir + "mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.25.txt"
mcl2 = mcl_dir + "mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.3.txt"
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
all_similarities = []
for i, cluster1 in enumerate(mcl1_clusters):
  similarities = {}
  for j, cluster2 in enumerate(mcl2_clusters):
    max_matches = max(len(cluster1), len(cluster2))
    intersection = np.intersect1d(cluster1, cluster2)
    if not intersection.size == 0:
      similarity = intersection.size / max_matches
      name = str(i + 1) + '-' + str(j + 1)
      similarities[name] = similarity
  all_similarities.append(similarities)
  print(f"cluster {i+1} barcode comparison done...")
  if i == 19:
    break

for similarities in all_similarities:
  best_cluster = max(similarities)
  best_max = similarities[best_cluster]
  print(best_cluster)
  print(best_max)
  print("")
