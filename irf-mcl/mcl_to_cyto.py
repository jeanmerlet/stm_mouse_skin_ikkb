from os import listdir
from os.path import isfile, join

data_dir = "/Users/6j9/projects/mouse/irf-mcl"
total_edge_file = "/Users/6j9/projects/mouse/irf-mcl/mcl_edgefiles/mouseUPenn_iRF-LOOP_normalized_sorted.txt"
# Files of mcl clusters - 1 cluster per line, tab delimited spaces between items in a cluster
mcl_dir = join(data_dir, "mcl_clusters")
mcl_files = [f for f in listdir(mcl_dir) if isfile(join(mcl_dir, f))]
mcl_files = ["/Users/6j9/projects/mouse/irf-mcl/mcl_clusters/mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.20_reordered.txt"]
# Files of edge lists for the given clusters
out_dir = join(data_dir, "mcl_edgefiles")
out_files = [mcl_file.replace("mcl", "cyto") for mcl_file in mcl_files]
out_files = ["mouseUPenn_iRF-LOOP_normalized_0.04_cyto_clusters_1.20_reordered.txt"]

processed_mcls = []
for mcl_file in mcl_files:
  cluster_count = 0
  mcl_cluster_ids = {}
  with open(join(mcl_dir, mcl_file), 'rt') as in_file:
    for line in in_file:
      cluster_count += 1
      line_list = line.strip().split('\t')
      for item in line_list:
        mcl_cluster_ids[item] = cluster_count
  processed_mcls.append(mcl_cluster_ids)

# Important!! - this code assumes the first line of the original edge file is a header
max_edges = 11084**2
for i, mcl_cluster_ids in enumerate(processed_mcls):
  line_count = 0
  with open(join(out_dir, out_files[i]), 'wt') as out_file:
    with open(total_edge_file, 'rt') as in_file:
      for line in in_file:
        line_count += 1
        if line_count != 1:
          line_list = line.strip().split('\t')
          if mcl_cluster_ids.get(line_list[0]) == mcl_cluster_ids.get(line_list[1]):
            out_file.write(line)
        # Need to cut off at original mcl input file cutoff %
        cutoff = round(max_edges * 0.0004)
        if line_count == cutoff:
          break
  print(f"finished {out_files[i]}")
