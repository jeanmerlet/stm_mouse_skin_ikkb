import numpy as np

mcl_path = "/Users/6j9/projects/mouse/irf-mcl/mcl_clusters/mouseUPenn_iRF-LOOP_normalized_0.04_mcl_clusters_1.30.txt"
out_path = "/Users/6j9/projects/mouse/cytoscape/coloring/cluster_mean_exp_order_cols_1.30.tsv"

cluster_count = 0 
mcl_cluster_ids = {}
bc_col = []
with open(mcl_path, 'rt') as in_file:
  for line in in_file:
    line_list = line.strip().split('\t')
    for barcode in line_list:
      mcl_cluster_ids[barcode] = cluster_count
      bc_col.append(barcode)
    cluster_count += 1
bc_col = np.asarray(bc_col)

# 1.20
# cluster_order = np.asarray([10, 41, 56, 68, 70, 73, 35, 23, 15, 0, 46, 53, 4, 2, 1, 48, 51, 59, 62, 67, 69, 71, 72, 16, 52, 55, 24, 26, 9, 5, 31, 40, 13, 32, 44, 14, 6, 3, 29, 36, 39, 60, 20, 18, 37, 8, 65, 27, 25, 19, 28, 49, 7, 61, 43, 11, 66, 34, 42, 58, 45, 17, 54, 33, 64, 12, 50, 30, 38, 63, 22, 57, 47, 21])

# 1.30
cluster_order = np.asarray([16, 19, 22, 24, 26, 27, 12, 30, 31, 32, 23, 35, 36, 39, 40, 38, 41, 42, 10, 43, 44, 29, 45, 46, 9, 47, 48, 49, 50, 13, 17, 34, 51, 20, 52, 53, 14, 55, 56, 11, 25, 57, 33, 58, 0, 37, 59, 18, 60, 8, 61, 5, 21, 28, 62, 54, 63, 6, 64, 65, 3, 7, 66, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 147, 148, 149, 150, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 1, 2, 4, 15, 67, 68, 145, 146, 151, 152, 153, 231, 232, 245, 246, 278, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 279])

start_clusters, end_clusters = [], []
for barcode in bc_col:
  start_clusters.append(mcl_cluster_ids[barcode] + 1)
  end_clusters.append(np.nonzero(cluster_order == mcl_cluster_ids[barcode])[0][0] + 1)
start_clusters = np.asarray(start_clusters)
end_clusters = np.asarray(end_clusters)
cols = np.vstack((bc_col, start_clusters, end_clusters)).T
#cols = np.vstack((bc_col, start_clusters)).T

np.savetxt(out_path, cols, delimiter="\t", fmt="%s")
