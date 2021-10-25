import pandas as pd
import numpy as np
import os
import re

inf = "1.20"
n = "10"

base_colfile_dir = "/Users/6j9/projects/mouse/cytoscape/coloring/"
in_file_path = base_colfile_dir + "cluster_mean_exp_order_cols_" + inf + ".tsv"
out_file_path = base_colfile_dir + "top_" + n + "_feature_cynotations_" + inf + ".tsv"

top_features_dir = "/Users/6j9/projects/mouse/irf-mcl/barcodes_by_cluster/top_features_" + inf

df = pd.read_csv(in_file_path, sep="\t", index_col=0, header=0)
col_name = "top_" + n + "_genes_" + inf
df.columns = [df.columns[0], col_name]

def get_top_n_genes_per_cluster(n, inf):
  top_n_genes = {}
  for r, d, f in os.walk(top_features_dir):
    for feat_file in f:
      if 'symb' in feat_file:
        cluster_num = re.match("symb_(\d+)_", feat_file).groups()[0]
        feat_path = os.path.join(r, feat_file)
        genes = []
        with open(feat_path) as in_file:
          for i, line in enumerate(in_file):
            if i < n:
              tmp = line.strip("\n").split("\t")[-2:]
              tmp[-1] = str(round(100 * float(tmp[-1]), 2))
              tmp = ": ".join(tmp)
              genes.append(tmp)
        top_n_genes[cluster_num] = " ".join(genes)
  return top_n_genes

top_n_genes = get_top_n_genes_per_cluster(int(n), inf)

def write_top_n_genes(df, top_n_genes):
  cluster_col = df.columns[0]
  for bc in df.index:
    cluster_num = str(df.loc[bc, cluster_col])
    value = ("* cluster " + cluster_num.ljust(3) + "* ")
    if len(cluster_num) == 1:
      cluster_num = "0" + cluster_num
    value = value + top_n_genes[cluster_num]
    df.loc[bc, col_name] = value
  df.to_csv(out_file_path, sep="\t", columns=[col_name])

write_top_n_genes(df, top_n_genes)
