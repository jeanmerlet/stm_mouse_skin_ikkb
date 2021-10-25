import pandas as pd
import numpy as np
import os

inf = "1.30"

convert_path = "/Users/6j9/projects/mouse/irf-mcl/barcodes_by_cluster/ens_to_symb.tsv"
target_dir = "/Users/6j9/projects/mouse/irf-mcl/barcodes_by_cluster/top_features_" + inf
out_path = target_dir + "/top_genes_full_" + inf + ".tsv"

convert = pd.read_csv(convert_path, index_col=0, header=None, sep="\t")
convert.iloc[:, 0] = convert.iloc[:, 0].str.upper()

feat_list = []
for r, d, f in os.walk(target_dir):
  for feat_file in f:
    if 'topFeatures' in feat_file:
      feat_list.append(os.path.join(r, feat_file))

def round_score(score):
  return str(round(100 * score, 2))

for i, f in enumerate(feat_list):
  print(i)
  df = pd.read_csv(f, index_col=None, header=None, sep="\t")
  ens = df.iloc[:, 0]
  symb = np.full(len(ens), "", dtype='U20')
  for j, gene in enumerate(ens):
    symb[j] = convert.loc[gene].values[0] + ": " + round_score(df.iloc[:, 2][j])
  df.iloc[:, 0] = symb
  df = df.iloc[:, [1, 0]]
  df.columns = ['cluster', 'gene-score']
  df = df.groupby(['cluster']).agg(', '.join)
  if i == 0:
    df.to_csv(out_path, sep="\t", index=True, header=True)
  else:
    df.to_csv(out_path, sep="\t", index=True, header=None, mode='a')
