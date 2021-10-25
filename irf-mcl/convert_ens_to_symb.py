import pandas as pd
import numpy as np
import os
import re

inf = "1.30"

in_file_dir = "/Users/6j9/projects/mouse/irf-mcl/feature_effect/top_features_" + inf
out_file_dir = in_file_dir.replace("top_features", "top_gene_symbols")
ens_to_symb_path = "/Users/6j9/projects/mouse/irf-mcl/feature_effect/ens_to_symb.tsv"

ens_to_symb = pd.read_csv(ens_to_symb_path, sep="\t", index_col=0, header=None)
ens_to_symb.columns = ['symbol']
ens_to_symb['symbol'] = ens_to_symb['symbol'].str.upper()

feature_files = []
symbol_files = []
for r, d, f in os.walk(in_file_dir):
  for feature_file in f:
    if 'topFeatures' in feature_file:
      feature_files.append(os.path.join(r, feature_file))
      symbol_file = feature_file.replace("topFeatures", "top_gene_symbols")
      symbol_file = symbol_file.replace(".txt", ".tsv")
      symbol_files.append(os.path.join(out_file_dir, symbol_file))

for i, f in enumerate(feature_files):
  print(i)
  df = pd.read_csv(f, sep="\t", index_col=0, header=None)
  df.columns = ["cluster_number", "importance_score"]
  df["symbols"] = ens_to_symb.loc[df.index]
  df = df.set_index("symbols")
  df.to_csv(symbol_files[i], sep="\t")

