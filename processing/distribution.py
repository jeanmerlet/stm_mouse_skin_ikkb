import csv
import gzip
import os
import scipy
import scipy.io
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt

experiment = 'FGC2091_7_92847'
matrix_dir = './data/' + experiment + '/raw_feature_bc_matrix'
m = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))

barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode='rt'), delimiter="\t")]

genes_path = os.path.join(matrix_dir, "features.tsv.gz")
gene_ids = [row[0] for row in csv.reader(gzip.open(genes_path, mode='rt'), delimiter="\t")]

barcodes = np.asarray(barcodes)
gene_ids = np.asarray(gene_ids)

nrows, ncols = np.shape(m)[0], np.shape(m)[1]

# number of cells expressing at least 1 gene (y axis):
nnz_cols = m.getnnz(axis=0)
nonzero_cols = nnz_cols[nnz_cols > 1]
#print(len(nonzero_cols))

# max number of genes expressed by one cell (x axis):
max_genes = max(nnz_cols)

def cells_lesseq_genes(x):
  return len(nnz_cols[nnz_cols <= x])

def cells_greateq_genes(x):
  return len(nnz_cols[nnz_cols >= largest_hund[x] - x])

#for i in range(50):
  #x1 = len(nnz_cols) - cells_lesseq_genes(i)
  #print(f"{i}: {x1}")
  #if i > 0:
    #x2 = len(nnz_cols) - cells_lesseq_genes(i-1)
    #print(f"{i} - {i-1}: {x2 -x1}")

#print("")


split, split1, ceil = 24, 300, 1001
cols = nonzero_cols[nonzero_cols <= split]
cols1 = nonzero_cols[nonzero_cols <= split1]
cols1 = cols1[cols1 > split]
cols2 = nonzero_cols[nonzero_cols <= ceil]
cols2 = cols2[cols2 > split1]
cols3 = nonzero_cols[nonzero_cols > ceil]

fig = plt.figure()
fig.suptitle(experiment)

ax = fig.add_subplot(221)
ax.title.set_text(f"2 - {split} genes expressed ({len(cols)} cells)")
ax.set_ylabel('number of cells')
ax.set_xlabel('# of genes expressed')
ax.set_xticks(np.arange(2, split, step=4))
ax.hist(cols, bins=split-1)

ax1 = fig.add_subplot(222)
ax1.title.set_text(f"{split+1} - {split1} genes expressed ({len(cols1)} cells)")
ax1.set_ylabel('number of cells')
ax1.set_xlabel('# of genes expressed')
ax1.set_xticks(np.arange(split+1, split1, step=25))
ax1.hist(cols1, bins=(split1-split-1))
ax1.set_ylim([0, 90])

y = cols1
x = np.linspace(split+1, split1, 25)
#dist_names = ['gamma', 'beta', 'rayleigh', 'norm', 'pareto']
dist_names = ['pareto']
for dist_name in dist_names:
  distribution = getattr(st, dist_name)
  params = distribution.fit(y)
  args = params[:-2]
  mean = params[-2]
  std = params[-1]

  fitted_pdf = distribution.pdf(x, loc=mean, scale=std, *args) * len(cols1)
  ax1.plot(x, fitted_pdf, label=dist_name)
  #print(np.sum(np.power(y - fitted_pdf, 2.0)))
  #print(st.kstest(x, dist_name, args=args))
  #print(st.kstest(x, dist_name, args=params))

ax1.legend(loc='upper right')


ax2 = fig.add_subplot(223)
ax2.title.set_text(f"{split1+1} - {ceil} genes expressed ({len(cols2)} cells)")
ax2.set_ylabel('number of cells')
ax2.set_xlabel('# of genes expressed')
ax2.set_xticks(np.arange(split1, ceil, step=100))
ax2.hist(cols2, bins=175)

ax3 = fig.add_subplot(224)
ax3.title.set_text(f"{ceil+1} - {max_genes} genes expressed ({len(cols3)} cells)")
ax3.set_ylabel('number of cells')
ax3.set_xlabel('# of genes expressed')
ax3.set_xticks(np.arange(ceil, max_genes, step=500))
ax3.hist(cols3, bins='auto')

plt.show()
