import csv
import gzip
import os
import scipy.io
import scipy.interpolate as interp
import scipy.optimize as optim
import numpy as np
import matplotlib.pyplot as plt

experiment = 'FGC2091_7_92848'
matrix_dir = '/Users/6j9/projects/mouse/data/' + experiment + '/raw_feature_bc_matrix'

# each element of m represents feature UMI count (rows) / unique barcode (cols)
m = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx"))

# sum barcode umi counts
umi_counts = np.asarray(m.sum(axis=0))[0]
#print(np.max(umi_counts))

# remove barcodes with 0 UMI count for all features from m
nonzero_umi_idx = np.nonzero(umi_counts != 0)
nonzero_umi_counts = umi_counts[nonzero_umi_idx]

# calculate ambient rna profile
#ambient_cutoff = 100
#cutoff_umi_idx = np.nonzero(umi_counts != 0 & umi_counts <= ambient_cutoff)
#print(len(cutoff_umi_idx))
#ambient_x = np.arange(1, len(ambient_y) + 1, 1)

# order UMI counts by decreasing total count
y = np.flip(np.sort(nonzero_umi_counts))
#print(np.shape(y))
#print(np.shape(y[np.nonzero(y > 1300)]))
test = y[85]
#print(test)
x = np.arange(1, len(y) + 1, 1)

# fit a spline around the knee bend range
# spline is cubic so as to ensure double differentiability
#def get_spline(startx, endx, y):
#  y = y[startx:endx]
#  x = np.arange(startx, endx, 1)
#  spline_f = interp.CubicSpline(x, y)
#  return spline_f
# cutoff point is signed curvature of spline
#def signed_curvature(x):
#  return (splinef2(x) / ((1+splinef1(x)**2)**1.5))
#scurv = lambda x: splinef2(x) / (1 + splinef1(x)**2)**1.5
#start_x, end_x = 100, 10000
#spline_f = get_spline(start_x, end_x, y)
#splinef1 = spline_f.derivative(nu=1)
#splinef2 = spline_f.derivative(nu=2)
#spline_x = np.arange(start_x, end_x, 1)
#spline_y = spline_f(spline_x)
#sc = signed_curvature(spline_x)
#knee_point = spline_y[np.argmin(sc)]
#print(knee_point)
#print(spline_y[np.argmin(scurv(spline_x))])
#spline_y2 = splinef2(spline_x)
#inflection_points_idx = np.nonzero(spline_y2 == 0)
#inflection_point = spline_y[np.min(inflection_points_idx)]

fig = plt.figure()
fig.suptitle('92848 - filtered')
ax = fig.add_subplot(111)
ax.title.set_text(f"Barcode UMI Count by Decreasing Rank")
ax.set_xlabel('Rank')
ax.set_ylabel('UMI Count')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1, 500000)
ax.set_ylim(1, 50000)
ax.plot(y, marker='o', markeredgecolor='black', markerfacecolor='none', linestyle='')
#ax.plot(spline_x, spline_y, label='spline')
ax.axhline(500, color='C2', label='UMI count = 500')
ax.axhline(test/10, color='C3', label='99th pctl / 10 = 1300')
ax.legend(loc='upper right')
plt.show()
