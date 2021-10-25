import os
import re
import glob
import subprocess

root_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/star_raw_mtx'
target_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/aligned'

mtx_paths = glob.glob(root_dir + '/*barcodes.tsv', recursive=True)

for mtx_path in mtx_paths:
    acc_id = re.search('(SRR.*)_', mtx_path).groups()[0]
    target_path = os.path.join(target_dir, acc_id + '_1._Solo.out', 'Gene', 'filtered', 'barcodes.tsv')
    os.rename(mtx_path, target_path)

