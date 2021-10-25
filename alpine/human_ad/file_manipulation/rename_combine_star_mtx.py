import os
import re
import glob
import subprocess

root_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/aligned'
target_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/he/star_raw_mtx'

mtx_paths = glob.glob(root_dir + '/**/filtered/barcodes.tsv', recursive=True)

for mtx_path in mtx_paths:
    acc_id = re.search('(SRR.*)_1.', mtx_path).groups()[0]
    head, tail = os.path.split(mtx_path)
    os.rename(mtx_path, os.path.join(target_dir, acc_id + '_' + tail))

