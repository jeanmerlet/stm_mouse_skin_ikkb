import pandas as pd
import scprep
import os
import re
import glob

root_mtx_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/all_aligned'
meta_path = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/disease_condition.tsv'
out_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/prrx1_matrices'

conditions = pd.read_csv(meta_path, sep='\t', index_col=0, header=0)

def get_non_psoriasis_sample_dirs(dir_path):
    acc_dirs = glob.glob(dir_path + '/**/Gene/filtered')
    mtx_paths = []
    sample_ids = []
    for acc_dir in acc_dirs:
        acc_id = re.search('aligned_\d+/(.*)_S1', acc_dir).groups()[0]
        if conditions.loc[acc_id, 'condition'] != 'psoriasis':
            mtx_paths.append(os.path.join(dir_path, acc_id + '_S1_L001_R1_001._Solo.out', 'Gene', 'filtered'))
            sample_ids.append(acc_id)
    return mtx_paths, sample_ids

def combine_samples(sample_dirs, sample_ids):
    samples = []
    for i, d in enumerate(sample_dirs):
        print(i)
        samples.append(scprep.io.load_10X(d, gene_labels='both'))
    mtx, labels = scprep.utils.combine_batches(data=samples, batch_labels=sample_ids, append_to_cell_names=True)
    del(samples)
    return mtx, labels

experiments = os.listdir(root_mtx_dir)
experiments.sort()
print(experiments)
for i, experiment in enumerate([os.path.join(root_mtx_dir, exp) for exp in experiments]):
    if 'aligned_4820' in experiment:
        print(f'combining experiment {experiment}')
        sample_dirs, sample_ids = get_non_psoriasis_sample_dirs(experiment)
        mtx, labels = combine_samples(sample_dirs, sample_ids)
        prrx1_iloc = ['PRRX1' in x for x in mtx.columns.values]
        prrx1_mtx_idx = mtx.loc[:, prrx1_iloc] > 0
        mtx = mtx.loc[prrx1_mtx_idx.values, :]
        mtx = scprep.filter.remove_rare_genes(mtx, cutoff=0, min_cells=10)
        #mtx = scprep.normalize.library_size_normalize(mtx)
        #mtx = scprep.transform.log(mtx)
        print(mtx.shape)
        out_path = os.path.join(out_dir, experiments[i] + '_unprocessed.tsv')
        mtx = mtx.sparse.to_dense()
        mtx.to_csv(out_path, sep='\t', header=True, index=True)
        print(f'finished writing experiment {experiment}!')
