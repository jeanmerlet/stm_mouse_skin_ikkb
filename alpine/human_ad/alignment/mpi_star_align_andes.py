from mpi4py import MPI
import subprocess
import os

prefix = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet'
star = os.path.join(prefix, 'rna-seq_tools/alignment/andes/STAR-2.7.9a/source/STAR')
genome = os.path.join(prefix, 'rna-seq_tools/alignment/ref_gens/human/GRCh38/star_index_2.7.9a')
bc_whitelist = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/rna-seq_tools/alignment/bc_wl/737K-august-2016.txt'
data = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/fastq/8090'
out_dir = '/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/data/human_ad/reynolds/aligned_809_reverse/'

def sc_star_align(fastq1, fastq2, prefix):
  out_prefix = out_dir + prefix + '_'
  subprocess.run([star, '--runThreadN', '21', '--genomeDir', genome, '--soloType CB_UMI_Simple', '--soloCBwhitelist', bc_whitelist, '--soloStrand', 'Reverse', '--outFileNamePrefix', out_prefix, '--readFilesIn', fastq2, fastq1, '--soloUMIfiltering', 'MultiGeneUMI', '--readFilesCommand', 'zcat', '--soloBarcodeReadLength', '0'])

def bulk_star_align_single(fastq, prefix):
  out_prefix = out_dir + prefix + '_'
  subprocess.run([star, '--runThreadN', '21', '--genomeDir', genome, '--outFileNamePrefix', out_prefix, '--readFilesIn', fastq])

def bulk_star_align_paired(fastq1, fastq2, prefix):
  out_prefix = out_dir + prefix + '_'
  subprocess.run([star, '--runThreadN', '21', '--genomeDir', genome, '--outFileNamePrefix', out_prefix, '--readFilesIn', fastq1, fastq2])

fastq_list = []
for r, d, f in os.walk(data):
  for fastq in f:
    if '.fastq' in fastq:
      fastq_list.append(os.path.join(r, fastq))
fastq_list.sort()

paired_fastq_list = []
pair = []
for i, fastq in enumerate(fastq_list):
  pair.append(fastq)
  if i % 2 == 1:
    paired_fastq_list.append(pair)
    pair = []

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

for i, pair in enumerate(paired_fastq_list):
  if i % size != rank: continue
  fastq_path1 = pair[0]
  fastq_path2 = pair[1]
  head, fastq = os.path.split(fastq_path1)
  prefix = fastq[:-8]
  print(f'{prefix} (task number {i}) being aligned by {rank} of {size}')
  sc_star_align(fastq_path1, fastq_path2, prefix)
