#!/bin/bash 

#SBATCH -A SYB105
#SBATCH -J combine_mtx
#SBATCH -o ./logs/combine_mtx.%J.out
#SBATCH -e ./logs/combine_mtx.%J.err
#SBATCH -t 8:00:00
#SBATCH -N 1

source activate andes_base

PYSCRIPT="/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/scripts/human_ad/tsne_clustering/reynolds/write_combined_processed_mtx.py"

srun -N 1 -n 1 -c 16 python $PYSCRIPT
