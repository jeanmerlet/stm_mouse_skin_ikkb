#!/bin/bash 
  
#SBATCH -A SYB105
#SBATCH -t 4:00:00
#SBATCH -J comb_mtx
#SBATCH -o ./logs/combine_mtx%J.out
#SBATCH -e ./logs/combine_mtx%J.err
#SBATCH -N 1

module load python
source activate andes_base
which python
conda list

PYSCRIPT="/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/scripts/human_ad/tsne_clustering/reynolds/write_mtxs_to_combined_prrx1.py"

srun -N 1 -n 1 -c 16 python $PYSCRIPT

