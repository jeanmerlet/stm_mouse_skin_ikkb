#!/bin/bash 

#SBATCH -A SYB105
#SBATCH -J batch
#SBATCH -o ./logs/batch.%J.out
#SBATCH -e ./logs/batch.%J.err
#SBATCH -t 8:00:00
#SBATCH -N 1

source activate andes_base

RSCRIPT="/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/scripts/human_ad/tsne_clustering/reynolds/seurat_cluster_combined_matrix.r"

srun -N 1 -n 1 -c 16 Rscript $RSCRIPT
