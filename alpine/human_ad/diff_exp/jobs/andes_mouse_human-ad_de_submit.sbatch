#!/bin/bash

#SBATCH -N 1
#SBATCH -A syb105
#SBATCH -J deseq2
#SBATCH -t 8:00:00
#SBATCH -e ./logs/deseq2_%J.err
#SBATCH -o ./logs/deseq2_%J.out

srun -n 1 -c 8 Rscript /gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/scripts/human_ad/diff_exp/human_ad_de_reynolds.R
