#!/bin/bash

#SBATCH -N 1
#SBATCH -A syb105
#SBATCH -J deseq2
#SBATCH -t 4:00:00
#SBATCH -e ./logs/deseq2_%J.err
#SBATCH -o ./logs/deseq2_%J.out

Rscript /gpfs/alpine/syb105/proj-shared/Personal/jmerlet/rna-seq_tools/de/deseq2/mouse_de_cre-pos_v_cre-neg.R
