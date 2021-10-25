#!/bin/bash

#SBATCH -N 1
#SBATCH -A syb105
#SBATCH -J deseq2
#SBATCH -t 24:00:00

Rscript \
/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/rna-seq_tools/de/deseq2/mouse_de.R
