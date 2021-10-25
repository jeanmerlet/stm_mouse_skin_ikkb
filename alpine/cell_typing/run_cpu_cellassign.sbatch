#!/bin/bash

#SBATCH -A SYB106
#SBATCH -t 24:00:00
#SBATCH -J cellassign_cpu
#SBATCH -o ./logs/mouse_ct.%J.o
#SBATCH -e ./logs/mouse_ct.%J.e
#SBATCH -N 1
#SBATCH --mem=0

celltyping_path='/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/scripts/cell_typing/mouse_cellassign.R'

Rscript $celltyping_path
