#!/bin/bash 

#SBATCH -A SYB105
#SBATCH -J para_STAR_align
#SBATCH -o ./logs/align.%J.out
#SBATCH -e ./logs/align.%J.err
#SBATCH -t 24:00:00
#SBATCH -N 52

MPIFILE="/gpfs/alpine/syb105/proj-shared/Personal/jmerlet/projects/mouse/scripts/human_ad/alignment/mpi_star_align_andes.py"

srun -n 52 -c 16 python $MPIFILE
