#!/bin/bash
#SBATCH -A anirudh
#SBATCH -p long
#SBATCH -c 4
#SBATCH --mem=10240
#SBATCH -t 25:00:00

/scratch/matlab/R2013b/bin/matlab -nodesktop -nosplash -r run_sim_par
