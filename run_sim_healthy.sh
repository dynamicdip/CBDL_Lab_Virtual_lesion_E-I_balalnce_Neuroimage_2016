#!/bin/bash
#SBATCH -A anirudh
#SBATCH -p long
#SBATCH -c 2
#SBATCH --mem=25600
#SBATCH -t 25:00:00

/scratch/matlab/R2013b/bin/matlab -nodesktop -nosplash -r run_sim
