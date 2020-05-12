#!/bin/bash
#SBATCH -A anirudh
#SBATCH -p long
#SBATCH -c 5
#SBATCH --mem=30720
#SBATCH -t 25:00:00

/scratch/matlab/R2013b/bin/matlab -nodesktop -nosplash -r run_optimal_G

