#!/bin/bash 
#SBATCH -J ctloop_array # A single job name for the array 
#SBATCH -n 1 # Number of cores 
#SBATCH -N 1 # All cores on one machine
#SBATCH --mem 4000 # Memory request (4Gb) 
#SBATCH -t 0-2:00 # Maximum execution time (D-HH:MM) 
#SBATCH -o ctloop_array_%A_%a.out # Standard output 
#SBATCH -e ctloop_array_%A_%a.err # Standard error  

module load Python/3.9.6-GCCcore-11.2.0
module load Biopython/1.78-foss-2020b

python ctloop "${SLURM_ARRAY_TASK_ID}"

exit 0