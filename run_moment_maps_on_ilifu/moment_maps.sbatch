#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --exclude=slwrk-[102,103,109,125,126]
#SBATCH --cpus-per-task 1
#SBATCH -p Main
#SBATCH --mem=100GB
#SBATCH -J moment_maps
#SBATCH -o moment_maps.%j.stdout
#SBATCH -e moment_maps.%j.stderr
#SBATCH --time=1:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

singularity exec /idia/software/containers/python-3.6.img python Tutorial_moment_maps_sbatch.py

cp variables.py variables.py-copy
mv variables.py-copy variables.${SLURM_JOB_NAME}-${SLURM_JOB_ID}.py