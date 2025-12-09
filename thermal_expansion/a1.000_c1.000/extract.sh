#!/bin/bash
# --- SLURM Directives ---
#SBATCH --job-name=extract                # job name
#SBATCH --mail-user=qhossain@ufl.edu      # Email for job notifications
#SBATCH --mail-type=FAIL,END              # Send email on job failure or completion
#SBATCH --output=out.%j                   # Standard output file (with Job ID)
#SBATCH --error=err.%j                    # Standard error file (with Job ID)
#SBATCH --ntasks=1                        # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=1                 # Number of CPU cores per task
#SBATCH --mem-per-cpu=100G                  # Memory per CPU core
#SBATCH --time=10:00:00                   # Maximum job run time
#SBATCH --partition=hpg-turin


extract_forceconstants -rc2 5 --readforcemap
