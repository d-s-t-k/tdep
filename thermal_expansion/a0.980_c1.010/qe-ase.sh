#!/bin/bash
# --- SLURM Directives ---
#SBATCH --job-name=ASE-QE
#SBATCH --mail-user=qhossain@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --output=%A_%a.out 
#SBATCH --error=%A_%a.err
#SBATCH --partition=hpg-b200
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:10:00

module purge
module load nvhpc/25.3 openmpi/5.0.7 espresso/7.5.0

# ---- keep UCX/UD out entirely; use shm+tcp only ----
export OMPI_MCA_pml=ob1
export OMPI_MCA_pml_base_exclude=ucx
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_coll=^hcoll
export OMPI_MCA_mtl=^ofi,psm2,psm,portals4
unset UCX_TLS UCX_NET_DEVICES

# threads (GPU QE usually prefers 1)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

# One config per array task:
python qe-ase.py
