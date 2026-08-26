#!/bin/bash
#SBATCH -p debug
#SBATCH -J tsr_v2_inst
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH -o v2/logs/install_%j.out
#SBATCH -e v2/logs/install_%j.err
set -euo pipefail
module purge
module load apps/R/4.3.0-gcc-11.3.1
export R_LIBS_USER=/arf/scratch/hbulut/TSRClust/Rlibs_4.3.0
cd /arf/scratch/hbulut/TSRClust/TSRClust_SciReports_TRUBA
mkdir -p "$R_LIBS_USER" v2/logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript v2/scripts/install_v2_packages.R
echo V2_INSTALL_OK
