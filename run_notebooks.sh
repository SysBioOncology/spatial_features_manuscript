#!/usr/bin/env bash
#SBATCH -J run_notebooks.R
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out


source "${CONDA_PREFIX}/bin/activate" "r-4.3.3"

Rscript "${PWD}/analysis/run_notebooks.R" 
