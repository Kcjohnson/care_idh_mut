#!/bin/bash

#SBATCH --job-name=rerun_NL26_infercnv
#SBATCH --chdir=/vast/palmer/pi/verhaak/kcj28/care_mut/
#SBATCH --output=/vast/palmer/pi/verhaak/kcj28/care_mut/logs/infercnv/astrocytomas/rerun_NL26_1_infercnv.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kevin.c.johnson@yale.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=23:00:00
#SBATCH --partition=priority
#SBATCH -A prio_verhaak

# Activate the conda environment for running inferCNV, if not currently active.
module load miniconda
conda activate infercnv_env

# Run infercnv script.
Rscript --verbose /vast/palmer/pi/verhaak/kcj28/care_mut/scripts/infercnv/astrocytomas/run_infercnv_NL26_1.R

### END ###
