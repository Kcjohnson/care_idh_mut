#!/bin/bash

#SBATCH --job-name=inferCNV-%a
#SBATCH --chdir=/vast/palmer/pi/verhaak/kcj28/care_mut/
#SBATCH --output=/vast/palmer/pi/verhaak/kcj28/care_mut/logs/infercnv/astrocytomas/infercnv_astrocytomas-%a.log
#SBATCH --error=/vast/palmer/pi/verhaak/kcj28/care_mut/logs/infercnv/astrocytomas/infercnv_astrocytomas-%a.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kevin.c.johnson@yale.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --array=1-45
#SBATCH --partition=priority
#SBATCH -A prio_verhaak

# Run infercnv script.
bash /vast/palmer/pi/verhaak/kcj28/care_mut/scripts/infercnv/astrocytomas/run_infercnv_astrocytomas.sh ${SLURM_ARRAY_TASK_ID}

### END ###
