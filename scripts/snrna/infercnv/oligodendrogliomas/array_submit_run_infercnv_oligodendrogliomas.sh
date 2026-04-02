#!/bin/bash

#SBATCH --job-name=oligodendroglioma_infercnv-%a
#SBATCH --chdir=/vast/palmer/pi/verhaak/kcj28/care_mut/
#SBATCH --output=/vast/palmer/pi/verhaak/kcj28/care_mut/logs/infercnv/oligodendrogliomas/infercnv_oligodendrogliomas-%a.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kevin.c.johnson@yale.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=128G
#SBATCH --time=23:59:00
#SBATCH --array=1-30
#SBATCH --partition=priority
#SBATCH -A prio_verhaak


# Run infercnv script.
bash /vast/palmer/pi/verhaak/kcj28/care_mut/scripts/infercnv/oligodendrogliomas/run_infercnv_oligodendrogliomas.sh ${SLURM_ARRAY_TASK_ID}

### END ###
