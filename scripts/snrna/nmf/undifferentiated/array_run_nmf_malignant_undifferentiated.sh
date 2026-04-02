#!/bin/bash

#SBATCH --job-name=nmf-%a
#SBATCH --chdir=/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/nmf_res_caremut/undifferentiated
#SBATCH --output=/vast/palmer/pi/verhaak/kcj28/care_idh_mut/logs/nmf/undifferentiated/undifferentiated_nmf_%a.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kevin.c.johnson@yale.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=8-23:59:59
#SBATCH --partition=priority
#SBATCH -A prio_verhaak
#SBATCH --array=1-72


# Run NMF script.
bash /vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/snrna/nmf/undifferentiated/run_nmf_malignant_undifferentiated.sh ${SLURM_ARRAY_TASK_ID}

### END ###