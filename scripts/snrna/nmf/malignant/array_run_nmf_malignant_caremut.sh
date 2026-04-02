#!/bin/bash

#SBATCH --job-name=nmf-%a
#SBATCH --chdir=/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/nmf_2026/
#SBATCH --output=/vast/palmer/pi/verhaak/kcj28/care_mut/logs/nmf_2026/malignant/nmf_malignant_n10_2026-%a.log
#SBATCH --error=/vast/palmer/pi/verhaak/kcj28/care_mut/logs/nmf_2026/malignant/nmf_malignant_n10_2026-%a.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kevin.c.johnson@yale.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=10-23:59:59
#SBATCH --partition=priority
#SBATCH -A prio_verhaak
#SBATCH --array=1-74


# Run NMF script.
bash /vast/palmer/pi/verhaak/kcj28/care_mut/scripts/nmf_2026/malignant/run_nmf_malignant_caremut.sh ${SLURM_ARRAY_TASK_ID}

### END ###