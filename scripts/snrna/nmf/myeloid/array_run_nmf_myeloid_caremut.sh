#!/bin/bash

#SBATCH --job-name=nmf-myeloid-%a
#SBATCH --chdir=/vast/palmer/scratch/verhaak/kcj28/nmf_myeloid
#SBATCH --output=/vast/palmer/pi/verhaak/kcj28/care_mut/logs/nmf_2026/myeloid/nmf_2026_myeloid_caremut-%a.log
#SBATCH --error=/vast/palmer/pi/verhaak/kcj28/care_mut/logs/nmf_2026/myeloid/nmf_2026_myeloid_caremut-%a.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kevin.c.johnson@yale.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=6-00:00:00
#SBATCH --partition=priority
#SBATCH -A prio_verhaak
#SBATCH --array=1-73


# Run NMF script.
bash /vast/palmer/pi/verhaak/kcj28/care_mut/scripts/nmf_2026/myeloid/run_nmf_myeloid_caremut.sh ${SLURM_ARRAY_TASK_ID}

### END ###