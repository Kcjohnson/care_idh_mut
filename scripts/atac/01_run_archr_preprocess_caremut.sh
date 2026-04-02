#!/bin/bash
#SBATCH --job-name=caremut-archr
#SBATCH --chdir=/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac
#SBATCH --output=/vast/palmer/pi/verhaak/kcj28/care_idh_mut/logs/archr/caremut-create-arrow.log
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=kevin.c.johnson@yale.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=23:00:00
#SBATCH --partition=priority
#SBATCH -A prio_verhaak

# Load conda env to perform ArchR analysis of snATACseq data.
module load miniconda
conda activate ARCHRenv

STARTTIME=`date`
echo $STARTTIME
echo "conda env is ${CONDA_PREFIX}"
echo ""

# Run Rscript
Rscript --verbose /vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/atac/01_run_archr_preprocess_caremut.R


### END ####

