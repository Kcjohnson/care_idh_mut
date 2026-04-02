#!/bin/bash

### This script runs infercnv on each sample in the CARE astrocytoma dataset ###

# Activate the conda environment for running inferCNV, if not currently active.
module load miniconda
conda activate infercnv_env

### Input arguments ###
ARRAYID="`expr $1`"

# Sample ID list.
SAMPLE_ID_FILE="/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/infercnv_inputs/astrocytomas/astrocytoma_dataset_infercnv_samples.txt"

# Target sample ID.
SAMPLE_ID=$(sed "${ARRAYID}q;d" ${SAMPLE_ID_FILE})

# Temp directory (will be created if not existing).
TEMP_DIR=/vast/palmer/scratch/verhaak/kcj28/infercnv/
mkdir -p ${TEMP_DIR}

# Output directory (will be created if not existing).
OUTDIR=/vast/palmer/pi/verhaak/kcj28/care_mut/results/infercnv/astrocytoma_samples/
mkdir -p ${OUTDIR}

### Print out relevant information.
STARTTIME=`date`
echo $STARTTIME
echo "Analyzing $SAMPLE_ID"
echo ""

# Run Rscript.
#ulimit -s unlimited
Rscript --verbose /vast/palmer/pi/verhaak/kcj28/care_mut/scripts/infercnv/astrocytomas/run_infercnv_astrocytomas.R \
			-i ${SAMPLE_ID} \
			-t ${TEMP_DIR} \
			-o ${OUTDIR} \
			-n ${SLURM_CPUS_PER_TASK}


### END ####
