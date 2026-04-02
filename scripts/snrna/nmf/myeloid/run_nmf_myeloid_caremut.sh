#!/bin/bash

### This script runs NMF on each sample's myeloid cells in the CAREmut dataset ###

# Activate the conda environment for running NMF, if not currently active.
module load miniconda
conda activate NMFenv


### Input arguments ###
ARRAYID="`expr $1`"

# Identify where the expression matrices are stored
IN_PATH=/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/nmf_2026/Myeloid/

# Output directory (will be created if not existing).
OUT_PATH=/vast/palmer/pi/verhaak/kcj28/care_mut/results/nmf_res_2026/myeloid_n73/
mkdir -p $OUT_PATH

# Set other variables for NMF
RANK_LB=3
RANK_UB=10
NRUN=10

# Print out relevant information.
STARTTIME=`date`
echo $STARTTIME
echo "Analyzing indexed sample: $ARRAYID"
echo ""


# Run Rscript.
Rscript --verbose /vast/palmer/pi/verhaak/kcj28/care_mut/scripts/nmf_2026/myeloid/run_nmf_myeloid_caremut.R $IN_PATH $ARRAYID $OUT_PATH $RANK_LB $RANK_UB $NRUN

ENDTIME=`date`
echo $ENDTIME

### END ####