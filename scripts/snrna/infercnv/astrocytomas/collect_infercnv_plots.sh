#!/bin/bash

# Source and destination directories
SRC_DIR="/vast/palmer/pi/verhaak/kcj28/care_mut/results/infercnv/astrocytoma_samples"
DEST_DIR="/vast/palmer/pi/verhaak/kcj28/care_mut/results/infercnv/astrocytoma_plots"

# Create destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Loop over each sample subdirectory
for sample_dir in "$SRC_DIR"/*/; do
    sample_name=$(basename "$sample_dir")
    src_img="${sample_dir}infercnv.png"

    if [[ -f "$src_img" ]]; then
        cp "$src_img" "${DEST_DIR}/${sample_name}_infercnv.png"
        echo "Copied: ${sample_name}_infercnv.png"
    else
        echo "WARNING: No infercnv.png found in ${sample_name}"
    fi
done

echo "Done. Files copied to $DEST_DIR"
