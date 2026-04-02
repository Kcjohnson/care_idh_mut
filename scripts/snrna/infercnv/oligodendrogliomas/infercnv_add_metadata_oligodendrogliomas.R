##################################
# Add metadata for infercnv results for CARE oligodendroglioma cohort.
# Author: Kevin Johnson
# Date Updated: 2026.03.18
##################################

## Run in infercnv_env conda environment

## Extract infercnv metadata to identify malignant cells.
library(tidyverse)
library(infercnv) # infercnv_1.14.2

## Specify directories:
parent_outdir <- "/vast/palmer/pi/verhaak/kcj28/care_mut/results/infercnv/oligodendroglioma_samples"

## List all samples in each directory.
samples_dir <- dir(parent_outdir, pattern = "infercnv_")

# Create seurat object that adds metadata for these infercnv runs.
seurat_list_add_meta <- lapply(samples_dir, function(sample_name){
  print(sample_name)
  tryCatch(
    {
      seurat_obj <- infercnv::add_to_seurat(infercnv_output_path = paste0(parent_outdir, "/", sample_name))
      return(seurat_obj) # Return the result if successful
    },
    error = function(e) {
      # Print the error message
      cat("Error occurred for", sample_name, ":")
      message(e)
      return(NULL) # Return NULL to indicate failure
    }
  )
})

# Note: Since no seurat object is actually provided, this will write a metadata matrix for each sample folder: map_metadata_from_infercnv.txt

### END ###