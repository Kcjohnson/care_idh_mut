##################################
# Derive IDH-mutant hierarchy coordinates for IDH-mutant malignant cells
# Author: Kevin Johnson
# Date Updated: 2026.04.10
##################################

library(tidyverse)
library(RColorBrewer)
library(viridis)
library(Seurat)
library(Matrix)
library(openxlsx)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
out_data_dir <- file.path(proj_dir, "processed_data/rna/")
script_dir  <- file.path(proj_dir, "scripts")

setwd(proj_dir)

# Avishay Spitzer provided these helper scripts. Use the score_within_samples*() functions for malignant cells.
source(file.path(script_dir, "utils", "plot_theme.R"))
source(file.path(script_dir, "utils", "caremut_utils.R"))


### ### ### ### ### ### ### ###
# Set-up
### ### ### ### ### ### ### ###
# Load in the CARE IDH-mutant metadata and count matrices.
md <- read.table(paste0(proj_dir, "/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20260320.txt"), sep = "\t", row.names = 1, header = TRUE)
md_trim <- md %>% 
  dplyr::select(SampleID, case_barcode, idh_codel_subtype, care_id) %>% 
  distinct()
rownames(md_trim) <- NULL

# Load in the CAREmut UMI produced at the beginning of the project.
umi_data_all <- readRDS("data/snrna/care_mut_umi_data_all_20230729.RDS")
names(umi_data_all)


# Use top 50 genes from each signature, when possible, to be consistent with features used

# Public IDH-A and IDH-O signatures from Tirosh Nature 2016 and Venteicher 2017 Science papers (SmartSeq2 single cell data)
# Venteicher Astrocytoma
venteicher <- readWorkbook("data/misc/venteicher_table_s3.xlsx", startRow = 5, colNames = TRUE)
venteicher_signatures <- venteicher %>% 
  dplyr::select("Venteicher_OC_scRNA"=`Oligo-program.(Fig..2C)`, "Venteicher_AC_scRNA"=`Astro-program.(Fig..2C)`,
                "Venteicher_Stemness_scRNA"=`Stemness.program.(Fig..3C)`)
venteicher_signatures$Venteicher_OC_scRNA <- trimws(venteicher_signatures$Venteicher_OC_scRNA)
venteicher_signatures$Venteicher_AC_scRNA <- trimws(venteicher_signatures$Venteicher_AC_scRNA)
venteicher_signatures$Venteicher_Stemness_scRNA <- trimws(venteicher_signatures$Venteicher_Stemness_scRNA)
venteicher_sig_list <- as.list(venteicher_signatures[1:50,1:3])
venteicher_sig_list <- lapply(venteicher_sig_list, function(x) x[!is.na(x)])

# Tirosh Oligodendroglioma
tirosh <- readWorkbook("data/misc/tirosh_nature_2016_supplementary_table_1.xlsx", startRow = 9, colNames = TRUE)
tirosh_signatures <- tirosh %>% 
  dplyr::select("Tirosh_OC_scRNA"=`OC.(PCA-only)`,
                "Tirosh_AC_scRNA"=`AC.(PCA-only)`,
                "Tirosh_Stemness_scRNA"=`stemness`)
tirosh_signatures$Tirosh_OC_scRNA <- trimws(tirosh_signatures$Tirosh_OC_scRNA)
tirosh_signatures$Tirosh_AC_scRNA <- trimws(tirosh_signatures$Tirosh_AC_scRNA)
tirosh_signatures$Tirosh_Stemness_scRNA <- trimws(tirosh_signatures$Tirosh_Stemness_scRNA)
tirosh_signatures_list <- as.list(tirosh_signatures[1:50,1:3])

sigs_list <- c(venteicher_sig_list, tirosh_signatures_list)

# Extract the malignant cells. I am removing SJ02-3, which had only a few malignant cells post-QC
md_malignant <- md %>% 
  filter(CellType_final=="Malignant", SampleID!="SJ02-3") 

# Define a function that subsets columns based on CellID
subset_columns <- function(mat, cell_ids) {
  col_idx <- which(colnames(mat) %in% cell_ids)
  mat[, col_idx, drop = FALSE]
}

# Apply the function to each element of the list
umi_data_all_malignant <- map(umi_data_all, subset_columns, cell_ids = md_malignant$CellID)
names(umi_data_all_malignant)

### ### ### ### ### ### ### ###
# IDH-mutant hierarchy coordinates
### ### ### ### ### ### ### ###
# This is a separate scoring method than what was used for the malignant metaprograms
set.seed(123)
mp_hierarchy_scores <- score_within_samples_log2_hierarchy(umi_data_all_malignant, md = md_malignant, sigs = sigs_list)

# Compute the lineage score
mp_hierarchy_scores$Lineage <- apply(mp_hierarchy_scores[, c("Tirosh_AC_scRNA", "Tirosh_OC_scRNA")], 1, max)

# Compute the stemness score
mp_hierarchy_scores$Stemness <- apply(mp_hierarchy_scores[, c("Tirosh_Stemness_scRNA", "Lineage")], 1, function(x) x[1] - x[2])

# Compute the lineage score which will be used for plotting (i.e. with the random component for the cells with negative lineage scores)
mp_hierarchy_scores$LineagePlot <- apply(mp_hierarchy_scores[, c("Tirosh_AC_scRNA", "Tirosh_OC_scRNA")], 1, function(x) {
  res <- max(x[1], x[2])
  if (res < 0)
    res <- runif(1, min = 0, max = 0.15)
  else if(x[1] > x[2])
    res <- -1 * res
  return (res)
})

# Write out the x-y coordinates for this plot
write.table(mp_hierarchy_scores, file = paste0(proj_dir, "/results/scoring/malignant_malignant_hierarchy_signature_scores.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
