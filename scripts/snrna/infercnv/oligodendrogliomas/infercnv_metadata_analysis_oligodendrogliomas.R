##################################
# Analyze infercnv metadata to revise uncertain malignant classification for CARE oligodendroglioma samples
# Author: Kevin Johnson
# Date: 2026.03.20
##################################

library(tidyverse)
library(purrr)
library(ggdist)
library(viridis)
library(Seurat)
library(harmony)
library(Matrix)
library(ggpubr)

# Define directories
infercnv_data_dir <- "/vast/palmer/pi/verhaak/kcj28/care_mut/results/infercnv/oligodendroglioma_samples/"
out_data_dir <- "/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna"
fig_dir <- "/vast/palmer/pi/verhaak/kcj28/care_mut/results/figures/rna/"
source("/vast/palmer/pi/verhaak/kcj28/care_mut/scripts/misc/plot_theme.R")

################################################
## Step 1: Load the results from Seurat pre-processing:
################################################
# This represents the metadata as processed through Seurat for removing singlets, harmony batch correction based on case barcode, and annotating clusters based on previously determined marker genes.
caremut_md <- read.table("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/oligodendroglioma_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_20260318.txt", sep = "\t", header = TRUE)

colnames(caremut_md)[colnames(caremut_md)=="scna_burden"] <- "scna_burden_old"
colnames(caremut_md)[colnames(caremut_md)=="malignant_cor"] <- "malignant_cor_old"
caremut_md$RNA_snn_res.0.1 <- NULL
caremut_md$seurat_clusters <- NULL
all(caremut_md$sample_barcode==caremut_md$sample_barcode.y, na.rm = TRUE)
caremut_md$sample_barcode.y <- NULL
caremut_md$sample_barcode.x <- NULL


################################################
## Step 2: Load the results from infercnv run
################################################

# Get a list of file paths for files specific to the CAREmut
file_paths <- dir(infercnv_data_dir, full.names = TRUE)

# Initialize an empty list to store the data frames
data_list <- list()

# Define cohort-wide malignant cells
malignant_rows <- caremut_md %>% filter(CellType_clusters == "Malignant")
malignant_cells <- malignant_rows$CellID

# Extract last element from a strsplit.
last <- function(x) { return( x[length(x)] ) }

# Iterate over each file
for (file_path in file_paths) {
  # Get the directory/sample name
  dir_name <- sapply(strsplit(file_path, "/"), last)
  
  # Attempt to read the file, handling missing file errors due to infercnv not completing
  tryCatch(
    {
      # Read the file
      data <- read.table(paste0(file_path, "/map_metadata_from_infercnv.txt"), header = TRUE, sep = "\t")
      
      # Add a column with the directory name
      print(dir_name)
      data$dir_name <- dir_name
      
      # Calculate each cell's CNV correlation with the average tumor cell. NOTE: That this is looking at the proportion of a chromosome with amp/deletion - not specific segments.
      malignant_data <- data[rownames(data)%in%malignant_cells, ]
      malignant_data_filt <- malignant_data %>% 
        dplyr::select(starts_with("proportion_loss_"), starts_with("proportion_dupli_"))
      
      # Assumes the matrix is named 'malignant_data_filt'
      column_averages <- colMeans(malignant_data_filt)
      
      # Restricting to the features that are proportional loss or duplicated.
      data_filt <- data  %>% 
        dplyr::select(starts_with("proportion_loss_"), starts_with("proportion_dupli_"))
      
      # Compute correlations between each row of 'malignant_data_filt' and 'column_averages'. Spearman correlation.
      row_correlations <- apply(data_filt, 1, function(row) cor(row, column_averages, method = "s"))
      all(rownames(data)%in%names(row_correlations))
      
      data$malignant_cor <- row_correlations
      
      # Append the data frame to the list
      data_list[[dir_name]] <- data
    },
    error = function(e) {
      # Print an error message for the missing file
      cat("Error occurred while processing", file_path, ":")
      message(e)
    }
  )
}

# warnings() - In cor(row, column_averages, method = "s") : the standard deviation is zero

# Combine all data frames into a single data frame
combined_data <- do.call(rbind, data_list)
combined_data$CellID <- sapply(strsplit(rownames(combined_data), "\\."), last)
rownames(combined_data) <- sapply(strsplit(rownames(combined_data), "\\."), last)

# Retrieve the continuous scores. For example, "proportion_dupli_chr1" is how many of the genes kept after filtering on chr1 are part of a duplication. 
infercnv_continuous <- combined_data %>% 
  dplyr::select(CellID, sample_id = dir_name, malignant_cor, starts_with("proportion_loss_"), starts_with("proportion_dupli_")) %>% 
  inner_join(caremut_md, by="CellID")

# There will be some cells that were not considered for infercnv because their celltype. 49 cells where the cells were too few to be run by celltype in infercnv.
dim(caremut_md)[1]-dim(infercnv_continuous)[1]

# Define a somatic copy number alteration (scna) burden metric across these cells by taking the row average:
infercnv_merge <- combined_data %>% 
  dplyr::select(CellID, sample_id = dir_name, malignant_cor, starts_with("proportion_loss_"), starts_with("proportion_dupli_"))

# Make sure that these are the correct fields:
head(infercnv_merge[,4:47])
colnames(head(infercnv_merge[,4:47]))
infercnv_merge$scna_burden <-  rowMeans(infercnv_merge[,4:47])

# Add Azimuth normal brain cell type annotation:
oligodendroglioma_azimuth_md <- read.table("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/oligodendroglioma_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_azimuth_20260320.txt", header = TRUE, sep = "\t", row.names = 1)

# Load back in the Seurat object with these IDH-mutant oligodendroglioma data:
oligo_sobj_all_harmony_cleaned <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/oligodendroglioma_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_20260318.RDS")
oligo_sobj_all_harmony_cleaned@meta.data$scna_burden <- NULL
oligo_sobj_all_harmony_cleaned@meta.data$malignant_cor <- NULL
oligo_sobj_all_harmony_cleaned@meta.data$predicted.class <- NULL
oligo_sobj_all_harmony_cleaned@meta.data$predicted.subclass <- NULL


# Tidy up data to be merged. These will reflect the 49 cells with missing data.
infercnv_merge_comb <- caremut_md %>% 
  dplyr::select(-predicted.class, -predicted.subclass) %>% 
  inner_join(oligodendroglioma_azimuth_md, by="CellID") %>% 
  left_join(infercnv_merge, by="CellID") 

variables_to_keep <- colnames(infercnv_merge_comb)[!colnames(infercnv_merge_comb)%in%colnames(oligo_sobj_all_harmony_cleaned@meta.data)]
infercnv_merge_to_add <- infercnv_merge_comb %>% 
  dplyr::select(CellID, all_of(variables_to_keep))

# Add the metadata of the infercnv results
seurat_md <- oligo_sobj_all_harmony_cleaned@meta.data
add_data <- seurat_md %>% 
  left_join(infercnv_merge_to_add, by="CellID")
row.names(add_data) <- row.names(oligo_sobj_all_harmony_cleaned[[]])

# Sanity check to ensure that the data align.
all(row.names(add_data)==add_data$CellID)
oligo_sobj_all_harmony_cleaned <- AddMetaData(oligo_sobj_all_harmony_cleaned, metadata = add_data)

# Confirmed that all relevant data is there. This should be FALSE.
any(is.na(oligo_sobj_all_harmony_cleaned@meta.data$predicted.class))

# Define a co-deletion signal using the proportion of loss for both Chr1 and Chr19.
oligo_sobj_all_harmony_cleaned$proportion_loss_chr1chr19 = oligo_sobj_all_harmony_cleaned$proportion_loss_chr1+oligo_sobj_all_harmony_cleaned$proportion_loss_chr19

# The deletion signal is stronger for Chr1 presumably because it is a larger section of the genome
FeaturePlot(oligo_sobj_all_harmony_cleaned, features = "proportion_loss_chr1", pt.size = 1)
FeaturePlot(oligo_sobj_all_harmony_cleaned, features = "proportion_loss_chr19", pt.size = 1)
FeaturePlot(oligo_sobj_all_harmony_cleaned, features = "proportion_loss_chr1chr19", pt.size = 1)
FeaturePlot(oligo_sobj_all_harmony_cleaned, features = "scna_burden", pt.size = 1)
FeaturePlot(oligo_sobj_all_harmony_cleaned, features = "malignant_cor", pt.size = 1)

# The malignant signature is strongest for Chr1 alone. The excitatory neurons because they express more genes also have a signal
VlnPlot(oligo_sobj_all_harmony_cleaned, group.by = "CellType_clusters", features = "proportion_loss_chr1")
VlnPlot(oligo_sobj_all_harmony_cleaned, group.by = "CellType_clusters", features = "proportion_loss_chr19")
VlnPlot(oligo_sobj_all_harmony_cleaned, group.by = "CellType_clusters", features = "proportion_loss_chr1chr19")

oligo_infercnv_res <- oligo_sobj_all_harmony_cleaned@meta.data

# Investigate a reasonable cutoff for malignant definition based on chromosome 1 deletion.
oligo_infercnv_res_malignant <- oligo_infercnv_res %>% 
  filter(CellType_clusters=="Malignant")

# There are 1,045 cells (~1.5%) that are labelled as malignant that fall below 0.15 for both chromosomes
sum(oligo_infercnv_res_malignant$proportion_loss_chr1<0.15 & oligo_infercnv_res_malignant$proportion_loss_chr19<0.15)

# Low-quality calls across **MALIGNANT** nuclei are defined on the absence of these infercnv signals. 
# SJ02-3 and NL01-2 have the highest number of cells.
oligo_infercnv_res_malignant_lq <- oligo_infercnv_res_malignant %>% 
  filter(CellType_clusters=="Malignant", proportion_loss_chr1 < 0.15, proportion_loss_chr19<0.15)

table(oligo_infercnv_res_malignant_lq$SampleID) 

# Define the low-quality **NON-malignant** cell states based on high infercnv signal.
sum(is.na(oligo_infercnv_res$predicted.class))
oligo_infercnv_res_nonmalignant_lq <- oligo_infercnv_res %>% 
  filter(!CellType_clusters%in%c("Malignant") & predicted.class=="Non-Neuronal" & proportion_loss_chr1 > 0.15 & proportion_loss_chr19 > 0.15)

# Mostly removes Unresolved
table(oligo_infercnv_res_nonmalignant_lq$CellType_final_og)

### *****
# Define a malignant cell based on the following definition:
# 1. Gene expression based clustering of cell types is defined as malignant.
# 2. infercnv proportional loss of either chromosome 1 OR** chromosome 19 > 0.15 to make up for sparse signal.
# 3. Remove cells that may have clustered as non-malignant, but may be malignant: Cells not classified as neuronal by Azimuth with greater than 15% of Chr1 and Chr19
# genes with a copy number deletion
### *****

# Assign any non-malignant cell falling into the low-quality non-malignant cell bin as "cells with clonal cna" or "cells without clonal cna"
oligo_sobj_all_harmony_cleaned$CellType_infercnv <- ifelse(oligo_sobj_all_harmony_cleaned$CellType_clusters=="Malignant", "Cells with clonal CNAs", "Cells without clonal CNAs")
# Indicate the low-quality/suspicious **malignant cells.
oligo_sobj_all_harmony_cleaned$CellType_infercnv[oligo_sobj_all_harmony_cleaned$CellType_clusters=="Malignant" & oligo_sobj_all_harmony_cleaned$CellID%in%oligo_infercnv_res_malignant_lq$CellID] <- "Cells without clonal CNAs"
# Indicate the low-quality/suspicious **non-malignant cells.
oligo_sobj_all_harmony_cleaned$CellType_infercnv[oligo_sobj_all_harmony_cleaned$CellType_clusters!="Malignant" & oligo_sobj_all_harmony_cleaned$CellID%in%oligo_infercnv_res_nonmalignant_lq$CellID] <- "Cells with clonal CNAs"

png(paste0(fig_dir, "oligodendroglioma_snrna_post_dfinder_cleaned_umap_celltype_infercnv_20260320.png"), width=10, height=8, res=300, units='in')
DimPlot(oligo_sobj_all_harmony_cleaned, group.by = "CellType_infercnv", pt.size = 1) +
  labs(x="UMAP1", y="UMAP2",color="Inferred copy number status", title="") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "top", legend.justification = "center")
dev.off()

###
# Create a final metadata classification for these cells
oligo_sobj_all_harmony_cleaned$CellType_final <- oligo_sobj_all_harmony_cleaned$CellType_clusters
oligo_sobj_all_harmony_cleaned$CellType_final[oligo_sobj_all_harmony_cleaned$CellType_clusters=="Malignant" & oligo_sobj_all_harmony_cleaned$CellType_infercnv=="Cells without clonal CNAs"] <- "Unresolved"
oligo_sobj_all_harmony_cleaned$CellType_final[oligo_sobj_all_harmony_cleaned$CellType_clusters!="Malignant" & oligo_sobj_all_harmony_cleaned$CellType_infercnv=="Cells with clonal CNAs"] <- "Unresolved"

# For any neuron labelled by a cluster identification, convert it to Unresolved if it very confidently maps to a Non-Neuronal population.
oligo_sobj_all_harmony_cleaned$CellType_final[oligo_sobj_all_harmony_cleaned$CellType_clusters%in%c("ExcNeuron", "InhNeuron") & oligo_sobj_all_harmony_cleaned$predicted.class=="Non-Neuronal" & oligo_sobj_all_harmony_cleaned$predicted.class.score > 0.99] <- "Unresolved"


png(paste0(fig_dir, "oligodendroglioma_snrna_post_dfinder_cleaned_umap_celltype_post_infercnv_20260320.png"), width=10, height=8, res=300, units='in')
DimPlot(oligo_sobj_all_harmony_cleaned, group.by = "CellType_final", pt.size = 1,
  cols =  c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
            "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
            "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7", "Unresolved" = "gray70")) +
  labs(x="UMAP1", y="UMAP2",color="Cell type", title="") +
  theme(legend.position = "top", legend.justification = "center")
dev.off()

# See where the Unresolved cells are located in Harmony UMAP space.
unresolved_cells <- rownames(oligo_sobj_all_harmony_cleaned@meta.data)[oligo_sobj_all_harmony_cleaned@meta.data$CellType_final=="Unresolved"]

png(paste0(fig_dir, "oligodendroglioma_snrna_post_dfinder_cleaned_umap_celltype_post_infercnv_unresolved_cells_20260320.png"), width=10, height=8, res=300, units='in')
DimPlot(oligo_sobj_all_harmony_cleaned, label=F,  cells.highlight = list(unresolved_cells), cols.highlight = c("darkblue"), cols= "grey", pt.size = 1,  raster = FALSE)
dev.off()

## Write out metadata table.
oligo_md_out <- oligo_sobj_all_harmony_cleaned@meta.data

write.table(oligo_md_out, 
            file = paste0(out_data_dir, "/oligodendroglioma_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_infercnv_20260320.txt"),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# Quick sanity check. Both were TRUE.
#oligo_md_out_test <- read.table("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/oligodendroglioma_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_infercnv_20260320.txt", sep = "\t", header = TRUE, row.names = 1)
all(oligo_md_out$CellType_final==oligo_md_out_test$CellType_final)
all(oligo_md_out$CellID==oligo_md_out_test$CellID)

## Load back in the metadata table from the first run.
codel_md_out <- read.delim(file = paste0(out_data_dir, "/codel_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_azimuth_infercnv_20240131.txt"), header = TRUE, sep = "\t", row.names = 1)
sum(codel_md_out$CellID%in%oligo_md_out$CellID)

# Only 50 malignant cells were dropped from the prior analysis.
codel_md_out_dropped <- codel_md_out %>% 
  filter(!CellID%in%oligo_md_out$CellID)
table(codel_md_out_dropped$CellType_final)

# Added 24 malignant cells.
oligo_md_out_added <- oligo_md_out %>% 
  filter(!CellID%in%codel_md_out$CellID)
table(oligo_md_out_added$CellType_final)


# Compare with prior approach for cell type identification.
comparison <- full_join(
  codel_md_out %>% select(CellID, CellType_final),
  oligo_md_out %>% select(CellID, CellType_final),
  by = "CellID",
  suffix = c("_original", "_rerun")
)

# Classify each row
comparison <- comparison %>%
  mutate(
    match_status = case_when(
      is.na(CellType_final_original) ~ "only_in_rerun",
      is.na(CellType_final_rerun) ~ "only_in_original",
      CellType_final_original == CellType_final_rerun ~ "agree",
      TRUE ~ "disagree"
    )
  )

# Summary counts
summary_counts <- comparison %>%
  count(match_status) %>%
  mutate(pct = n / sum(n) * 100)

print(summary_counts)

# Overall agreement % (agree / all cells including unmatched)
overall_agreement <- mean(comparison$match_status == "agree") * 100
cat(sprintf("Overall agreement (including unmatched): %.1f%%\n", overall_agreement))

# Agreement among matched cells only
matched_only <- comparison %>%
  filter(!is.na(CellType_final_original) & !is.na(CellType_final_rerun))

matched_agreement <- mean(matched_only$CellType_final_original == matched_only$CellType_final_rerun) * 100
cat(sprintf("Agreement among matched cells only: %.1f%%\n", matched_agreement))



### END ###