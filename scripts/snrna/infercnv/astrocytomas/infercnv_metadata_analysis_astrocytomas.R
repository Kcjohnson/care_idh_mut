##################################
# Analyze infercnv metadata to revise uncertain malignant classification for CARE IDH-mutant astrocytomas
# Author: Kevin Johnson
# Date Updated: 2026.03.20
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
infercnv_data_dir <- "/vast/palmer/pi/verhaak/kcj28/care_mut/results/infercnv/astrocytoma_samples"
fig_dir <- "/vast/palmer/pi/verhaak/kcj28/care_mut/results/figures/rna/"
source("/vast/palmer/pi/verhaak/kcj28/care_mut/scripts/misc/plot_theme.R")

################################################
## Step 1: Load the results from Seurat pre-processing:
################################################
# This represents the metadata as processed through Seurat for removing singlets, harmony batch correction based on case barcode + lab, and annotating clusters based on previously determined marker genes.
# I removed residual doublet clusters following an iterative approach of clustering, infercnv, and re-clustering with different data inputs.
caremut_md <- read.table("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/astrocytoma_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_20260317.txt", sep = "\t", header = TRUE)

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

# Define cohort-wide malignant cells. 
malignant_rows <- caremut_md %>% filter(CellType_clusters=="Malignant")
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
      
      # Calculate each cell's CNV correlation with the average tumor cell. NOTE: That this is looking at the proportion of a chromosome with amp/deletion based on genes - not specific segments.
      malignant_data <- data[rownames(data)%in%malignant_cells, ]
      malignant_data_filt <- malignant_data %>% 
        dplyr::select(starts_with("proportion_loss_"), starts_with("proportion_dupli_"))
      
      # Assumes the matrix is named 'malignant_data_filt'
      column_averages <- colMeans(malignant_data_filt)
      
      # Restricting to the features that are proportional loss or duplicated.
      data_filt <- data  %>% 
        dplyr::select(starts_with("proportion_loss_"), starts_with("proportion_dupli_"))
      
      # Compute correlations between each row of 'malignant_data_filt' and 'column_averages'
      row_correlations <- apply(data_filt, 1, function(row) cor(row, column_averages))
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

# warnings : In cor(row, column_averages) : the standard deviation is zero

# Combine all data frames into a single data frame
combined_data <- do.call(rbind, data_list)
combined_data$CellID <- sapply(strsplit(rownames(combined_data), "\\."), last)

# There are duplicated cell IDs because I had to use NL26-0 normal cells as the reference cells for NL26-1 since it had no non-malignant cells
duplicated_cells <- which(duplicated(combined_data$CellID))
combined_data <- combined_data[-duplicated_cells, ]
rownames(combined_data) <- sapply(strsplit(rownames(combined_data), "\\."), last)

# inferCNV can fail with too few cells per cluster/celltype so there were 125 non-malignant cells where infercnv was not run.
cells_without_infercnv <- caremut_md %>%
  group_by(SampleID, CellType_clusters) %>% 
  summarise(counts = n()) %>% 
  filter(counts < 5)

table(cells_without_infercnv$CellType_clusters)
sum(cells_without_infercnv$counts) # This is where the 125 cells without infercnv comes from
dim(caremut_md)[1]-dim(combined_data)[1]

# Retrieve the continuous scores. For example, "proportion_dupli_chr7" is how many of the genes kept after filtering on chr7 are part of a duplication. 
infercnv_continuous <- combined_data %>% 
  dplyr::select(CellID, sample_id = dir_name, malignant_cor, starts_with("proportion_loss_"), starts_with("proportion_dupli_")) %>% 
  right_join(caremut_md, by="CellID")
# Correlation between two runs of infercnv with slightly different inputs. Correlation coefficient = 0.949
cor.test(infercnv_continuous$malignant_cor, infercnv_continuous$malignant_cor_old, na.action = "na.omit")

# Define a somatic copy number alteration (scna) burden metric across these cells by taking the row average:
infercnv_merge <- combined_data %>% 
  dplyr::select(CellID, sample_id = dir_name, malignant_cor, starts_with("proportion_loss_"), starts_with("proportion_dupli_"))
# Confirm that the appropriate columns are used.
head(infercnv_merge[,4:47])
infercnv_merge$scna_burden <-  rowMeans(infercnv_merge[,4:47])

# Add Azimuth normal brain cell type annotation:
astrocytoma_azimuth_md <- read.table("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/astrocytoma_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_azimuth_20260320.txt", header = TRUE, sep = "\t", row.names = 1)

# Tidy up the data to be merged.
infercnv_merge_trim <- caremut_md %>% 
  # Drop the prior classes/subclasses
  dplyr::select(-predicted.class, -predicted.subclass) %>% 
  inner_join(astrocytoma_azimuth_md, by="CellID") %>% 
  left_join(infercnv_merge, by="CellID") %>% 
  # Create a shorthand ID for cells that are likely malignant or non-malignant based on expression clustering.
  mutate(CellType_binary = ifelse(CellType_clusters=="Malignant", "Malignant", "Non-malignant"))

# Load back in the Seurat object with these IDH astrocytoma data:
astro_sobj_all_harmony_cleaned <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/astrocytoma_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_20260317.RDS")
astro_sobj_all_harmony_cleaned@meta.data$scna_burden <- NULL
astro_sobj_all_harmony_cleaned@meta.data$malignant_cor <- NULL
astro_sobj_all_harmony_cleaned@meta.data$predicted.class <- NULL
astro_sobj_all_harmony_cleaned@meta.data$predicted.subclass <- NULL

# Avoid duplicating variable names already present in the metadata
variables_to_keep <- colnames(infercnv_merge_trim)[!colnames(infercnv_merge_trim)%in%colnames(astro_sobj_all_harmony_cleaned@meta.data)]
infercnv_merge_trim_to_add <- infercnv_merge_trim %>% 
  dplyr::select(CellID, all_of(variables_to_keep))

# Add the metadata of the infercnv results
seurat_md <- astro_sobj_all_harmony_cleaned@meta.data
add_data <- seurat_md %>% 
  left_join(infercnv_merge_trim_to_add, by="CellID")

# Double-check that all cellular identifiers match.
all(add_data$CellID==row.names(astro_sobj_all_harmony_cleaned[[]]))
row.names(add_data) <- row.names(astro_sobj_all_harmony_cleaned[[]])
astro_sobj_all_harmony_cleaned <- AddMetaData(astro_sobj_all_harmony_cleaned, metadata = add_data)


# Inspect the UMAP based on the labels previously analyzed.
# Expression-based.
DimPlot(astro_sobj_all_harmony_cleaned, group.by = "CellType_binary", pt.size = 1)
# The predicted classes and subclasses do align with what's expected.
DimPlot(astro_sobj_all_harmony_cleaned, group.by = "predicted.class", pt.size = 1)
DimPlot(astro_sobj_all_harmony_cleaned, group.by = "predicted.subclass", pt.size = 1)
DimPlot(astro_sobj_all_harmony_cleaned, group.by = "CellType_clusters", label = TRUE, pt.size = 1)
DimPlot(astro_sobj_all_harmony_cleaned, group.by = "RNA_snn_res.0.6", label = TRUE, pt.size = 1)
table(astro_sobj_all_harmony_cleaned$SampleID, astro_sobj_all_harmony_cleaned$RNA_snn_res.0.6)

# Define a CNA signal using the proportion of duplication for both Chr7, Chr8 plus loss for Chr10, Chr13 (common events in astrocytomas)
astro_sobj_all_harmony_cleaned$proportion_chr78dupli_1013loss  = astro_sobj_all_harmony_cleaned$proportion_dupli_chr7+astro_sobj_all_harmony_cleaned$proportion_dupli_chr8+astro_sobj_all_harmony_cleaned$proportion_loss_chr10+astro_sobj_all_harmony_cleaned$proportion_loss_chr13

# Visualize the copy number differences across the UMAP
FeaturePlot(astro_sobj_all_harmony_cleaned, features = "proportion_dupli_chr7", pt.size = 1)
FeaturePlot(astro_sobj_all_harmony_cleaned, features = "proportion_dupli_chr8", pt.size = 1)
FeaturePlot(astro_sobj_all_harmony_cleaned, features = "proportion_loss_chr10", pt.size = 1)
FeaturePlot(astro_sobj_all_harmony_cleaned, features = "proportion_loss_chr13", pt.size = 1)
FeaturePlot(astro_sobj_all_harmony_cleaned, features = "proportion_loss_chr22", pt.size = 1)
FeaturePlot(astro_sobj_all_harmony_cleaned, features = "proportion_chr78dupli_1013loss", pt.size = 1)
FeaturePlot(astro_sobj_all_harmony_cleaned, features = "scna_burden", pt.size = 1)
FeaturePlot(astro_sobj_all_harmony_cleaned, features = "malignant_cor", pt.size = 1)

cor.test(astro_sobj_all_harmony_cleaned$malignant_cor_old, astro_sobj_all_harmony_cleaned$malignant_cor, method = "p")
cor.test(astro_sobj_all_harmony_cleaned$scna_burden_old, astro_sobj_all_harmony_cleaned$scna_burden, method = "p")

# Examine the distribution of these values across Seurat clusters
VlnPlot(astro_sobj_all_harmony_cleaned, group.by = "CellType_clusters", features = "proportion_dupli_chr7")
VlnPlot(astro_sobj_all_harmony_cleaned, group.by = "CellType_clusters", features = "scna_burden")
VlnPlot(astro_sobj_all_harmony_cleaned, group.by = "CellType_clusters", features = "malignant_cor")

# Export the data.frame to manipulate.
astro_infercnv_res <- astro_sobj_all_harmony_cleaned@meta.data

# Investigate a reasonable cutoff for malignant definition based on scna burden and malignant correlation.
astro_infercnv_res_malignant <- astro_infercnv_res %>% 
  filter(CellType_clusters%in%c("Malignant"))

pdf(paste0(fig_dir, "astrocytoma_malignant_scna_correlations.pdf"), width = 5, height = 4, useDingbats = FALSE)
ggplot(astro_infercnv_res, aes(x=malignant_cor)) +
  geom_histogram(binwidth=0.05) +
  plot_theme +
  facet_grid(CellType_clusters~., scales = "free") 
dev.off()


# Setting a threshold off of the malignant correlation OR the somatic copy number alteration burden. Don't want to penalize samples that have subclonal CNAs or very few CNAs and may not strongly correlate
# This also recovers samples that may have a low burden estimate, but high correlation with other malignant cells.
# 93% of these cells would also be present in a stricter analysis where malignant correlation was set to 0.5 and/or higher burden requirement
sum(astro_infercnv_res_malignant$malignant_cor>0.5 | astro_infercnv_res_malignant$scna_burden>0.4, na.rm = TRUE)/dim(astro_infercnv_res_malignant)[1]
# 97% at this level:
sum(astro_infercnv_res_malignant$malignant_cor>0.5 | astro_infercnv_res_malignant$scna_burden>0.15, na.rm = TRUE)/dim(astro_infercnv_res_malignant)[1]
# 99% at this level:
sum(astro_infercnv_res_malignant$malignant_cor>0.3 | astro_infercnv_res_malignant$scna_burden>0.15, na.rm = TRUE)/dim(astro_infercnv_res_malignant)[1]

# More stringent option:
stringent_threshold <- astro_infercnv_res_malignant %>% 
  filter(CellType_clusters=="Malignant", malignant_cor > 0.5 | scna_burden>0.15)
# Less stringent option requiring detectable burden or correlation
less_stringent_threshold <- astro_infercnv_res_malignant %>% 
  filter(malignant_cor > 0.3 | scna_burden>0.15)

cells_in_limbo <- anti_join(less_stringent_threshold, stringent_threshold, by="CellID")
# This analysis indicates that 36/45 samples would have cells in this category suggesting that it's not dominated by a few low quality samples.
n_distinct(cells_in_limbo$SampleID)

# Defining a set of low quality *MALIGNANT* cells that we will assign "unresolved" status.
astro_infercnv_res_malignant_lq <- astro_infercnv_res_malignant %>% 
  filter(CellType_clusters%in%c("Malignant"), malignant_cor < 0.3 & scna_burden < 0.15)

# Define the low-quality *NON-malignant* cell states based on high infercnv signal. Including the Azimuth Non-Neuronal classification because neurons tend to have high inferred CNA signal.
astro_infercnv_res_nonmalignant_lq <- astro_infercnv_res %>% 
  filter(!CellType_clusters%in%c("Malignant"), predicted.class=="Non-Neuronal", malignant_cor > 0.3, scna_burden > 0.15)

# Most common cell type conflation appears to be Astrocytes and Mural cells. Dominated by a few samples.
table(astro_infercnv_res_nonmalignant_lq$CellType_final_og)
table(astro_infercnv_res_nonmalignant_lq$CellType_clusters)
table(astro_infercnv_res_nonmalignant_lq$SampleID)

### *****
# Define a malignant cell based on the following definition:
# 1. Gene expression based clustering of cell types is defined as malignant.
# 2. infercnv CNA signal (scna burden) scna_burden > 0.15 OR** malignant_cor > 0.3.
# 3. Remove cells that may have clustered as non-malignant, but may be malignant: Cells not classified as neuronal by Azimuth with greater than 0.15 scna burden AND** a mean malignant cell correlation of greater than 0.3.
### *****

# Assign any non-malignant cell falling into the low-quality non-malignant cell bin as "cells with clonal cna" or "cells without clonal cna".
# Set base group to be revised by alternative fields:
astro_sobj_all_harmony_cleaned$CellType_infercnv <- ifelse(astro_sobj_all_harmony_cleaned$CellType_clusters%in%c("Malignant"), "Cells with clonal CNAs", "Cells without clonal CNAs")
# Indicate the low-quality/suspicious malignant cells.
astro_sobj_all_harmony_cleaned$CellType_infercnv[astro_sobj_all_harmony_cleaned$CellType_clusters%in%c("Malignant") & astro_sobj_all_harmony_cleaned$CellID%in%astro_infercnv_res_malignant_lq$CellID] <- "Cells without clonal CNAs"
# Indicate the low-quality/suspicious **non-malignant cells.
astro_sobj_all_harmony_cleaned$CellType_infercnv[!astro_sobj_all_harmony_cleaned$CellType_clusters%in%c("Malignant") & astro_sobj_all_harmony_cleaned$CellID%in%astro_infercnv_res_nonmalignant_lq$CellID] <- "Cells with clonal CNAs"

png(paste0(fig_dir, "astrocytoma_snrna_post_dfinder_cleaned_umap_celltype_infercnv_20260319.png"), width=10, height=8, res=300, units='in')
DimPlot(astro_sobj_all_harmony_cleaned, group.by = "CellType_infercnv", pt.size = 1) +
  labs(x="UMAP1", y="UMAP2",color="Inferred copy number status", title="") +
  scale_color_brewer(palette = "Set2") +
  theme(legend.position = "top", legend.justification = "center")
dev.off()


#### 
# Create a final metadata classification for these cells - default to the initial CellType assessment for those 125 cells that were not submitted to infercnv, where infercnv estimates could be unstable.
astro_sobj_all_harmony_cleaned$CellType_final <- astro_sobj_all_harmony_cleaned$CellType_clusters
# Two definitions for the "Unresolved" category.
astro_sobj_all_harmony_cleaned$CellType_final[astro_sobj_all_harmony_cleaned$CellType_clusters%in%c("Malignant") & astro_sobj_all_harmony_cleaned$CellType_infercnv=="Cells without clonal CNAs"] <- "Unresolved"
astro_sobj_all_harmony_cleaned$CellType_final[!astro_sobj_all_harmony_cleaned$CellType_clusters%in%c("Malignant") & astro_sobj_all_harmony_cleaned$CellType_infercnv=="Cells with clonal CNAs"] <- "Unresolved"
# For any neuron labelled by cluster identification, convert it to Unresolved if it VERY confidently maps to a Non-Neuronal population. These might reflect overclustering due to batch correction.
astro_sobj_all_harmony_cleaned$CellType_final[astro_sobj_all_harmony_cleaned$CellType_clusters%in%c("ExcNeuron", "InhNeuron") & astro_sobj_all_harmony_cleaned$predicted.class=="Non-Neuronal" & astro_sobj_all_harmony_cleaned$predicted.class.score > 0.99] <- "Unresolved"

# astro_sobj_all_harmony_cleaned$CellType_final[astro_sobj_all_harmony_cleaned$CellType_final=="Malignant" & is.na(astro_sobj_all_harmony_cleaned$scna_burden)] <- "Unresolved"

table(astro_sobj_all_harmony_cleaned$CellType_final, astro_sobj_all_harmony_cleaned$predicted.class)

png(paste0(fig_dir, "astrocytoma_snrna_post_dfinder_cleaned_umap_celltype_post_infercnv_20260320.png"), width=10, height=8, res=300, units='in')
DimPlot(astro_sobj_all_harmony_cleaned, group.by = "CellType_final",
        cols =  c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                  "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
                  "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7", "Unresolved" = "gray70"),
        pt.size = 1) +
  labs(x="UMAP1", y="UMAP2",color="Cell type", title="") +
  theme(legend.position = "top", legend.justification = "center")
dev.off()

## Write out metadata table.
astro_md_out <- astro_sobj_all_harmony_cleaned@meta.data

# Reclassify the binary classification for malignant and non-malignant based on the CellType_final classification
astro_md_out <- astro_md_out %>% 
  mutate(CellType_binary = ifelse(CellType_final=="Malignant", "Malignant", "Non-malignant"))

# Use this output to create a final cell type metadata file.
write.table(astro_md_out, 
            file = paste0("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/astrocytoma_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_infercnv_20260320.txt"),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)


# Load in the original classification.
noncodel_md_out <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/noncodel_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_azimuth_infercnv_20240212.txt", header = TRUE, sep = "\t", row.names = 1)
sum(noncodel_md_out$CellID%in%astro_md_out$CellID)/length(noncodel_md_out$CellID)
sum(astro_md_out$CellID%in%noncodel_md_out$CellID)/length(astro_md_out$CellID)

# 113 malignant cells were dropped in the prior analysis. More concentrated in a few samples.
noncodel_md_out_dropped <- noncodel_md_out %>% 
  filter(!CellID%in%astro_md_out$CellID)
table(noncodel_md_out_dropped$CellType_final)

# Added 380 malignant cells. Largely spread out across samples.
astro_md_out_added <- astro_md_out %>% 
  filter(!CellID%in%noncodel_md_out$CellID)
table(astro_md_out_added$CellType_final)

# Compare with prior approach for cell type identification.
comparison <- full_join(
  noncodel_md_out %>% select(CellID, CellType_final),
  astro_md_out %>% select(CellID, CellType_final),
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


# Create a visualization based on the cutoffs set for malignant cell state classification.
pdf(paste0(fig_dir, "astrocytoma_snrna_malignant_cna_correlation_cna_burden_cutoffs.pdf"), width = 5, height = 4, useDingbats = FALSE)
ggplot(astro_md_out %>% 
         filter(CellType_final=="Malignant"), aes(x = scna_burden, y = malignant_cor, color=CellType_final)) + 
  geom_point(alpha=0.2) +
  ylab("CNA malignant correlation") +
  xlab("CNA signal (infercnv)") +
  labs(color = "Cell type") +
  plot_theme +
  theme(legend.position="bottom") +
  geom_vline(xintercept = 0.15, linetype="dotted", 
             color = "blue", size=1) +
  geom_hline(yintercept = 0.3, linetype="dotted", 
             color = "blue", size=1)
dev.off()

# Instead of individual points, take the average values for each cell type and plot by case + binary CellType.
# Filtering out celltypes with low numbers of cells, which introduce noise to the correlations.
case_plot <- astro_md_out %>% 
  filter(CellType_final!="Unresolved") %>% 
  group_by(case_barcode, CellType_final, CellType_binary) %>% 
  summarise(avg_malignant_cor = median(malignant_cor, na.rm=T),
            avg_burden = median(scna_burden, na.rm=T),
            counts = n()) %>% 
  ungroup() %>% 
  # Restricting to higher confidence metrics rather than a cell type being influenced by an outlier
  filter(counts > 40) %>% 
  group_by(case_barcode, CellType_binary) %>% 
  summarise(meta_malignant_cor = median(avg_malignant_cor, na.rm=T),
            meta_burden = median(avg_burden, na.rm=T)) 

pdf(paste0(fig_dir, "astrocytoma_snrna_correlation_snrna_malignant_infercnv_cna_by_case.pdf"), width = 5, height = 4, useDingbats = FALSE)
ggplot(case_plot, aes(x = meta_burden, y = meta_malignant_cor, color=CellType_binary)) + 
  geom_point(size = 3) +
  ylab("Median CNA malignant correlation\n(infercnv w/ rna malignant cells") +
  xlab("Median CNA signal (infercnv)") +
  labs(color = "Cell type") +
  plot_theme +
  scale_color_manual(values=c("Non-malignant" = "#D9D9D9", 
                              "Malignant" = "#FB8072")) +
  theme(legend.position="bottom") 
dev.off()


sample_plot_matched <- astro_md_out %>% 
  filter(CellType_final!="Unresolved") %>% 
  group_by(SampleID, CellType_final, CellType_binary) %>% 
  summarise(avg_malignant_cor = median(malignant_cor, na.rm=T),
            avg_burden = median(scna_burden, na.rm=T),
            counts = n()) %>% 
  ungroup() %>% 
  filter(counts > 40) %>% 
  group_by(SampleID, CellType_binary) %>% 
  summarise(meta_malignant_cor = median(avg_malignant_cor, na.rm=T),
            meta_burden = median(avg_burden, na.rm=T)) 

pdf(paste0(fig_dir, "astrocytoma_snrna_correlation_snrna_malignant_infercnv_cna_by_sample.pdf"), width = 5, height = 4, useDingbats = FALSE)
ggplot(sample_plot_matched, aes(x = meta_burden, y = meta_malignant_cor, color=CellType_binary)) + 
  geom_point(size = 3) +
  ylab("Median CNA malignant correlation") +
  xlab("Median CNA signal (infercnv)") +
  labs(color = "Cell type") +
  plot_theme +
  scale_color_manual(values=c("Non-malignant" = "#D9D9D9", 
                              "Malignant" = "#FB8072")) +
  theme(legend.position="bottom")
dev.off()

### END ###