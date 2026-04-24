##################################
# Create UMAP plots for the two tumor types
# Author: Kevin Johnson
##################################

library(tidyverse)
library(RColorBrewer)
library(viridis)
library(Seurat)
library(harmony)
library(Matrix)
library(ggpubr)
library(openxlsx)
library(DoubletFinder)

# Specify output directory and use the umap theme
fig_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/figures/rna/"
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/umap_theme.R")

# Final metadata with essential information about CellType_final and copy number based cell type assignment
caremut_md <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20260320.txt", sep = "\t", row.names = 1, header = TRUE)

# Need to load in the seurat object and then add the extra metadata
oligo_sobj_all_harmony_cleaned <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/oligodendroglioma_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_20260318.RDS")

# Extract metadata so that we can annotate
seurat_md <- oligo_sobj_all_harmony_cleaned@meta.data

# Identify the variables to add to the seurat object that are missing to avoid duplicate column names
variables_to_keep <- colnames(caremut_md)[!colnames(caremut_md)%in%colnames(seurat_md)]
data_to_add_filt <- caremut_md %>% 
  dplyr::select(CellID, all_of(variables_to_keep))

add_data <- seurat_md %>% 
  left_join(data_to_add_filt, by="CellID")

# Sanity checks to make sure the CellIDs align before adding the metadata
all(row.names(oligo_sobj_all_harmony_cleaned[[]])==add_data$CellID)
row.names(oligo_sobj_all_harmony_cleaned[[]])==add_data$CellID

row.names(add_data) <- row.names(oligo_sobj_all_harmony_cleaned[[]])

oligo_sobj_all_harmony_cleaned <- AddMetaData(oligo_sobj_all_harmony_cleaned, metadata = add_data)

# Determine plotting order for the legend,
cell_state_order <- c("Malignant",
                      "Lymphocyte" ,
                      "Myeloid",
                      "Endothelial",
                      "Mural",
                      "Oligodendrocyte",
                      "Astrocyte",
                      "ExcNeuron",
                      "InhNeuron",
                      "Unresolved")
oligo_sobj_all_harmony_cleaned$CellType_final <- factor(oligo_sobj_all_harmony_cleaned$CellType_final, levels=rev(cell_state_order))


oligo_malignant <- DimPlot(oligo_sobj_all_harmony_cleaned, group.by = "CellType_infercnv",
                           cols =  c("Cells with clonal CNAs" = "black", "Cells without clonal CNAs" = "gray80"),
                           # Manually add point size for larger number of cells
                           pt.size = 1) +
  labs(x="UMAP1", y="UMAP2",color="Inferred copy number status", title="") +
  theme(legend.position = "none") 

ggsave(paste0(fig_dir, "oligo_snrna_malignant_manuscript_umap.pdf"), oligo_malignant, width = 3, height = 3, dpi = 300)



umap_plot <- DimPlot(oligo_sobj_all_harmony_cleaned, group.by = "CellType_final",
                     cols =  c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                               "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
                               "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7", "Unresolved" = "gray70"),
                     # Manually add point size for larger number of cells
                     pt.size = 1) +
  labs(x="UMAP1", y="UMAP2",color="Cell type", title="") +
  theme(legend.position = "none") 

# The pdf version created some shadows that couldn't be resolved by increasing the image size. Not as noticeable in the gray vs. black malignant status
ggsave(paste0(fig_dir, "oligo_snrna_manuscript_umap.pdf"), umap_plot, width = 3, height = 3, dpi = 300)

ggsave(paste0(fig_dir, "oligo_snrna_manuscript_umap.png"), umap_plot, width = 3, height = 3, dpi = 300)


codel_malignant_legend <- DimPlot(oligo_sobj_all_harmony_cleaned, group.by = "CellType_infercnv",
                                  cols =  c("Cells with clonal CNAs" = "black", "Cells without clonal CNAs" = "gray80"),
                                  # Manually add point size for larger number of cells
                                  pt.size = 1) +
  labs(x="UMAP1", y="UMAP2",color="Inferred copy number status", title="") +
  theme(legend.position = "top", legend.justification = "center")


umap_plot <- DimPlot(oligo_sobj_all_harmony_cleaned, group.by = "CellType_final",
                     cols =  c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                               "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
                               "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7", "Unresolved" = "gray70"),
                     # Manually add point size for larger number of cells
                     pt.size = 1) +
  labs(x="UMAP1", y="UMAP2",color="Cell type", title="") +
  theme(legend.position = "top", legend.justification = "center")


g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend_grob <- g_legend(umap_plot)

ggsave(paste0(fig_dir, "oligo_snrna_manuscript_umap_legend.pdf"), legend_grob, width = 8, height = 4, dpi = 300)

legend_grob_malignant <- g_legend(codel_malignant_legend)
ggsave(paste0(fig_dir, "oligo_snrna_malignantmanuscript_umap_legend.pdf"), legend_grob_malignant, width = 8, height = 4, dpi = 300)


### ### ### ### ### ### ### 
### Astro. cohort
### ### ### ### ### ### ### 
astro_sobj_all_harmony_cleaned <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/astrocytoma_seuratobj_filt_singlets_cleaned_mdata_celltype_annotated_20260317.RDS")

seurat_md <- astro_sobj_all_harmony_cleaned@meta.data

# Identify the variables to add to the seurat object that are missing to avoid duplicates
variables_to_keep <- colnames(caremut_md)[!colnames(caremut_md)%in%colnames(seurat_md)]
data_to_add_filt <- caremut_md %>% 
  dplyr::select(CellID, all_of(variables_to_keep))

add_data <- seurat_md %>% 
  left_join(data_to_add_filt, by="CellID")

all(row.names(astro_sobj_all_harmony_cleaned[[]])==add_data$CellID)
row.names(astro_sobj_all_harmony_cleaned[[]])==add_data$CellID

row.names(add_data) <- row.names(astro_sobj_all_harmony_cleaned[[]])

astro_sobj_all_harmony_cleaned <- AddMetaData(astro_sobj_all_harmony_cleaned, metadata = add_data)


astro_sobj_all_harmony_cleaned$CellType_final <- factor(astro_sobj_all_harmony_cleaned$CellType_final, levels=rev(cell_state_order))


astro_malignant <- DimPlot(astro_sobj_all_harmony_cleaned, group.by = "CellType_infercnv",
                           cols =  c("Cells with clonal CNAs" = "black", "Cells without clonal CNAs" = "gray80"),
                           # Manually add point size for larger number of cells
                           pt.size = 1) +
  labs(x="UMAP1", y="UMAP2",color="Inferred copy number status", title="") +
  theme(legend.position = "none") 

ggsave(paste0(fig_dir, "astro_snrna_malignant_manuscript_umap.pdf"), astro_malignant, width = 3, height = 3, dpi = 300)


umap_plot <- DimPlot(astro_sobj_all_harmony_cleaned, group.by = "CellType_final",
                     cols =  c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                               "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
                               "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7", "Unresolved" = "gray70"),
                     # Manually add point size for larger number of cells
                     pt.size = 1) +
  labs(x="UMAP1", y="UMAP2",color="Cell type", title="") +
  theme(legend.position = "none") 

# The pdf version created some shadows that couldn't be resolved by increasing the image size. Not as noticeable in the gray vs. black malignant status
ggsave(paste0(fig_dir, "astro_snrna_manuscript_umap.pdf"), umap_plot, width = 3, height = 3, dpi = 300)

ggsave(paste0(fig_dir, "astro_snrna_manuscript_umap.png"), umap_plot, width = 3, height = 3, dpi = 300)


codel_malignant_legend <- DimPlot(astro_sobj_all_harmony_cleaned, group.by = "CellType_infercnv",
                                  cols =  c("Cells with clonal CNAs" = "black", "Cells without clonal CNAs" = "gray80"),
                                  # Manually add point size for larger number of cells
                                  pt.size = 1) +
  labs(x="UMAP1", y="UMAP2",color="Inferred copy number status", title="") +
  theme(legend.position = "top", legend.justification = "center")


umap_plot <- DimPlot(astro_sobj_all_harmony_cleaned, group.by = "CellType_final",
                     cols =  c("Malignant" = "#FB8072", "astrodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                               "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
                               "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7", "Unresolved" = "gray70"),
                     # Manually add point size for larger number of cells
                     pt.size = 1) +
  labs(x="UMAP1", y="UMAP2",color="Cell type", title="") +
  theme(legend.position = "top", legend.justification = "center")


g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend_grob <- g_legend(umap_plot)

ggsave(paste0(fig_dir, "astro_snrna_manuscript_umap_legend.pdf"), legend_grob, width = 8, height = 4, dpi = 300)

legend_grob_malignant <- g_legend(codel_malignant_legend)
ggsave(paste0(fig_dir, "astro_snrna_malignantmanuscript_umap_legend.pdf"), legend_grob_malignant, width = 8, height = 4, dpi = 300)


### END ###