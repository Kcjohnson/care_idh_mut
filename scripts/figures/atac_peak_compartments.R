##################################
# Inspect the longitudinal differentially accessible peaks across major cell types
# Author: Kevin Johnson
# Date Updated: 2026.04.07
##################################

library(tidyverse)
library(ggpubr)
library(EnvStats)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
setwd(proj_dir)
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")

# Longitudinal differences in ATAC peak chromatin accessibility. CellType controlled longidutinal analysis (e.g., P100T1-Myeloid vs P100T2-Myeloid).
atac_patient_celtype_long <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/archr_care_tme_celltype_longitudinal_peaks_min40cells_per_celltype.txt", sep = "\t", header = TRUE)

atac_patient_celtype_long <- atac_patient_celtype_long %>% 
  mutate(compartment = ifelse(`cell_state.j.`=="Malignant", "Malignant","TME"),
         `Cell Type` = `cell_state.j.`,
         all_peaks = up_peaks+down_peaks) %>% 
  filter(use_cells > 49, bgd_cells > 49)

pdf(file = "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/figures/longitudinal_atac_peak_accessibility.pdf", height = 2.5, width = 3.25, bg = "transparent", useDingbats = FALSE)
ggplot(atac_patient_celtype_long, aes(x=compartment, y=all_peaks)) +
  geom_boxplot(
    outlier.shape = NA) +
  geom_point(aes(color = `Cell Type`), position = position_jitter(width = 0.1, seed = 42), 
             size = 0.8) +
  labs(x = "Compartment", y = "Per patient longitudinal differential\npeak accessibility per patient") +
  scale_color_manual(values=c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                             "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
                             "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7")) +
  stat_compare_means(method = "wilcox", label = "p.format", size = 2.25) +
  plot_theme +
  theme(
    legend.position = "right"
  ) +
    stat_n_text(size = 2.25) 
dev.off()

### END ###