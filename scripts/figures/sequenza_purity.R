#################################
# Plot the longitudinal difference in Sequenza estimated tumor purity - initial vs recurrence.
# Author: Kevin Johnson
#################################


library(tidyverse)
library(ggplot2)
library(ggpubr)
library(EnvStats)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
setwd(proj_dir)
source("scripts/utils/plot_theme.R")

# When multiple Sequenza cellularity estimates were available, the preferred option was WGS
seqz_input <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/genetics/sequenza_longitudinal_cellularity_estimates.txt", sep = "\t", header = TRUE)

# This is restricted to two timepoints per patient.
pdf(file = paste0(proj_dir, "/figures/sequenza_dna_tumor_purity.pdf"), height = 3,  width = 2.5, useDingbats = FALSE, bg = "transparent")
ggplot(seqz_input, aes(x = timepoint, y = seqz_cellularity*100)) + 
  geom_boxplot() +
  geom_line(aes(group=case_barcode), color="gray70", linetype=2, size = 0.5) +
  geom_point(aes(color=tumor_type), size = 0.5) +
  scale_color_manual(values = c("Astro." = "#800074", "Olig." = "#298C8C")) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  labs(x="", y="Malignant cell abundance\nSequenza DNA tumor purity (%)", color="Tumor\ntype") +
  theme(strip.background = element_blank()) +
  stat_n_text(size = 2.25) +
  ylim(-10, 100)
dev.off()


### END ###