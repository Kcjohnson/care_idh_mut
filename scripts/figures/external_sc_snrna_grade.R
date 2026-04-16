##################################
# Plot the malignant cell state differences across tumor grade by tumor type
# Author: Kevin Johnson
##################################

library(tidyverse)
library(ggpubr)
library(EnvStats)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/"
fig_dir     <- file.path(proj_dir, "figures/")
setwd(proj_dir)

source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")

# Malignant cell state assignment was performed separately per cohort.
mal_longitudinal_abundance <- read.delim(paste0(proj_dir, "data/public/public_all_scrna_snrna_malignant_cell_state_abundance.tsv"), header = TRUE, sep = "\t")

mal_longitudinal_abundance$State <- factor(mal_longitudinal_abundance$State, levels=c("AC-like", "MES-like", "OPC-like", "NPC-like", "Undifferentiated", "Cycling"))
mal_longitudinal_abundance$idh_codel_subtype <- factor(mal_longitudinal_abundance$idh_codel_subtype, levels=c("Oligo.", "Astro."))


pdf(paste0(fig_dir, "public_all_samples_grade_n139_kruskal.pdf"), width = 4, height = 3, useDingbats = FALSE, bg = "transparent")
ggplot(mal_longitudinal_abundance %>% 
         filter(!is.na(grade)), aes(x = as.factor(grade), y = freq*100)) + 
  geom_boxplot(aes(fill=State), outlier.shape = NA) +
  geom_point(size = 0.6) +
  plot_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "OPC-like"="#E8F5A3",
                             "NPC-like" = "#7fbf7b",
                             "Undifferentiated" = "gray90")) +
  stat_compare_means(method = "kruskal", paired = FALSE, size = 2.25, label="p.format") +
  facet_grid(idh_codel_subtype~State, scales="free") +
  stat_n_text(size = 2.25) +
  labs(x="Publicly available sc/snRNA (n = 139)", y="Relative malignant cell proportion (%)")
dev.off()

# I couldn't determine how to split the Kruskal and Wilcox tests across the facets so using kruskal for Astro. and Wilcox for Oligo.
pdf(paste0(fig_dir, "public_all_samples_grade_n139_wilcox.pdf"), width = 4, height = 3, useDingbats = FALSE, bg = "transparent")
ggplot(mal_longitudinal_abundance %>% 
         filter(!is.na(grade)), aes(x = as.factor(grade), y = freq*100)) + 
  geom_boxplot(aes(fill=State), outlier.shape = NA) +
  geom_point(size = 0.6) +
  plot_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "OPC-like"="#E8F5A3",
                             "NPC-like" = "#7fbf7b",
                             "Undifferentiated" = "gray90")) +
  stat_compare_means(method = "wilcox", paired = FALSE, size = 2.25, label="p.format") +
  facet_grid(idh_codel_subtype~State, scales="free") +
  stat_n_text(size = 2.25) +
  labs(x="Publicly available sc/snRNA (n = 139)", y="Relative malignant cell proportion (%)")
dev.off()

