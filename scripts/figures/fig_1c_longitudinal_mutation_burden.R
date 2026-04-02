#################################
# Plot the longitudinal change in mutation burden for the cohort cases with matched normal - initial vs recurrence.
# Author: Kevin Johnson
# Updated: 2026.04.01
#################################


library(tidyverse)
library(ggplot2)
library(ggpubr)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
setwd(proj_dir)
source("scripts/utils/plot_theme.R")

# Input data including whether the tumor had evidence of treatment-associated hypermutation. 
mut_burden_pairs <- read.table("data/genetics/matched_normal_blood_longitudinal_mutation_burden.txt", sep = "\t", header = TRUE)

# This is restricted to two timepoints per patient.
pdf(file = paste0(proj_dir, "/figures/fig1c_longitudinal_mut_burden_paper.pdf"), height = 3.5, width = 3,5, useDingbats = FALSE)
ggplot(mut_burden_pairs, aes(x = timepoint, y = mut_burden)) + 
  geom_boxplot() +
  geom_line(aes(group = case_barcode, 
                color = ifelse(hypermutant_case == 1, "hypermutant", "non-hm")), 
            linetype = 2) +
  geom_point(aes(color = subtype)) +
  scale_y_log10() +
  scale_color_manual(
    values = c(
      "Astro."      = "#800074",
      "Oligo."      = "#298C8C",
      "hypermutant" = "red",
      "non-hm"      = "gray70")) +
  plot_theme +
  theme(strip.background = element_blank(),
        legend.position  = "none") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 4, label = "p.format") +
  labs(x = "", y = "Mutation per megabase")
dev.off()


### END ###