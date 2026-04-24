#################################
# Plot the longitudinal change in mutation burden for the cohort cases with matched normal - initial vs recurrence.
# Author: Kevin Johnson
#################################


library(tidyverse)
library(ggplot2)
library(ggpubr)
library(EnvStats)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
setwd(proj_dir)
source("scripts/utils/plot_theme.R")

# Input data including whether the tumor had evidence of treatment-associated hypermutation. 
# When both WGS and WXS estimates were available, WXS mutation burden was selected as the representative measurement due to higher coverage across gene regions.
mut_burden_pairs <- read.table("data/genetics/matched_normal_blood_longitudinal_mutation_burden.txt", sep = "\t", header = TRUE)

# This is restricted to two timepoints per patient.
pdf(file = paste0(proj_dir, "/figures/fig1c_longitudinal_mut_burden_paper.pdf"), width = 3.25, height = 2.75, bg = "transparent", useDingbats = FALSE)
ggplot(mut_burden_pairs, aes(x = timepoint, y = mut_burden)) + 
  geom_boxplot() +
  geom_line(aes(group = case_barcode, 
                color = ifelse(hypermutant_case == 1, "hypermutant", "non-hm")), 
            linetype = 2,
            size = 0.5) +
  geom_point(aes(color = subtype), size = 0.5) +
  scale_y_log10() +
  scale_color_manual(
    values = c(
      "Astro."      = "#800074",
      "Oligo."      = "#298C8C",
      "hypermutant" = "red",
      "non-hm"      = "gray70")) +
  plot_theme +
  theme(strip.background = element_blank(),
        legend.position  = "right") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label = "p.format") +
  labs(x = "", y = "Mutation per megabase", color = "Type") +
  stat_n_text(size = 2.25)
dev.off()

pdf(file = paste0(proj_dir, "/figures/fig1c_longitudinal_mut_burden_paper_non_hm.pdf"),  width = 3.25, height = 3, bg = "transparent", useDingbats = FALSE)
ggplot(mut_burden_pairs %>% 
         filter(hypermutant_case!=1), aes(x = timepoint, y = mut_burden)) + 
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
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label = "p.format") +
  labs(x = "", y = "Mutation per megabase") +
  stat_n_text(size = 2.25)
dev.off()


### END ###