##################################
# Purpose: Analyze all publicly available IDH-mutant glioma sc/snRNAseq data for malignant state abundance across tumor grade
# Author: Kevin Johnson
# Date: 2026-03-11
# Input: data/public_all_scrna_snrna_malignant_cell_state_abundance.tsv
# Output: figures/edf7f_public_tumor_type_grade.pdf
##################################

library(tidyverse) # v1.3.1
library(ggpubr) # v0.4.0
library(EnvStats) # v2.7.0

source("utils/plot_theme.R")
all_malignant_abundance <- read.table("data/public_all_scrna_snrna_malignant_cell_state_abundance.tsv", sep = "\t", header = TRUE)

# Remove any sample without grade-defined.
all_malignant_abundance_grade <- all_malignant_abundance %>% 
  filter(!is.na(grade))

# Determine how many unique samples:
sample_num <- n_distinct(all_malignant_abundance_grade$SampleID)

all_malignant_abundance_grade$State <- factor(all_malignant_abundance_grade$State, levels=c("AC-like", "MES-like", "OPC-like", "NPC-like", "Undifferentiated", "Cycling"))
all_malignant_abundance_grade$idh_codel_subtype <- factor(all_malignant_abundance_grade$idh_codel_subtype, levels = c("Oligo.", "Astro."))

# Visualize the data using a Kruskal-Wallis and Wilcoxon since Oligo. samples only have two grade possibilities.
# Apply manual corrections to any minor differences in Kruskal vs. Wilcoxon in Illustrator.
pdf("figures/edf7f_public_tumor_type_grade.pdf", width = 9, height = 6, useDingbats = FALSE, bg = "transparent")
ggplot(all_malignant_abundance_grade, aes(x = as.factor(grade), y = freq*100)) + 
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
  labs(x = paste0("Publicly available sc/snRNA (n = ", sample_num, ")\nKruskal-Wallis"), 
       y = "Relative malignant cell proportion (%)")
dev.off()

ggplot(all_malignant_abundance_grade, aes(x = as.factor(grade), y = freq*100)) + 
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
  labs(x = paste0("Publicly available sc/snRNA (n = ", sample_num, ")\nWilcoxon"), 
       y = "Relative malignant cell proportion (%)")


# Compute statistical results to confirm.
stat_results <- all_malignant_abundance_grade %>%
  group_by(idh_codel_subtype, State) %>%
  summarise(
    p_value = tryCatch({
      if (first(idh_codel_subtype) == "Oligo.") {
        wilcox.test(freq ~ as.factor(grade))$p.value
      } else {
        kruskal.test(freq ~ as.factor(grade))$p.value
      }
    }, error = function(e) NA_real_),
    test_used = if_else(first(idh_codel_subtype) == "Oligo.", "Wilcoxon", "Kruskal-Wallis"),
    .groups = "drop"
  ) %>%
  mutate(signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ ""
  ))


### END ###