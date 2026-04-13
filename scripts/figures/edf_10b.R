##################################
# Purpose: Analyze publicly available longitudinal snRNAseq data for IDH-mutant glioma
# Author: Kevin Johnson
# Date: 2026.03.11
##################################

library(tidyverse) # v1.3.1
library(ggpubr) # v0.4.0
library(EnvStats) # v2.7.0

source("utils/plot_theme.R")
snrna_longitudinal_df <- read.table("data/public_snrna_longitudinal_malignant_cell_state_abundance.tsv", sep = "\t", header = TRUE)

# Reformat to feature only the Mesenchymal state and relabel T1/T2.
snrna_longitudinal_df_plot <- snrna_longitudinal_df %>% 
  filter(State%in%c("MES-like")) %>% 
  mutate(new_timepoint = recode(timepoint_t1t2, `T1` = "Initial",
                                `T2` = "Recur."))

pdf("figures/edf10b_public_longitudinal_mes.pdf", width = 2.5, height = 2.25, useDingbats = FALSE, bg = "transparent")
ggplot(snrna_longitudinal_df_plot, aes(x = new_timepoint, y = freq*100)) + 
  geom_line(aes(group=PatientID), color="gray70", linetype=2) +
  geom_boxplot(outlier.shape = NA, aes(fill=State)) +
  geom_point(size = 0.6) +
  scale_linetype_manual(values="dashed") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "OPC-like"="#E8F5A3",
                             "NPC-like" = "#7fbf7b",
                             "Undifferentiated" = "gray90")) +
  guides(fill = "none") +
  plot_theme +
  stat_compare_means(method = "wilcox",
                     paired = TRUE, # PatientID for T1/T2
                     label="p.format",
                     size = 2.25) +
  stat_n_text(size = 2.25) +
  labs(x="External cohorts snRNA", y="MES-like malignant cell abundance (%)")
dev.off()

### END ###