##################################
# Plot the longitudinal shift in MES-like malignant cells in publicly available snRNAseq
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
mal_longitudinal_abundance <- read.delim(paste0(proj_dir, "data/public/public_snrna_longitudinal_malignant_cell_state_abundance.tsv"), header = TRUE, sep = "\t")

mal_longitudinal_abundance$State <- factor(mal_longitudinal_abundance$State, levels=c("AC-like", "MES-like", "OPC-like", "NPC-like", "Undifferentiated", "Cycling"))

# Restrict to Astrocytomas since there are 12 astrocytomas and 1 oligodendroglioma
pdf(paste0(fig_dir, "public_longitudinal_mes_abundance.pdf"), width = 2.5, height = 2.25, useDingbats = FALSE, bg = "transparent")
ggplot(mal_longitudinal_abundance %>% 
         filter(idh_codel_subtype=="Astro.", State=="MES-like") %>% 
         mutate(new_timepoint = recode(timepoint, `T1` = "Initial",
                                       `T2` = "Recur.")), aes(x = new_timepoint, y = freq*100)) + 
  geom_boxplot(outlier.shape = NA, aes(fill=State)) +
  geom_line(aes(group=PatientID), color="gray70", linetype=2) +
  geom_point(size = 0.6) +
  facet_grid(.~State) +
  scale_linetype_manual(values="dashed") +
  scale_fill_manual(values=c("MES-like"="#F77D58")) +
  guides(fill=FALSE) +
  plot_theme +
  stat_compare_means(method = "wilcox",
                     paired = TRUE,
                     label="p.format",
                     size = 2.25) +
  stat_n_text(size = 2.25) +
  labs(x="External cohorts snRNA", y="MES-like malignant cell abundance (%)")
dev.off()


