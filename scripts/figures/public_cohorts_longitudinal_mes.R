##################################
# Purpose: Analyze publicly available longitudinal snRNAseq malignant data for IDH-mutant glioma
# Author: Kevin Johnson
##################################

library(tidyverse) # v1.3.1
library(ggpubr) # v0.4.0
library(EnvStats) # v2.7.0

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/"
fig_dir     <- file.path(proj_dir, "figures/")
setwd(proj_dir)

source(paste0(proj_dir, "scripts/utils/plot_theme.R"))


snrna_longitudinal_df <- read.table(paste0(proj_dir, "data/public/public_snrna_longitudinal_malignant_cell_state_abundance.tsv"), sep = "\t", header = TRUE)

# Reformat to feature only the Mesenchymal state and relabel T1/T2.
snrna_longitudinal_df_plot <- snrna_longitudinal_df %>% 
  filter(State%in%c("MES-like")) %>% 
  mutate(new_timepoint = recode(timepoint, `T1` = "Initial",
                                `T2` = "Recur."))

# There are 13 longitudinal samples. 1 sample is an oligodendroglioma. It doesn't really matter whether we remove it in terms of statistical significance.
# Cleaner for it to be 12 astrocytomas.
pdf(paste0(fig_dir, "edf10_public_longitudinal_mes.pdf"), width = 2.5, height = 2.25, useDingbats = FALSE, bg = "transparent")
ggplot(snrna_longitudinal_df_plot %>% 
         filter(idh_codel_subtype=="Astro."), aes(x = new_timepoint, y = freq*100)) + 
  geom_boxplot(outlier.shape = NA, aes(fill=State)) +
  geom_line(aes(group=PatientID), color="gray70", linetype=2) +
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