##################################
# Genetic and genomic longitudinal features across EDF1d, EDF1e, and EDF1f
# Author: Kevin Johnson
# Date Updated: 2026.03.30
##################################

library(tidyverse)
library(ggpubr)
library(EnvStats)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
source(paste0(proj_dir, "/scripts/utils/plot_theme.R"))
setwd(proj_dir)


sample_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240608.txt",  sep="\t", header = TRUE)
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt",  sep="\t", header = TRUE)

sample_md_seqz <- sample_md %>% 
  filter(timepoint!="T3") %>% 
  dplyr::select(case_barcode:timepoint, seqz_purity_wgs, seqz_purity_wxs) %>% 
  mutate(seqz_purity = ifelse(!is.na(seqz_purity_wgs), seqz_purity_wgs, seqz_purity_wxs)) %>% 
  dplyr::select(case_barcode, timepoint, seqz_purity) %>% 
  pivot_wider(names_from = timepoint, values_from = seqz_purity)

int_vs_recur_seqz <- sample_md_seqz %>% 
  pivot_longer(cols = c(T1, T2), names_to ="timepoint", values_to = "seqz_purity") %>% 
  inner_join(patient_md, by="case_barcode") %>% 
  mutate(timepoint = recode(timepoint, `T1` = "Init.",
                            `T2` = "Recur."),
         type = recode(idh_codel_subtype, `IDHmut-codel` = "Olig.",
                       `IDHmut-noncodel` = "Astro."))

pdf(paste0(fig_dir, "edf1e_dna_purity"), width = 4, height = 4, useDingbats = FALSE)
ggplot(int_vs_recur_seqz, aes(x = timepoint, y = seqz_purity*100)) + 
  geom_boxplot() +
  geom_line(aes(group=case_barcode), color="gray70", linetype=2) +
  geom_point(aes(color=type)) +
  scale_color_manual(values = c("Astro." = "#800074", "Olig." = "#298C8C")) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 4, label="p.format") +
  labs(x="", y="Malignant cell abundance\nSequenza DNA tumor purity (%)", color="Tumor\ntype") +
  theme(strip.background = element_blank()) +
  stat_n_text() +
  ylim(-10, 100)
dev.off()
