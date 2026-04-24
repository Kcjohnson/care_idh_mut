##################################
# Report quality metrics for snRNA that ended up in the final analyses
# Author: Kevin Johnson
##################################

library(tidyverse)
library(EnvStats)
library(cowplot)
library(ggpubr)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/"
fig_dir     <- file.path(proj_dir, "figures/")
setwd(proj_dir)

source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")


md <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20260320.txt", sep = "\t", row.names = 1, header = TRUE)
md <- md %>% 
  mutate(tumor_type = case_when(idh_codel_subtype=="IDHmut-codel" ~ "Oligo.",
                              idh_codel_subtype=="IDHmut-noncodel" ~ "Astro.")) %>% 
  mutate(patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2)))

md_trim <- md %>% 
  dplyr::select(SampleID, case_barcode, tumor_type, care_id, patient_id, timepoint) %>% 
  distinct()
rownames(md_trim) <- NULL

# Inspect the median abundance for each cell type
top_cell_prop <- cell_prop %>% 
  group_by(CellType_final) %>% 
  summarise(median_freq = median(freq)) %>% 
  arrange(desc(median_freq))

cellnum_by_time <- md %>% 
  group_by(case_barcode, timepoint) %>% 
  summarise(cell_number = n()) %>% 
  ungroup() %>% 
  mutate(timepoint = recode(timepoint, `T1` = "I",
                            `T2` = "R"))

feature_by_time <- md %>% 
  group_by(case_barcode, timepoint) %>% 
  summarise(median_nFeature = median(nFeature_RNA)) %>% 
  ungroup() %>% 
  mutate(timepoint = recode(timepoint, `T1` = "I",
                            `T2` = "R"))

mt_by_time <- md %>% 
  group_by(case_barcode, timepoint) %>% 
  summarise(median_percent_mt = median(percent.mt)) %>% 
  ungroup() %>% 
  mutate(timepoint = recode(timepoint, `T1` = "I",
                            `T2` = "R"))

all_md <- cellnum_by_time %>% 
  inner_join(mt_by_time, by=c("case_barcode", "timepoint")) %>% 
  inner_join(feature_by_time, by=c("case_barcode", "timepoint")) %>% 
  pivot_longer(cols=c(cell_number, median_percent_mt, median_nFeature), names_to = "metric", values_to = "values")


cellnum_plot <- ggplot(cellnum_by_time %>% filter(timepoint!="T3"), aes(x = timepoint, y = cell_number)) + 
  geom_line(aes(group=case_barcode), color="gray70", linetype=1, alpha = 0.5) +
  geom_point(size = 0.5) +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  labs(x="Time point", y="Median cell number") +
  theme(strip.background = element_blank()) +
  ylim(0, 12500)

pdf(paste0(fig_dir, "longitudinal_cell_number.pdf"), width = 1.5, height = 2.5, useDingbats = FALSE, bg = "transparent")
cellnum_plot
dev.off()

nfeatures <- ggplot(feature_by_time %>% filter(timepoint!="T3"), aes(x = timepoint, y = median_nFeature)) + 
  geom_line(aes(group=case_barcode), color="gray70", linetype=1, alpha = 0.5) +
  geom_point(size = 0.5) +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  labs(x="Time point", y="Median gene count\n(nFeature RNA)") +
  theme(strip.background = element_blank()) +
  ylim(0, 5500)

pdf(paste0(fig_dir, "longitudinal_nfeatures.pdf"), width = 1.5, height = 2.5, useDingbats = FALSE, bg = "transparent")
nfeatures
dev.off()

mt_rna <- ggplot(mt_by_time %>% filter(timepoint!="T3"), aes(x = timepoint, y = median_percent_mt)) + 
  geom_line(aes(group=case_barcode), color="gray70", linetype=1, alpha = 0.5) +
  geom_point(size = 0.5) +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  labs(x="Time point", y="Median percent counts\nfrom mitochondrial RNA") +
  theme(strip.background = element_blank()) +
  ylim(0, 5)

pdf(paste0(fig_dir, "longitudinal_mito.pdf"), width = 1.5, height = 2.5, useDingbats = FALSE, bg = "transparent")
mt_rna 
dev.off()

pdf(paste0(fig_dir, "longitudinal_quality_metrics.pdf"), width = 4, height = 2.75, useDingbats = FALSE, bg = "transparent")
plot_grid(cellnum_plot, nfeatures, mt_rna, ncol = 3)
dev.off()
