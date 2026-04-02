##################################
# Report quality metrics for snRNA that ended up in the final analyses
# Author: Kevin Johnson
# Date Updated: 2024.06.12
##################################

library(tidyverse)
library(RColorBrewer)
library(viridis)
library(Seurat)
library(harmony)
library(presto)
library(EnvStats)
library(ggpubr)
source("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/scripts/misc/plot_theme.R")
source("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/scripts/misc/caremut_utils.R")
source("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/scripts/misc/umap_theme.R")

fig_dir <- "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/figures/rna/"

# Load in the CAREmut metadata - updated 2024.02.12
md <- read.table("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20240212.txt", sep = "\t", row.names = 1, header = TRUE)
md$subtype <- ifelse(md$idh_codel_subtype=="IDHmut-codel", "IDH-O", "IDH-A")
md$subtype <- factor(md$subtype, levels=c("IDH-O", "IDH-A"))

md_trim <- md %>% 
  dplyr::select(SampleID, case_barcode, subtype, care_id, patient_id, timepoint) %>% 
  distinct()
rownames(md_trim) <- NULL

cell_prop <- md %>% 
  group_by(SampleID, CellType_final) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  dplyr::select(SampleID, CellType_final, counts, freq) %>% 
  complete(SampleID, CellType_final,
           fill = list(counts = 0, freq = 0)) %>%
  distinct() 

top_cell_prop <- cell_prop %>% 
  group_by(CellType_final) %>% 
  summarise(median_freq = median(freq)) %>% 
  arrange(desc(median_freq))

cellnum_by_time <- md %>% 
  group_by(case_barcode, timepoint, subtype) %>% 
  summarise(cell_number = n()) %>% 
  ungroup() %>% 
  mutate(timepoint = recode(timepoint, `T1` = "I",
                            `T2` = "R"))

feature_by_time <- md %>% 
  group_by(case_barcode, timepoint, subtype) %>% 
  summarise(median_nFeature = median(nFeature_RNA)) %>% 
  ungroup() %>% 
  mutate(timepoint = recode(timepoint, `T1` = "I",
                            `T2` = "R"))

mt_by_time <- md %>% 
  group_by(case_barcode, timepoint, subtype) %>% 
  summarise(median_percent_mt = median(percent.mt)) %>% 
  ungroup() %>% 
  mutate(timepoint = recode(timepoint, `T1` = "I",
                            `T2` = "R"))

all_md <- cellnum_by_time %>% 
  inner_join(mt_by_time, by=c("case_barcode", "timepoint", "subtype")) %>% 
  inner_join(feature_by_time, by=c("case_barcode", "timepoint", "subtype")) %>% 
  pivot_longer(cols=c(cell_number, median_percent_mt, median_nFeature), names_to = "metric", values_to = "values")

ggplot(all_md %>% filter(timepoint!="T3"), aes(x = timepoint, y = values)) + 
  geom_line(aes(group=case_barcode), color="gray70", linetype=1, alpha = 0.5) +
  geom_point(aes(color=subtype)) +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size =4, label="p.format") +
  scale_color_manual(values=c("IDH-A" = "#67A9CF",
                              "IDH-O"="#EF8A62"),
                     labels=c("IDH-O", "IDH-A")) +
  labs(x="Time point", y="Detected gene number (nFeature RNA)", color="IDHmut subtype") +
  theme(strip.background = element_blank()) +
  facet_wrap(.~metric, scales="free") 


pdf(paste0(fig_dir, "seurat_mitochondrial_rna_paired_wilcox.pdf"), width = 2, height = 2.75, useDingbats = FALSE)
ggplot(mt_by_time %>% filter(timepoint!="T3"), aes(x = timepoint, y = median_percent_mt)) + 
  geom_line(aes(group=case_barcode), color="gray70", linetype=1, alpha = 0.5) +
  geom_point() +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size =4, label="p.format") +
  labs(x="Time point", y="Median percent counts\nfrom mitochondrial RNA") +
  theme(strip.background = element_blank()) +
  ylim(0, 5)
dev.off()

pdf(paste0(fig_dir, "seurat_gene_count_paired_wilcox.pdf"), width = 2, height = 2.75, useDingbats = FALSE, bg = "transparent")
ggplot(feature_by_time %>% filter(timepoint!="T3"), aes(x = timepoint, y = median_nFeature)) + 
  geom_line(aes(group=case_barcode), color="gray70", linetype=1, alpha = 0.5) +
  geom_point() +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size =4, label="p.format") +
  labs(x="Time point", y="Median gene count\n(nFeature RNA)") +
  theme(strip.background = element_blank())
dev.off()

pdf(paste0(fig_dir, "seurat_cell_number_paired_wilcox.pdf"), width = 2, height = 2.75, useDingbats = FALSE)
ggplot(cellnum_by_time %>% filter(timepoint!="T3"), aes(x = timepoint, y = cell_number)) + 
  geom_line(aes(group=case_barcode), color="gray70", linetype=1, alpha = 0.5) +
  geom_point() +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size =4, label="p.format") +
  labs(x="Time point", y="Median cell number") +
  theme(strip.background = element_blank()) +
  ylim(0, 12500)
dev.off()
  

### END ### 