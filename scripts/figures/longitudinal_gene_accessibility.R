##################################
# Inspect the longitudinal differentially accessible genes across malignant cell states
# Author: Kevin Johnson
##################################

library(tidyverse)
library(ggpubr)
library(EnvStats)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
setwd(proj_dir)
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")


# Load sample-specific metadata
sample_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240608.txt", sep="\t", header = TRUE)
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt", sep="\t", header = TRUE)

# Linker for snATAC samples
atac_mapping_df <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/misc/verhaak_lab_patient_id_linker_for_atac.tsv", sep = "\t", header = T)

patient_list_results_filt <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/archr_malignant_cell_state_longitudinal_gene_score_min25cells_perstate.txt", sep="\t", header = TRUE)
patient_list_results_filt <- patient_list_results_filt %>% 
  inner_join(atac_mapping_df, by=c("patient_id.i." = "patient_id_linker"))
  
res <- patient_list_results_filt %>% 
  left_join(patient_md, by="patient_id")
res$down_genes <- ifelse(res$down_genes!=0, res$down_genes*-1, res$down_genes)

# Restrict to a minimum of 100 cells per case vs control group
comparisons_to_keep <- patient_list_results_filt %>% 
  filter(use_cells>=100, bgd_cells>=100)

malignant_atac_res <- res %>% 
  # Keep analyses where at least 100 cells were measured per time point
  filter(use_cells>=100, bgd_cells>=100) %>% 
  dplyr::select(patient_id, cell_state = cell_state.j., acquired_genetic_alt_t1t2, up_genes, down_genes) %>% 
  pivot_longer(cols= c(up_genes:down_genes),
               names_to = "Type",
               values_to = "Values")  %>% 
  # Relabel features
  mutate(cell_state = recode(cell_state, `AC` = "AC-like",
                             `OPC` = "OPC-like",
                             `NPC` = "NPC-like",
                             `MES` = "MES-like",
                             `Undifferentiated` = 'Undifferentiated'))  
malignant_atac_res$acquired_genetic_alt_t1t2 <- ifelse(is.na(malignant_atac_res$acquired_genetic_alt_t1t2), "Did not acquire", "Acquired genetic alt.")
malignant_atac_res$Type <- ifelse(malignant_atac_res$Type=="up_genes", "Increased accessibility", "Decreased accessibility")
malignant_atac_res$Type <- factor(malignant_atac_res$Type, levels=c("Increased accessibility", "Decreased accessibility"))


pdf(paste0(fig_dir, "longitudinal_dag_barplot_by_genetics.pdf"), width = 3.5, height = 2.75, useDingbats = FALSE, bg = "transparent")
ggplot(malignant_atac_res %>% 
         # Remove the one MES-like change since plotting space is limited
         filter(cell_state%in%c("AC-like", "OPC-like", "NPC-like", "Undifferentiated")) %>% 
         mutate(
           cell_state = recode(cell_state, `Undifferentiated` = "Undiff."),
           Type = recode(Type, 
                         `Decreased accessibility` = "Decreased", 
                         `Increased accessibility` = "Increased")
         ), aes(x = patient_id, y = Values, fill = Type)) +
  geom_bar(stat = 'identity', position = 'identity', color = 'black') +
  labs(x = 'Patients', y = 'Longitudinal differentially accessible\ngenes state-controlled', fill="") +
  scale_fill_manual(values = c("Decreased" = 'dodgerblue', "Increased" = 'red')) +
  plot_theme +
  facet_grid(cell_state~acquired_genetic_alt_t1t2, scales="free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
                legend.position = "right",
                legend.key.size = unit(0.2, "cm"),
                legend.text = element_text(size = 6),
                legend.title = element_text(size = 6),
                legend.spacing.y = unit(0.1, "cm"),
                legend.margin = margin(0, 0, 0, 0),
                legend.box.margin = margin(0, 0, 0, -5))
dev.off()


# Compare the number of longitudinal differentially accessible genes across these recurrence associated genetic alterations
malignant_res_compare <- malignant_atac_res %>% filter(cell_state%in%c("AC-like", "OPC-like", "Undifferentiated","MES-like", "NPC-like")) %>% 
  mutate(all_dag = ifelse(Type=="Decreased accessibility", -1*Values, Values)) %>% 
  group_by(patient_id, cell_state, acquired_genetic_alt_t1t2) %>% 
  summarise(total_dag = sum(all_dag))

pdf(paste0(fig_dir, "longitudinal_dag_by_genetics_wilcoxon.pdf"), width = 2, height = 2.75, useDingbats = FALSE, bg = "transparent")
ggplot(malignant_res_compare %>% 
         mutate(cell_state = recode(cell_state, `Undifferentiated` = "Undiff.")), aes(x=acquired_genetic_alt_t1t2, y=total_dag)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color=cell_state), size = 0.5) +
  geom_point(aes(color=cell_state), position = position_jitter(width = 0.1, seed = 42), size = 0.5) +
  stat_compare_means(method = "wilcox", paired = FALSE, size =2.25, label="p.format") +
  plot_theme + 
  stat_n_text(size = 2.25) + 
  scale_color_manual(values=c("Cycling" = "#6BAED6", 
                              "AC-like" = "#AA2756", 
                              "MES-like"="#F77D58",
                              "NPC-like" = "#7fbf7b",
                              "OPC-like"="#E8F5A3",
                              "Undiff." = "gray90")) +
  labs(x="", y="Longitudinal differentially\naccessible genes", color="Malignant\ncell state") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "right",
        legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.spacing.y = unit(0.1, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, -5)) 
dev.off()

