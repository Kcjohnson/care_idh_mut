##################################
# Inspect the longitudinal differentially accessible genes across malignant cell states
# Author: Kevin Johnson
# Date Updated: 2026.04.09
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

#patient_list_results_filt <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/archr_care_cell_type_longitudinal_gene_score_min40cells_perstate_20240608.txt", sep="\t", header = TRUE)

mapping_df <- read.table("/vast/palmer/pi/verhaak/shared/monitor/metadata/glass_sn_bulk_master_mapping_file_20230729.txt", sep = "\t", header = T)
mapping_df_trim <- mapping_df %>% 
  mutate(patient_id_linker = sapply(strsplit(multiome_id, "-"), "[[", 1)) %>% 
  dplyr::select(patient_id, patient_id_linker) %>% 
  distinct()

patient_list_results_filt <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/archr_malignant_cell_state_longitudinal_gene_score_min25cells_perstate.txt", sep="\t", header = TRUE)
patient_list_results_filt <- patient_list_results_filt %>% 
  inner_join(mapping_df_trim, by=c("patient_id.i." = "patient_id_linker"))
  
res <- patient_list_results_filt %>% 
  left_join(patient_md, by="patient_id")
res$down_genes <- ifelse(res$down_genes!=0, res$down_genes*-1, res$down_genes)

patients_to_include <- patient_list_results_filt %>% 
  filter(use_cells>100, bgd_cells>100)


malignant_res <- res %>% 
  dplyr::select(patient_id, cell_state = cell_state.j., grade_change_t1t2, acquired_genetic_alt_t1t2, idh_codel_subtype, up_genes, down_genes) %>% 
  pivot_longer(cols= c(up_genes:down_genes),
               names_to = "Type",
               values_to = "Values")  %>% 
  mutate(cell_state = recode(cell_state, `AC` = "AC-like",
                             `OPC` = "OPC-like",
                             `Undifferentiated` = 'Undifferentiated'))
malignant_res$acquired_genetic_alt_t1t2 <- ifelse(is.na(malignant_res$acquired_genetic_alt_t1t2), "No acquired alt.", "Acquired genetic alt.")
malignant_res$Type <- ifelse(malignant_res$Type=="up_genes", "Increased accessibility", "Decreased accessibility")
malignant_res$Type <- factor(malignant_res$Type, levels=c("Increased accessibility", "Decreased accessibility"))


#pdf(paste0(fig_dir, "longitudinal_dag_barplot_by_genetics_20240616.pdf"), width = 7, height = 5, useDingbats = FALSE, bg = "transparent")
ggplot(malignant_res %>% filter(cell_state%in%c("AC-like", "OPC-like", "Undifferentiated"),
                                patient_id%in%patients_to_include$patient_id), aes(x = patient_id, y = Values, fill = Type)) +
  geom_bar(stat = 'identity', position = 'identity', color = 'black') +
  labs(x = 'Patients', y = 'Longitudinal differential accessibility\nstate-conrolled', fill="Longitudinal change") +
  scale_fill_manual(values = c("Decreased accessibility" = 'dodgerblue', "Increased accessibility" = 'red')) +
  plot_theme +
  facet_grid(cell_state~acquired_genetic_alt_t1t2, scales="free") +
  ylim(-4500,2000) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 
#dev.off()


# Compare the number of differentially accessible genes
malignant_res_compare <- malignant_res %>% filter(cell_state%in%c("AC-like", "OPC-like", "Undifferentiated"),
                                                  patient_id%in%patients_to_include$patient_id) %>% 
  mutate(all_dag = ifelse(Type=="Decreased accessibility", -1*Values, Values)) %>% 
  group_by(patient_id, cell_state, acquired_genetic_alt_t1t2) %>% 
  summarise(total_dag = sum(all_dag))

#pdf(paste0(fig_dir, "longitudinal_dag_by_genetics_20240616.pdf"), width = 4, height = 5, useDingbats = FALSE, bg = "transparent")
ggplot(malignant_res_compare, aes(x=acquired_genetic_alt_t1t2, y=total_dag)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color=cell_state)) +
  stat_compare_means(method = "wilcox", paired = FALSE, size =4, label="p.format") +
  plot_theme + 
  stat_n_text() + 
  scale_color_manual(values=c("Cycling" = "#6BAED6", 
                              "AC-like" = "#AA2756", 
                              "MES-like"="#F77D58",
                              "NPC-like" = "#7fbf7b",
                              "OPC-like"="#E8F5A3",
                              "Undifferentiated" = "gray90")) +
  labs(x="", y="Longitudinal differentially\naccessible genes", color="Malignant\ncell state") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 
#dev.off()


### END ###