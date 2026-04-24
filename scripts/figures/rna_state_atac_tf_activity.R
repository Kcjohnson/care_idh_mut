##################################
# Visualize the snATAC-derived transcription factor motif activity scores across RNA cell states
# Author: Kevin Johnson
# Date Updated: 2026.04.10
##################################

library(tidyverse)
library(ggpubr)
library(EnvStats)
library(purrr)
library(broom)
library(scales)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
setwd(proj_dir)
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")

md_all <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20260320.txt", sep = "\t", row.names = 1, header = TRUE)
caremut_md_case <- md_all %>% 
  mutate(patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2))) %>% 
  dplyr::select(SampleID, lab, sample_barcode, idh_codel_subtype, care_id, patient_id, timepoint) %>% 
  distinct()
rownames(caremut_md_case) <- NULL


# Curated genomic information at the patient-level (changes) and sample-level.
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt", sep="\t", header = TRUE)
sample_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240608.txt", sep="\t", header = TRUE)

sample_grade <- sample_md %>% 
  mutate(grade = paste0("G", grade_num)) %>% 
  dplyr::select(care_id, grade, idh_codel_subtype) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."))

# Read in the cell state classification based on p-value assignment.
care_state_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/scoring/caremut_malignant_cell_state_assignment.txt", sep="\t", header = TRUE)

# Relabel some of the features.
care_state_md <- care_state_md %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."),
         patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2)), 
         cell_state = recode(State, `MP_AC1_MUT` = "AC-like",
                             `MP_OPC_MUT` = "OPC-like",
                             `MP_NPC_MUT` = "NPC-like",
                             `MP_MES_MUT` = "MES-like",
                             `MP_AC2_MUT` = "AC-like",
                             "Undifferentiated" = "Undifferentiated")) 

# Examination of TF activity from ATAC data
tf_motif_activity <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/atac/archr_care_malignant_state_tf_motif_activity_zscore.txt", row.names = 1, header = TRUE, sep = "\t")

### ### ### ### ### ### ### ### ###
### Compare ATAC TF activity with RNA states and metaprogram expression
### ### ### ### ### ### ### ### ###
malignant_tf_annot <- care_state_md %>% 
  inner_join(tf_motif_activity, by=c("CellID")) %>% 
  inner_join(sample_grade, by="care_id")


# Calculate the correlations between each RNA metaprogram and each ATAC TF Zscores
mp_score_rna <- malignant_tf_annot[, c(3:8)]
tf_zscores_atac <- malignant_tf_annot[, c(17:886)]

correlations <- sapply(mp_score_rna, function(individual_mp) {
  cor_values <- cor(tf_zscores_atac, individual_mp, method = "s")
  return(cor_values)
})

rownames(correlations) <- colnames(tf_zscores_atac)
module_tf_zscore_cor <- as.data.frame(correlations)
module_tf_zscore_cor$TF <- rownames(module_tf_zscore_cor)
module_tf_zscore_cor_rna <- module_tf_zscore_cor[, c(1:6, 7)]

# Filter rows where any of the specified columns have a correlation coefficient greater than 0.4
# Ignore MP_AC2
columns_to_check <- c(1,3:6)
module_tf_zscore_cor_rna_filt <- module_tf_zscore_cor_rna[apply(module_tf_zscore_cor_rna[, columns_to_check], 1, function(row) any(row > 0.4)), ]

module_tf_zscore_cor_rna_filt_long <- module_tf_zscore_cor_rna_filt %>% 
  dplyr::select(-MP_AC2_MUT) %>% 
  pivot_longer(cols = c(MP_AC1_MUT, MP_MES_MUT:MP_CC_MUT), names_to = "metaprogram", values_to = "correlation") %>% 
  mutate(MP = recode(metaprogram, `MP_AC1_MUT` = "AC MP",
                     `MP_OPC_MUT` = "OPC MP",
                     `MP_MES_MUT` = "MES MP",
                     `MP_NPC_MUT` = "NPC MP",
                     `MP_CC_MUT` = "Cell Cycle"),
         TF = sapply(strsplit(TF, "_"), "[[", 1)) 


# For plotting purposes, it can be beneficial to simply plot the top hits
top_10_cor <- module_tf_zscore_cor_rna_filt_long %>%
  group_by(MP) %>%
  top_n(10, correlation) 


terms_to_plot <- unique(top_10_cor$TF)
Res_plot <- module_tf_zscore_cor_rna_filt_long %>%
  filter(TF%in%terms_to_plot)

mp_order <- c("OPC MP", "NPC MP", "MES MP", "AC MP", "Cell Cycle")

Res_plot$MP <- factor(Res_plot$MP, levels=mp_order)
Res_plot$TF <- as.factor(Res_plot$TF)

max_correlations <- Res_plot %>%
  group_by(TF) %>%
  summarize(max_correlation = max(correlation)) %>%
  arrange(desc(max_correlation))

# Reorder the TF factor based on the highest correlation values
Res_plot$TF <- factor(Res_plot$TF, levels =rev(max_correlations$TF))

pdf(paste0(fig_dir, "malignant_mp_scores_atac_tf_activity_correlation.pdf"), width = 2.25, height = 2.25, bg = "transparent", useDingbats = FALSE)
ggplot(data = Res_plot %>% 
         filter(TF%in%c("TCF12", "ASCL1", "NHLH1", "STAT3", "FOSL1", "JUND")), aes(y=TF,x=MP, fill = correlation)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "dodgerblue", high = "#FF0000", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="ATAC-RNA per cell\ncorrelation") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 7, hjust = 1), 
        text=element_text(size=7),
        axis.text.y = element_text(size=7),
        legend.position="top") +
  labs(y= "ATAC - TF activity", x= "RNA - metaprogram score") +
  guides(
    fill = guide_colorbar(
      title.theme = element_text(size = 6),
      label.theme = element_text(size = 6),
      barwidth  = unit(0.6, "in"),
      barheight = unit(0.08, "in"),
      ticks.colour = "black",
      frame.colour = "black"
    )
  )
dev.off()


###
### Compare TF activities across cell states
###
tf_zscores_atac_long <- malignant_tf_annot %>% 
  pivot_longer(cols=c(TFAP2B_1:TBX22_870), names_to = "TF", values_to = "z_score") 

tf_zscores_atac_long_summary <- tf_zscores_atac_long %>% 
  dplyr::group_by(SampleID, cell_state, TF) %>% 
  dplyr::summarise(counts = n()) %>% 
  # Set a minimum threshold for how many minimum cells per state, per sample need to be considered
  filter(counts >= 20) %>% 
  ungroup() %>% 
  mutate(sample_state = paste0(SampleID, "_", cell_state))

tf_zscores_atac_long_summary_filt <- tf_zscores_atac_long %>% 
  mutate(sample_state = paste0(SampleID, "_", cell_state)) %>% 
  filter(sample_state%in%tf_zscores_atac_long_summary$sample_state) %>% 
  dplyr::group_by(care_id, cell_state, TF) %>% 
  dplyr::summarise(median_zscore = median(z_score)) %>% 
  ungroup() %>% 
  inner_join(sample_grade, by="care_id") %>% 
  mutate(cell_state = recode(cell_state, `Undifferentiated` = "Undiff."))

tf_zscores_atac_long_summary_filt$cell_state <- factor(tf_zscores_atac_long_summary_filt$cell_state, levels= c("AC-like", "MES-like", "Undiff.", "NPC-like", "OPC-like"))
tf_zscores_atac_long_summary_filt$idh_codel_subtype <- factor(tf_zscores_atac_long_summary_filt$idh_codel_subtype, levels= c("Oligo.", "Astro."))
tf_zscores_atac_long_summary_filt$TF <- sapply(strsplit(tf_zscores_atac_long_summary_filt$TF, "_"), "[[", 1)

# Compare TF activities across cell states
zscore_res <- tf_zscores_atac_long_summary_filt %>%
  filter(TF%in%Res_plot$TF) %>% 
  group_by(TF) %>%
  nest() %>%
  mutate(
    kruskal_res = map(data, ~ kruskal.test(.x$median_zscore ~ .x$cell_state) %>% tidy()),
    n = map_int(data, nrow)
  ) %>%
  unnest(kruskal_res) %>%
  as.data.frame() %>% 
  mutate(adj_p_value = p.adjust(p.value, method = "fdr"))

tf_activity_avg <- tf_zscores_atac_long_summary_filt %>% 
  filter(TF%in%Res_plot$TF) %>% 
  group_by(TF,cell_state) %>% 
  summarise(avg_median_zscore = mean(median_zscore))

# Calculate a metric on which to sort the activities.
max_activity <- tf_activity_avg %>%
  group_by(TF) %>%
  summarize(max_zscore = max(avg_median_zscore)) %>%
  arrange(desc(max_zscore))

tf_activity_avg$TF <- factor(tf_activity_avg$TF, levels =rev(max_activity$TF))

#pdf(paste0(fig_dir, "tf_activity_malignant_cell_state_min20_cells_top_tfs.pdf"), width = 5, height = 6, useDingbats = FALSE, bg = "transparent")
ggplot(data = tf_activity_avg, aes(x=cell_state,y=TF, fill = avg_median_zscore)) +
  geom_tile(color = "white")+
  scale_fill_viridis_c(option = "plasma", name="Mean per sample\nTF activity z-score", limits=c(-5,5), oob=squish) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1), 
        text=element_text(size=10),
        axis.text.y = element_text(size=10),
        legend.position="top")+
  labs(x= "Malignant state", y= "Transcription factor") 
#dev.off()


pdf(paste0(fig_dir, "tf_activity_malignant_cell_state_min20_cells_ascl1_jund.pdf"), width = 3, height = 2.75, useDingbats = FALSE, bg = "transparent")
ggplot(tf_zscores_atac_long_summary_filt %>% filter(TF%in%c("JUND","ASCL1")), aes(x=cell_state, y=median_zscore)) +
  geom_boxplot(aes(fill=cell_state), outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1, seed = 42), 
             size = 0.8, alpha = 0.6) +
  facet_grid(TF~., scales="free") +
  labs(x="Malignant cell state",  y="TF activity\n(median deviation z-score)") +
  plot_theme +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undiff." = "gray90")) +
  stat_compare_means(method="kruskal", label="p.format", size = 2.25) +
  stat_n_text(size = 2.25) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(size = 7), 
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        strip.text = element_text(size = 7)) +
  guides(fill=FALSE) 
dev.off()

pdf(paste0(fig_dir, "tf_activity_malignant_cell_state_min20_cells_ascl1_jund_grade_combined_astro.pdf"), width = 3, height = 2.25, useDingbats = FALSE, bg = "transparent")
ggplot(tf_zscores_atac_long_summary_filt %>% filter(TF%in%c("JUND"), idh_codel_subtype=="Astro.", cell_state%in%c("AC-like", "Undiff.", "OPC-like")) %>% 
         mutate(comb_grade = ifelse(grade=="G4", "G4", "G2/G3")), aes(x=comb_grade, y=median_zscore)) +
  geom_boxplot(aes(fill=cell_state), outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1, seed = 42), 
             size = 0.8, alpha = 0.6) +
  facet_grid(.~cell_state, scales="free") +
  labs(x="Tumor grade, Astro.",  y="JUND TF activity\n(median deviation z-score)") +
  plot_theme +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undiff." = "gray90")) +
  stat_compare_means(method="wilcox", label="p.format", size = 2.25) +
  stat_n_text(size = 2.25) +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(size = 7), 
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        strip.text = element_text(size = 7)) +
  guides(fill=FALSE) 
dev.off()



### ### ### ### ### ### ### ### 
# Part 3: State controlled longitudinal analysis
### ### ### ### ### ### ### ###
tf_zscores_atac_long <- malignant_tf_annot %>% 
  pivot_longer(cols=c(TFAP2B_1:TBX22_870), names_to = "TF", values_to = "z_score") 

# Filter to the samples with at least 20 cells per state per TF.
tf_zscores_atac_long_summary <- tf_zscores_atac_long %>% 
  dplyr::group_by(SampleID, cell_state, TF) %>% 
  dplyr::summarise(counts = n()) %>% 
  filter(counts >= 20) %>% 
  ungroup() %>% 
  mutate(sample_state = paste0(SampleID, "_", cell_state))

# Use that filtering to calculate the median_zscore per state per TF.
tf_zscores_atac_long_summary_filt <- tf_zscores_atac_long %>% 
  mutate(sample_state = paste0(SampleID, "_", cell_state)) %>% 
  filter(sample_state%in%tf_zscores_atac_long_summary$sample_state) %>% 
  dplyr::group_by(SampleID, cell_state, TF) %>% 
  dplyr::summarise(median_zscore = median(z_score)) %>% 
  ungroup()


# Combine with T1 and T2 nomenclature and compare.
tf_zscores_df_wide <- tf_zscores_atac_long_summary_filt %>% 
  inner_join(caremut_md_case, by="SampleID") %>% 
  dplyr::select(patient_id, timepoint, cell_state, TF, median_zscore) %>% 
  filter(timepoint!="T3") %>% 
  pivot_wider(names_from = timepoint, values_from = median_zscore) %>% 
  filter(!is.na(T1), !is.na(T2)) 


# Performed a paired Wilcoxon test across each TF per cell state
longitudinal_zscore_res <- tf_zscores_df_wide %>%
  group_by(TF, cell_state) %>%
  nest() %>%
  mutate(
    paired_wilcox = map(data, ~ wilcox.test(.x$T1, .x$T2, paired = TRUE) %>% tidy()),
    median_diff = map_dbl(data, ~ median(.x$T2 - .x$T1)),
    n = map_int(data, nrow)
  ) %>%
  unnest(paired_wilcox) %>%
  as.data.frame() %>% 
  mutate(direction = if_else(median_diff > 0, "increase", "decrease"),
         adj_p_value = p.adjust(p.value, method = "fdr"))


sum(longitudinal_zscore_res$p.value<0.05)

longitudinal_zscore_res_filt <- longitudinal_zscore_res %>% 
  filter(p.value< 0.05, direction=="increase")
table(longitudinal_zscore_res_filt$TF)

caremut_md_case_filt <- caremut_md_case %>% 
  dplyr::select(patient_id, idh_codel_subtype) %>% 
  distinct()

tf_zscores_df_long <- tf_zscores_df_wide %>% 
  pivot_longer(cols=c(T1, T2), names_to = "timepoint") %>% 
  mutate(TF = sapply(strsplit(TF, "_"), "[[", 1)) %>% 
  filter(TF%in%c("JUND")) %>% 
  inner_join(caremut_md_case_filt, by="patient_id")

tf_zscores_df_long$cell_state <- factor(tf_zscores_df_long$cell_state, levels= c("AC-like", "MES-like", "Undifferentiated", "NPC-like", "OPC-like"))


pdf(paste0(fig_dir, "tf_activity_longitudinal_malignant_cell_state_min20_cells_jund_opc.pdf"), width = 2, height = 2.25, useDingbats = FALSE, bg = "transparent")
ggplot(tf_zscores_df_long %>% 
         mutate(timepoint = recode(timepoint, `T1` = "Init.",
                                   `T2` = "Recur.")) %>% 
         filter(cell_state%in%c("OPC-like")), aes(x = as.factor(timepoint), y = value)) + 
  geom_boxplot(aes(fill=cell_state), outlier.shape = NA) +
  geom_line(aes(group=patient_id), color="gray70", linetype=2) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  facet_grid(.~cell_state, scales="free") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  labs(x="Longitudinal state-controlled pairs\n(min. 20 ATAC nuclei per state)", y="JUND activity\n(median deviation z-score)") +
  plot_theme +
  theme(text = element_text(size = 7), 
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        strip.text = element_text(size = 7)) +
  stat_n_text(size = 2.25) +
  theme(legend.position = "none") 
dev.off()


# All cell states grouped together
ggplot(tf_zscores_df_long %>% 
         mutate(patient_state = paste0(patient_id, "_", cell_state)) %>% 
         filter(cell_state%in%c("AC-like", "OPC-like", "Undifferentiated")), aes(x = as.factor(timepoint), y = value)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group=patient_state, color=cell_state), linetype=2) +
  geom_point() +
  scale_color_manual(values=c("Cycling" = "#6BAED6", 
                              "AC-like" = "#AA2756", 
                              "MES-like"="#F77D58",
                              "NPC-like" = "#7fbf7b",
                              "OPC-like"="#E8F5A3",
                              "Undifferentiated" = "gray90")) +
  stat_compare_means(method = "wilcox", paired = TRUE, size =3, label="p.format") +
  labs(x="Longitudinal state-controlled pairs\n(min. 20 ATAC nuclei per state)", y="JUND activity\n(median deviation z-score)") +
  plot_theme +
  stat_n_text() +
  theme(legend.position = "none")


### END ###