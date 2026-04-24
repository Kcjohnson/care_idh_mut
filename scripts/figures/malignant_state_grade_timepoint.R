##################################
# Examine distributions of the CARE IDH-mutant cell states across tumor grade, tumor types, and time points
# Author: Kevin Johnson
##################################

# Reproduce the following figures: Fig3a, Fig3b, Fig3c, Fig3d, EDF7a, EDF8a, EDF8b, EDF8c

library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(EnvStats)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
source(paste0(proj_dir, "/scripts/utils/plot_theme.R"))
setwd(proj_dir)

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
  

md_trim <- care_state_md %>% 
  dplyr::select(SampleID, case_barcode, idh_codel_subtype, care_id, patient_id, timepoint) %>% 
  distinct()

# Curated genomic information at the patient-level (changes) and sample-level.
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt", sep="\t", header = TRUE)
patient_md <- patient_md %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                  `IDHmut-noncodel` = "Astro."))

sample_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240608.txt", sep="\t", header = TRUE)
sample_md <- sample_md %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."))

sample_grade <- sample_md %>%
  mutate(grade = paste0("G", grade_num)) %>% 
  dplyr::select(care_id, grade)

### ### ### ### ### ### ### ### ### ###
# Tumor type, grade, and time point differences
### ### ### ### ### ### ### ### ### ###
# cell_state reflects AC1 and AC2 collapsed into AC-like
care_state_freq <- care_state_md %>% 
  group_by(care_id, cell_state) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts/sum(counts)) %>% 
  ungroup() %>% 
  # If a state is not detected in a sample, add in zero values to indicate they were not observed
  complete(care_id, cell_state,
           fill = list(counts = 0, freq = 0)) %>%
  inner_join(md_trim, by="care_id") 

# Confirm that all sample sum 1 and that all cell types are measured.
care_state_freq %>%
  group_by(care_id) %>%
  summarise(freq_sum = sum(freq)) %>%
  mutate(pass = abs(freq_sum - 1) < 1e-10) %>%
  { if (all(.$pass)) {
    message("PASS: All samples frequencies sum to 1")
  } else {
    failing <- filter(., !pass)
    message("FAIL: ", nrow(failing), " samples do not sum to 1:")
    print(failing)
  }
  }
# Visual sanity check
ggplot(care_state_freq, aes(x = care_id, y=freq*100, fill = cell_state)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  plot_theme

care_state_freq_cc <- care_state_md %>% 
  group_by(care_id, isCC) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts/sum(counts)) %>% 
  ungroup() %>% 
  # If a state is not detected in a sample, add in zero values to indicate they were not observed
  complete(care_id, isCC,
           fill = list(counts = 0, freq = 0)) %>%
  filter(isCC==TRUE) %>% 
  mutate(cell_state = "Cycling") %>% 
  inner_join(md_trim, by="care_id") %>% 
  dplyr::select(-isCC)

# Combine malignant state and cycling percentages into one and relabel some of the features.
mal_care_md_state_summary <- care_state_freq %>% 
  bind_rows(care_state_freq_cc) %>% 
  inner_join(patient_md, by=c("case_barcode", "patient_id", "idh_codel_subtype"))  %>% 
  mutate(sample_time = recode(timepoint, "T1" = "Init.",
                       "T2" = "Recur.",
                       "T3" = "2nd Recur.")) %>% 
  inner_join(sample_grade, by="care_id")

mal_care_md_state_summary$cell_state <- factor(mal_care_md_state_summary$cell_state, levels=c("AC-like", "MES-like", "NPC-like", "OPC-like",  "Undifferentiated", "Cycling"))
mal_care_md_state_summary$idh_codel_subtype <- factor(mal_care_md_state_summary$idh_codel_subtype, levels=c("Oligo.", "Astro."))

# Confirm that each state has 74 observations
all(table(mal_care_md_state_summary$cell_state)==74)

# TUMOR TYPE
pdf(paste0(fig_dir, "edf7a_care_tumor_type_malignant_abundance.pdf"), width = 3.75, height = 2.75, useDingbats = FALSE, bg = "transparent")
ggplot(mal_care_md_state_summary, aes(x = idh_codel_subtype, y = freq*100)) + 
  geom_boxplot(aes(fill=cell_state), outlier.shape = NA) +
  geom_point(size = 0.6) +
  plot_theme +
  theme(legend.position = "none",
        axis.text.x= element_text(angle=45,hjust=1)) +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  stat_compare_means(method = "wilcox", paired = FALSE, size = 2.25, label="p.format") +
  facet_grid(.~cell_state, scales="free") +
  labs(x="Tumor type", y="snRNA malignant cell abundance (%)") +
  stat_n_text(size = 2.25)
dev.off()


# TUMOR GRADE n = 74
# KRUSKAL TEST for Astrocytoma
pdf(paste0(fig_dir, "fig3b_grade_malignant_abundance_kruskal.pdf"), width = 4, height = 3, useDingbats = FALSE, bg = "transparent")
ggplot(mal_care_md_state_summary %>% filter(cell_state%in%c("Undifferentiated", "MES-like", "AC-like", "Cycling")), aes(x = as.factor(grade), y = freq*100)) + 
  geom_boxplot(aes(fill=cell_state), outlier.shape = NA) +
  geom_point(size = 0.6) +
  plot_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  stat_compare_means(method = "kruskal", paired = FALSE, size = 2.25, label="p.format") +
  facet_grid(idh_codel_subtype~cell_state, scales="free") +
  labs(x="Tumor grade", y="Malignant cell abundance (%)") +
  stat_n_text(size = 2.25) 
dev.off()

# WILCOXON TEST for Oligodendroglioma
pdf(paste0(fig_dir, "fig3b_grade_malignant_abundance_wilcoxon.pdf"), width = 4, height = 3, useDingbats = FALSE, bg = "transparent")
ggplot(mal_care_md_state_summary %>% filter(cell_state%in%c("Undifferentiated", "MES-like", "AC-like", "Cycling")), aes(x = as.factor(grade), y = freq*100)) + 
  geom_boxplot(aes(fill=cell_state), outlier.shape = NA) +
  geom_point(size = 0.6) +
  plot_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  stat_compare_means(method = "wilcox", paired = FALSE, size = 2.25,  label="p.format") +
  facet_grid(idh_codel_subtype~cell_state, scales="free") +
  labs(x="Tumor grade", y="Malignant cell abundance (%)") +
  stat_n_text(size = 2.25) 
dev.off()


# INITIAL VS. RECURRENCE  n = 35 patients
pdf(paste0(fig_dir, "edf8a_init_v_recur_malignant_abundance.pdf"), width = 3.75, height = 3, useDingbats = FALSE, bg = "transparent")
ggplot(mal_care_md_state_summary %>% 
         filter(timepoint!="T3"), aes(x = sample_time, y = freq*100)) + 
  geom_boxplot(outlier.shape = NA, aes(fill=cell_state)) +
  geom_line(aes(group=patient_id), color="gray85", linetype=2, size = 0.5) +
  geom_point(size = 0.5) +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  facet_grid(.~cell_state, scales="free") +
  labs(x="Matched longitudinal analysis", y="Malignant cell abundance (%)", color="Subtype") +
  stat_n_text(size=2.25) +
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

mal_care_md_state_summary$acquired_genetic_alt_t1t2 <- ifelse(is.na(mal_care_md_state_summary$acquired_genetic_alt_t1t2), "Did not acquire", "Acquired genetic alt.")

pdf(paste0(fig_dir, "fig3d_genetic_alt_init_v_recur_malignant_abundance.pdf"), width = 4, height = 3, useDingbats = FALSE, bg = "transparent")
ggplot(mal_care_md_state_summary %>% 
         filter(timepoint!="T3", cell_state%in%c("AC-like", "MES-like", "Undifferentiated", "Cycling")), aes(x = sample_time, y = freq*100)) + 
  geom_boxplot(outlier.shape = NA, aes(fill=cell_state)) +
  geom_line(aes(group=patient_id), color="gray85", linetype=2, size = 0.5) +
  geom_point(size = 0.5) +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  facet_grid(acquired_genetic_alt_t1t2~cell_state, scales="free") +
  labs(x="Matched longitudinal analysis", y="Malignant cell abundance (%)", color="Subtype") +
  stat_n_text(size=2.25) +
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# INITIAL VS RECURRENCE stratified by tumor type.
pdf(paste0(fig_dir, "edf8b_init_v_recur_tumor_type_malignant_abundance.pdf"), width = 3.75, height = 3, useDingbats = FALSE, bg = "transparent")
ggplot(mal_care_md_state_summary %>% 
         filter(timepoint!="T3", cell_state%in%c("AC-like", "MES-like", "Undifferentiated", "Cycling")), aes(x = sample_time, y = freq*100)) + 
  geom_boxplot(outlier.shape = NA, aes(fill=cell_state)) +
  geom_line(aes(group=patient_id), color="gray85", linetype=2, size = 0.5) +
  geom_point(size = 0.5) +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  facet_grid(idh_codel_subtype~cell_state, scales="free") +
  labs(x="Matched longitudinal analysis", y="Malignant cell abundance (%)", color="Subtype") +
  stat_n_text(size=2.25) +
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# INITIAL VS RECURRENCE. Initial sample = first resection. Astrocytomas.
true_primaries <- patient_md %>% 
  filter(t1_clinical_time=="t1_primary")

pdf(paste0(fig_dir, "edf8c_true_primary_v_recur_malignant_abundance.pdf"), width = 3.75, height = 3, useDingbats = FALSE, bg = "transparent")
ggplot(mal_care_md_state_summary %>% 
         # Recode and restrict to first resection
         mutate(clinical_time = recode(sample_time, `Init.` = "Primary",
                                       `Recur.` = "Recur.")) %>% 
         filter(timepoint!="T3", case_barcode%in%true_primaries$case_barcode, idh_codel_subtype=="Astro.", cell_state%in%c("AC-like", "MES-like", "Undifferentiated", "Cycling")), aes(x = clinical_time, y = freq*100)) + 
  geom_boxplot(outlier.shape = NA, aes(fill=cell_state)) +
  geom_line(aes(group=patient_id), color="gray85", linetype=2, size = 0.5) +
  geom_point(size = 0.5) +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  facet_grid(.~cell_state, scales="free") +
  labs(x="Matched longitudinal - astrocytoma", y="Malignant cell abundance (%)") +
  stat_n_text(size=2.25) +
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# INITIAL VS RECURRENCE stratified by acquired genetic alterations
mal_care_md_state_summary <- mal_care_md_state_summary %>% 
  mutate(acquired_genetic_alt = ifelse(is.na(acquired_genetic_alt_t1t2), "Did not acquire alt.",  "Acquired genetic alt."))

pdf(paste0(fig_dir, "fig3d_init_v_recur_genetics_malignant_abundance.pdf"), width = 4, height = 3.25, useDingbats = FALSE, bg = "transparent")
ggplot(mal_care_md_state_summary %>% 
         filter(timepoint!="T3", cell_state%in%c("Undifferentiated", "AC-like", "MES-like", "Cycling")), aes(x = sample_time, y = freq*100)) + 
  geom_boxplot(outlier.shape = NA, aes(fill=cell_state)) +
  geom_line(aes(group=patient_id), color="gray85", linetype=2, size = 0.5) +
  geom_point(size = 0.5) +
  plot_theme +
  theme(legend.position = "bottom") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  facet_grid(acquired_genetic_alt~cell_state, scales="free") +
  labs(x="Matched longitudinal analysis", y="Malignant cell abundance (%)") +
  stat_n_text(size = 2.25) +
  guides(fill=FALSE)
dev.off()


# SS2 analysis of malignant states
care_ss2_state_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/snrna/caremut_ss2_malignant_state_assignment.txt", sep="\t", header = TRUE)

care_ss2_state_md <- care_ss2_state_md %>%
  mutate(cell_state = recode(State, `MP_AC1_MUT` = "AC-like",
                             `MP_AC2_MUT` = "AC-like",
                             `MP_OPC_MUT` = "OPC-like",
                             `MP_NPC_MUT` = "NPC-like",
                             `MP_MES_MUT` = "MES-like",
                             "Undifferentiated" = "Undifferentiated"))

ss2_care_state_freq <- care_ss2_state_md %>% 
  group_by(care_id, cell_state) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts/sum(counts)) %>% 
  ungroup() %>% 
  complete(care_id, cell_state,
           fill = list(counts = 0, freq = 0)) %>%
  inner_join(md_trim, by="care_id") 

# Confirm that all sample sum 1 and that all cell types are measured.
ss2_care_state_freq %>%
  group_by(care_id) %>%
  summarise(freq_sum = sum(freq)) %>%
  mutate(pass = abs(freq_sum - 1) < 1e-10) %>%
  { if (all(.$pass)) {
    message("PASS: All samples frequencies sum to 1")
  } else {
    failing <- filter(., !pass)
    message("FAIL: ", nrow(failing), " samples do not sum to 1:")
    print(failing)
  }
  }
# Visual sanity check
ggplot(ss2_care_state_freq, aes(x = care_id, y=freq*100, fill = cell_state)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  plot_theme

ss2_care_state_freq_cc <- care_ss2_state_md %>% 
  group_by(care_id, isCC) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts/sum(counts)) %>% 
  ungroup() %>% 
  complete(care_id, isCC,
           fill = list(counts = 0, freq = 0)) %>%
  filter(isCC==TRUE) %>% 
  mutate(cell_state = "Cycling") %>% 
  inner_join(md_trim, by="care_id") %>% 
  dplyr::select(-isCC)

ss2_mal_care_md_state_summary <- ss2_care_state_freq %>% 
  bind_rows(ss2_care_state_freq_cc) %>% 
  mutate(platform = "SS2")
ss2_mal_care_md_state_summary$cell_state <- factor(ss2_mal_care_md_state_summary$cell_state, levels=c("AC-like", "MES-like", "Undifferentiated", "NPC-like", "OPC-like", "Cycling"))

pdf(paste0(fig_dir, "edf8d_ss2_longitudinal_changes_twosided_wilcox.pdf"), width = 3.75, height = 3, useDingbats = FALSE, bg = "transparent")
ggplot(ss2_mal_care_md_state_summary %>% 
         filter(cell_state%in%c("AC-like", "MES-like", "Undifferentiated", "Cycling")), aes(x = timepoint, y = freq*100)) + 
  geom_boxplot(aes(fill=cell_state), outlier.shape = NA) +
  geom_line(aes(group=patient_id), color="gray70", linetype=2, size = 0.5) +
  geom_point(size = 0.5) +
  plot_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  facet_grid(.~cell_state, scales="free") +
  labs(x="SmartSeq2 matched longitudinal analysis", y="Malignant cell abundance (%)") +
  stat_n_text(size = 2.25) +
  scale_x_discrete(labels = c("T1" = "Init.",
                              "T2" = "Recur.")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# 
mp_hierarchy_scores <- read.delim(file = paste0(proj_dir, "/results/scoring/malignant_malignant_hierarchy_signature_scores.txt"), sep = "\t", row.names = 1, header = TRUE)

mp_scores_hierarchy_median <- mp_hierarchy_scores %>% 
  inner_join(md_trim, by="SampleID") %>% 
  group_by(patient_id, timepoint) %>% 
  summarise(median_lineage = median(LineagePlot), median_stemness = median(Stemness)) %>% 
  ungroup() %>% 
  filter(timepoint!="T3") %>% 
  left_join(patient_md, by=c("patient_id")) %>% 
  mutate(acquired_genetic_alt_t1t2 = ifelse(is.na(acquired_genetic_alt_t1t2), "Did not acquire", "Acquired genetic alt."),
         timepoint = recode(timepoint, `T1` = "Init.",
                            `T2` = "Recur."))

pdf(paste0(fig_dir, "edf8e_idh_mutant_hierarchy_stemness_score_wilcox.pdf"), width = 3.25, height = 2.5, useDingbats = FALSE, bg = "transparent")
ggplot(mp_scores_hierarchy_median, aes(x=as.factor(timepoint), y=median_stemness)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = idh_codel_subtype), size = 0.75) + 
  geom_line(aes(group=case_barcode), color="gray70", linetype=2, size = 0.5) +
  stat_compare_means(method="wilcox", paired = TRUE, label="p.format", size = 2.25) +
  scale_color_manual(values=c("Astro." = "#800074",
                             "Oligo." = "#298C8C")) +
  facet_grid(.~acquired_genetic_alt_t1t2, scales="free_x", space = "free") +
  plot_theme +
  stat_n_text(size = 2.25) +
  labs(x="Matched longitudinal analysis", y="Median hierarchy stemness score", color = "Tumor type") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7))
dev.off()

# Cell state abundance heatmap
library(grid)
library(scales)

ac_like_freq_order <- mal_care_md_state_summary %>%
  filter(cell_state == "AC-like") %>%
  arrange(desc(freq)) %>%
  pull(SampleID)

mal_care_md_state_summary_heatmap <- mal_care_md_state_summary %>%
  mutate(SampleID = factor(SampleID, levels = ac_like_freq_order)) 
mal_care_md_state_summary_heatmap$cell_state <- factor(mal_care_md_state_summary_heatmap$cell_state, levels=rev(c("AC-like", "OPC-like", "Undifferentiated", "NPC-like", "MES-like", "Cycling")))

gg_heatmap <- ggplot(mal_care_md_state_summary_heatmap, aes(x = SampleID, y = cell_state, fill=freq*100)) +
  geom_tile() +
  scale_fill_viridis(limits = c(0.0, 50), oob=squish,
                     guide = guide_colorbar(
                       barwidth  = unit(0.8, "in"),
                       barheight = unit(0.08, "in"),
                       title.position = "top",
                       title.hjust = 0.5,
                       ticks.colour = "black",
                       frame.colour = "black"
                     )) +
  facet_grid(.~idh_codel_subtype, scales = "free", space = "free") +
  labs(x="", y="Malignant state", fill="State\nabundance (%)") + 
  plot_theme +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 8),
        strip.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))


gg_grade <- ggplot(mal_care_md_state_summary_heatmap %>% 
                     filter(cell_state=="AC-like"), 
                   aes(x = SampleID, y = 1, fill = factor(grade))) +
  geom_tile() +
  labs(y="") +
  plot_theme +
  scale_fill_manual(name = "Grade", 
                    values = c("G2" = "#fee0d2", "G3" = "#fc9272", "G4" = "#de2d26")) +
  facet_grid(. ~ idh_codel_subtype, scales = "free", space = "free") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    legend.key.size = unit(4, "pt"),    
    legend.spacing.x = unit(1, "pt"),   
    legend.spacing.y = unit(1, "pt")
  )


pdf(paste0(fig_dir, "fig3a_cell_state_abundance_heatmap.pdf"), width = 4, height = 3, useDingbats = FALSE)
egg::ggarrange(gg_heatmap, gg_grade, nrow = 2, heights = c(0.8, 0.1))
dev.off()



# WATERFALL plot of longitudinal change
wide_summary <- mal_care_md_state_summary %>% 
  filter(timepoint!="T3") %>% 
  dplyr::select(cell_state, freq, patient_id, patient_id, timepoint) %>% 
  pivot_wider(names_from = timepoint, values_from = freq) %>% 
  mutate(dT = T2-T1)

res <- wide_summary %>% 
  dplyr::select(patient_id, cell_state, dT) %>% 
  pivot_wider(names_from = cell_state, values_from = dT) %>% 
  inner_join(patient_md, by="patient_id") %>% 
  mutate(genetic_change = ifelse(is.na(acquired_genetic_alt_t1t2), "Did not acquire", "Acquired"))

res_ordered <- res[order(res$Undifferentiated, decreasing=TRUE),]
res_ordered <- res_ordered %>% 
  mutate(patient_id = as_factor(patient_id)) %>%
  mutate(patient_id = fct_relevel(patient_id, res_ordered$patient_id))
res_ordered$fill <- res_ordered$Undifferentiated > 0

gg_top <- ggplot(res_ordered, aes(x=patient_id, y = Undifferentiated*100, fill = factor(fill))) +
  geom_bar(stat="identity") +
  plot_theme +
  scale_fill_manual(values=c("dodgerblue","#FF0000")) +
  labs(y = "Longitudinal Undiff.\nabundance change (%)") +
  theme(plot.title = element_text(size=8, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y= element_text(size=8),
        text = element_text(size = 8), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        legend.position = "none")


res$patient_id <- factor(res$patient_id, levels= levels(res_ordered$patient_id))

gg_middle_ac <- ggplot(res, aes(x=patient_id, y=1, fill=`AC-like`*100)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "dodgerblue", high = "#FF0000", mid = "white", 
                       midpoint = 0, space = "Lab",
                       guide = guide_colorbar(direction = "horizontal", barwidth = 3, barheight = 1)) +
  labs(y = "", fill="AC-like\nchange (%)") + 
  plot_theme +
  theme(plot.title = element_text(size=8, hjust = 0.5),
        axis.ticks.x = element_blank(),
        text = element_text(size = 8), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        # no ticks / labels / titles for this strip
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks   = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.position = "right")

gg_middle_cc <- ggplot(res, aes(x=patient_id, y=1, fill=`Cycling`*100)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "dodgerblue", high = "#FF0000", mid = "white", 
                       midpoint = 0, space = "Lab",
                       guide = guide_colorbar(direction = "horizontal", barwidth = 3, barheight = 1)) +
  labs(y = "", fill="Cycling\nchange (%)") + 
  plot_theme +
  theme(plot.title = element_text(size=8, hjust = 0.5),
        axis.ticks.x = element_blank(),
        text = element_text(size = 8), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        # no ticks / labels / titles for this strip
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks   = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.position = "right")

gg_genetics <- ggplot(res, aes(x=patient_id, y = 1, fill = genetic_change)) +
  geom_tile() +
  labs(y="", fill="Recur.-assoc.\ngenetic event") +
  scale_fill_manual(values=c("black", "white"),na.value="#E5E5E5") +
  plot_theme +
  theme(plot.title = element_text(size=8, hjust = 0.5),
        axis.ticks.x = element_blank(),
        text = element_text(size = 8), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        # no ticks / labels / titles for this strip
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks   = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.position = "right")


p <- egg::ggarrange(gg_top, gg_middle_cc, gg_middle_ac, gg_genetics,  nrow = 4, heights = c(0.6, 0.2, 0.2, 0.2))


pdf(paste0(fig_dir, "fig3c_longitudinal_undifferentiated_change.pdf"), width = 4, height = 3, useDingbats = FALSE)
grid.newpage()
grid.draw(p)
dev.off()

# Extract the associated p-values and correlation coefficients.
ggplot(res, aes(x=factor(genetic_change), y=`Undifferentiated`)) +
  geom_boxplot() +
  stat_compare_means(method="wilcox") +
  plot_theme
wilcox.test(res$Undifferentiated~res$genetic_change)

ggplot(res, aes(x=`AC-like`, y=Undifferentiated)) +
  geom_point() +
  stat_cor(method="pearson") +
  plot_theme
cor.test(res$`AC-like`, res$Undifferentiated, method = "p")

ggplot(res, aes(x=Cycling, y=Undifferentiated)) +
  geom_point() +
  stat_cor(method="pearson") +
  plot_theme
cor.test(res$Cycling, res$Undifferentiated, method = "p")


### END ###