##################################
# Purpose: Analyze the myeloid cell abundance across relevant features such as grade and time point
# Author: Kevin Johnson
##################################

library(tidyverse) 
library(ggpubr) 
library(EnvStats) 

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/"
fig_dir     <- file.path(proj_dir, "figures/")
setwd(proj_dir)

source(paste0(proj_dir, "scripts/utils/plot_theme.R"))

# Load in the overall snRNA data annotation
md <- read.table(paste0(proj_dir, "processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20260320.txt"), sep = "\t", row.names = 1, header = TRUE)
md <- md %>% 
  mutate(patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2)) )

caremut_md_case <- md %>% 
  dplyr::select(SampleID, lab, idh_codel_subtype, care_id, patient_id, timepoint) %>% 
  distinct()

# Malignant cell state assignment
md_malignant <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/scoring/caremut_malignant_cell_state_assignment.txt", sep="\t", header = TRUE)
md_malignant <- md_malignant %>% 
  mutate(patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2)))

# Assigned metaprogram/state score based on p-value calculation
myeloid_md <- read.delim(paste0(proj_dir, "data/snrna/myeloid_cell_classification.tsv"), header = TRUE, sep = "\t")

# Miller et al Nature myeloid program for single cell analysis applied to this single nucleus RNA sequencing data
miller_myeloid_program <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/snrna/miller_myeloid_program_mode.txt", sep = "\t", row.names = 1, header = TRUE)
miller_myeloid_program$CellID <- gsub("\\.", "-", rownames(miller_myeloid_program))


# Read in clinical information
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt", sep="\t", header = TRUE)
patient_md$received_rt_t1t2[patient_md$patient_id=="P107"] <- "0"
sample_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240608.txt", sep="\t", header = TRUE)
sample_md_grade <- sample_md %>% 
  mutate(grade = paste0("G", grade_num)) %>% 
  dplyr::select(care_id, patient_id, timepoint, grade, idh_codel_subtype)

# Restrict to individual sample levels across any time point
myeloid_samples_min25 <- md %>% 
  group_by(SampleID, CellType_final) %>% 
  summarise(counts = n()) %>% 
  filter(CellType_final=="Myeloid") %>% 
  ungroup() %>% 
  filter(counts >= 25) 

# Inspect the relative myeloid cell state abundance across different grade and by tumor types
sample_grade_myeloid_freq <- myeloid_md %>% 
  filter(SampleID%in%myeloid_samples_min25$SampleID) %>% 
  group_by(care_id, myeloid_state_collapsed) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(care_id, myeloid_state_collapsed,
           fill = list(counts = 0, freq = 0)) %>% 
  inner_join(sample_md_grade, by="care_id") %>% 
  mutate(tumor_type = case_when(idh_codel_subtype == "IDHmut-noncodel" ~ "Astro.",
                                idh_codel_subtype == "IDHmut-codel" ~ "Oligo."))

table(sample_grade_myeloid_freq$myeloid_state_collapsed)

sample_grade_myeloid_freq$tumor_type <- factor(sample_grade_myeloid_freq$tumor_type, levels=c("Oligo.","Astro."))
sample_grade_myeloid_freq$myeloid_state_collapsed <-  factor(sample_grade_myeloid_freq$myeloid_state_collapsed, levels= c("Microglia", "Inflammatory", "Macrophage", "Unresolved"))


pdf(paste0(fig_dir, "edf10_myeloid_grade_abundance_kruskal.pdf"), width = 3, height = 2.75, useDingbats = FALSE, bg = "transparent")
ggplot(sample_grade_myeloid_freq %>% 
         filter(myeloid_state_collapsed!="Unresolved"), aes(x = grade, y = freq*100)) + 
  geom_boxplot(aes(fill=myeloid_state_collapsed), outlier.shape = NA) +
  scale_fill_manual(values=c("Macrophage" = "#08519c",
                             "Inflammatory" = "#9ecae1",
                             "Microglia" ="#eff3ff",
                             "Unresolved" ="gray90")) +
  geom_point(position = position_jitter(width = 0.1, seed = 42), 
             size = 0.5, alpha = 0.5) +
  plot_theme +
  theme(legend.position = "none") +
  # Kruskal for Astro.
  stat_compare_means(method = "kruskal", size = 2.25, label="p.format") +
  labs(x="Tumor grade (min. 25 myeloid cells)", y="Relative myeloid cell abundance (%)") +
  theme(strip.background = element_blank()) +
  facet_grid(tumor_type~myeloid_state_collapsed, scales="free") +
  stat_n_text(size = 2.25)
dev.off()

pdf(paste0(fig_dir, "edf10_myeloid_grade_abundance_wilcox.pdf"), width = 3, height = 2.75, useDingbats = FALSE, bg = "transparent")
ggplot(sample_grade_myeloid_freq %>% 
         filter(myeloid_state_collapsed!="Unresolved"), aes(x = grade, y = freq*100)) + 
  geom_boxplot(aes(fill=myeloid_state_collapsed), outlier.shape = NA) +
  scale_fill_manual(values=c("Macrophage" = "#08519c",
                             "Inflammatory" = "#9ecae1",
                             "Microglia" ="#eff3ff",
                             "Unresolved" ="gray90")) +
  geom_point(position = position_jitter(width = 0.1, seed = 42), 
             size = 0.5, alpha = 0.5) +
  plot_theme +
  theme(legend.position = "none") +
  # Wilcoxon for Oligo.
  stat_compare_means(method = "wilcox", size = 2.25, label="p.format") +
  labs(x="Tumor grade (min. 25 myeloid cells)", y="Relative myeloid cell abundance (%)") +
  theme(strip.background = element_blank()) +
  facet_grid(tumor_type~myeloid_state_collapsed, scales="free") +
  stat_n_text(size = 2.25)
dev.off()


# Restrict to patients with sufficient myeloid cells at both time points (T1 and T2) for longitudinal analyses (n = 32)
longitudinal_myeloid_subset <- md %>% 
  group_by(SampleID, patient_id, timepoint, CellType_final) %>% 
  summarise(counts = n()) %>% 
  filter(CellType_final=="Myeloid", timepoint!="T3") %>% 
  ungroup() %>% 
  filter(counts >= 25) %>% 
  group_by(patient_id) %>% 
  summarise(samples_sufficient_cells = n()) %>% 
  filter(samples_sufficient_cells > 1)

n_distinct(longitudinal_myeloid_subset$patient_id)

myeloid_long_freq <- myeloid_md %>% 
  filter(patient_id%in%longitudinal_myeloid_subset$patient_id, timepoint!="T3") %>% 
  group_by(SampleID, myeloid_state_collapsed) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  dplyr::select(SampleID, myeloid_state_collapsed, counts, freq) %>% 
  complete(SampleID, myeloid_state_collapsed,
           fill = list(counts = 0, freq = 0)) %>%
  distinct() 
# These would be the samples either without sufficient myeloid cells at T1/T2 or are T3 time points 
caremut_md_case$SampleID[!caremut_md_case$SampleID%in%unique(myeloid_long_freq$SampleID)]

# This should be 64 samples and 32 patients
table(myeloid_long_freq$myeloid_state_collapsed)
n_distinct(myeloid_long_freq$SampleID)


myeloid_summary_state <- myeloid_long_freq %>% 
  inner_join(caremut_md_case, by = "SampleID") %>% 
  left_join(patient_md, by=c("patient_id", "idh_codel_subtype")) %>% 
  mutate(sample_time = recode(timepoint, "T1" = "Init.",
                                     "T2" = "Recur.")) %>% 
  mutate(treatment_t1t2 = recode(treatment_t1t2, `1` = "treated", 
                                 `0` = "not treated")) %>% 
  mutate(received_rt_t1t2 = recode(received_rt_t1t2, `1` = "RT", 
                                   `0` = "No RT")) %>% 
  mutate(received_alk_t1t2 = recode(received_alk_t1t2, `1` = "Alk.", 
                                    `0` = "No Alk.")) %>% 
  mutate(subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                          `IDHmut-noncodel` = "Astro."))

myeloid_summary_state$subtype <- factor(myeloid_summary_state$subtype, levels=c("Oligo.","Astro."))
myeloid_summary_state$myeloid_state_collapsed <-  factor(myeloid_summary_state$myeloid_state_collapsed, levels= c("Microglia", "Inflammatory", "Macrophage", "Unresolved"))

# Longitudinal change in relative myeloid cell abundance
pdf(paste0(fig_dir, "edf10_myeloid_longitudinal_abundance.pdf"), width = 2.25, height = 2.5, useDingbats = FALSE, bg = "transparent")
ggplot(myeloid_summary_state %>% 
         filter(myeloid_state_collapsed%in%c("Microglia", "Macrophage")), aes(x = sample_time, y = freq*100)) + 
  geom_boxplot(aes(fill=myeloid_state_collapsed)) +
  geom_line(aes(group=patient_id), color="gray70", linetype=2, size  = 0.5) +
  geom_point(size = 0.5) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  scale_fill_manual(values=c("Macrophage" = "#08519c",
                             "Inflammatory" = "#9ecae1",
                             "Microglia" ="#eff3ff",
                             "Unresolved" ="gray90")) +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  labs(x="Longitudinal pairs (min. 25 myeloid cells)", y="Relative myeloid cell abundance (%)") +
  theme(strip.background = element_blank()) +
  facet_grid(myeloid_state_collapsed~., scales="free") +
  stat_n_text(size = 2.25)
dev.off()

### ### ### ### ### ###
# Linear models testing association between relative myeloid abundance change and clinical variables
### ### ### ### ### ###
patient_md_trim <- patient_md %>% 
  dplyr::select(patient_id, idh_codel_subtype, grade_change_t1t2, surgical_interval_mo_t1t2, received_alk_t1t2, received_rt_t1t2) 
patient_md_trim$received_rt_t1t2[patient_md_trim$patient_id=="P17"] <- 0
patient_md_trim$received_rt_t1t2[patient_md_trim$patient_id=="P107"] <- 0

myeloid_care_md_summary_wide <- myeloid_summary_state %>% 
  dplyr::select(patient_id, timepoint, myeloid_state_collapsed, freq) %>% 
  pivot_wider(names_from = timepoint, values_from = freq) %>% 
  mutate(dT = T2-T1) %>% 
  left_join(patient_md_trim, by="patient_id") %>% 
  mutate(received_rt_t1t2 = recode(received_rt_t1t2, `1` = "RT", 
                                   `0` = "No RT")) %>% 
  mutate(received_alk_t1t2 = recode(received_alk_t1t2, `1` = "Alk.", 
                                    `0` = "No Alk.")) %>% 
  mutate(tumor_type = recode(idh_codel_subtype, `IDHmut-noncodel` = "Astro.", 
                                    `IDHmut-codel` = "Oligo.")) 
myeloid_care_md_summary_wide$tumor_type <- factor(myeloid_care_md_summary_wide$tumor_type, levels=c("Oligo.","Astro."))
myeloid_care_md_summary_wide$received_alk_t1t2 <- factor(myeloid_care_md_summary_wide$received_alk_t1t2, levels=c("No Alk.","Alk."))

# Microglia
myeloid_care_md_summary_wide_mg <- myeloid_care_md_summary_wide %>% 
  filter(myeloid_state_collapsed=="Microglia")

# Delta Microglia abundance is approximately normally distributed
hist(myeloid_care_md_summary_wide_mg$dT)

# MG-like shift at recurrence
fit <- lm(dT ~ grade_change_t1t2 + received_rt_t1t2 + surgical_interval_mo_t1t2, data = myeloid_care_md_summary_wide_mg)
summary(fit) # RT p-value = 0.004
fit <- lm(dT ~ grade_change_t1t2 + received_rt_t1t2 + surgical_interval_mo_t1t2 + tumor_type, data = myeloid_care_md_summary_wide_mg)
summary(fit) # RT p-value = 0.005
fit <- lm(dT ~ grade_change_t1t2 + received_rt_t1t2 + received_alk_t1t2 + surgical_interval_mo_t1t2 + tumor_type, data = myeloid_care_md_summary_wide_mg)
summary(fit) # RT p-value = 0.01

# Macrophage
myeloid_care_md_summary_wide_mac <- myeloid_care_md_summary_wide %>% 
  filter(myeloid_state_collapsed=="Macrophage")
# Delta Macrophage abundance is approximately normally distributed
hist(myeloid_care_md_summary_wide_mac$dT)

# Wilcoxon p-value = 0.009 for samples treated with RT versus those that were not.
wilcox.test(myeloid_care_md_summary_wide_mac$dT~as.factor(myeloid_care_md_summary_wide_mac$received_rt_t1t2))

# Macrophage-like shift at recurrence
fit <- lm(dT ~ grade_change_t1t2 + received_rt_t1t2 + surgical_interval_mo_t1t2, data = myeloid_care_md_summary_wide_mac)
summary(fit) # RT p-value = 0.003
fit <- lm(dT ~ grade_change_t1t2 + received_rt_t1t2 + surgical_interval_mo_t1t2 + tumor_type, data = myeloid_care_md_summary_wide_mac)
summary(fit) # RT p-value = 0.01
fit <- lm(dT ~ grade_change_t1t2 + received_rt_t1t2 + received_alk_t1t2 + surgical_interval_mo_t1t2 + tumor_type, data = myeloid_care_md_summary_wide_mac)
summary(fit) # RT p-value = 0.051. RT p-value = 0.04 when alk. NAs were classified as no alk. treatment

# What's the strength of the relationship between Microglia and Macrophage changes?
myeloid_care_md_summary_wide_mac$patient_id==myeloid_care_md_summary_wide_mg$patient_id
cor.test(myeloid_care_md_summary_wide_mac$dT, myeloid_care_md_summary_wide_mg$dT, method="s")

### ### ### ### ### ### ###
# Miller et al program annotation
### ### ### ### ### ### ###
care_myeloid_activity_avg <- miller_myeloid_program %>% 
  inner_join(myeloid_md, by="CellID") %>% 
  dplyr::select(CellID, care_id, Microglia:Monocyte) %>% 
  pivot_longer(cols = c(Microglia:Monocyte), names_to = "program", values_to = "activity") %>% 
  group_by(care_id, program) %>% 
  summarise(avg_activity = mean(activity),
            cell_counts = n()) %>% 
  ungroup() %>% 
  mutate(patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2))) %>% 
  filter(patient_id%in%longitudinal_myeloid_subset$patient_id, timepoint!="T3") %>% 
  left_join(patient_md_trim, by="patient_id") %>% 
  filter(program%in%c("Complement_Immunosuppressive", "Inflammatory_microglia", "Microglia", "Scavenger_Immunosuppressive")) %>% 
  arrange(program, patient_id, timepoint) 

care_myeloid_activity_avg$program <- factor(care_myeloid_activity_avg$program, levels=c("Microglia", "Inflammatory_microglia", "Complement_Immunosuppressive", "Scavenger_Immunosuppressive"))

pdf(paste0(fig_dir, "miller_myeloid_program_usage_rt_longitudinal.pdf"), width = 3,5, height = 2.75, useDingbats = FALSE, bg = "transparent")
ggplot(care_myeloid_activity_avg %>% 
         filter(timepoint!="T3", !is.na(received_rt_t1t2)) %>% 
         mutate(received_rt_t1t2 = recode(received_rt_t1t2, `1` = "RT",
                                          `0` = "No RT"),
                timepoint = recode(timepoint, `T1` = "I",
                                   `T2` = "R")), aes(x=timepoint, y=avg_activity)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group=patient_id), color="gray70", linetype=2, size  = 0.5) +
  geom_point(size = 0.5) +
  facet_grid(received_rt_t1t2~program, scales = "free", space="free") +
  plot_theme +
  stat_compare_means(method="wilcox", paired = TRUE, label = "p.format", size = 2.25) +
  stat_n_text(size = 2.25) +
  labs(x = "Longitudinal pairs (min. 25 myeloid cells)", y = "Mean % program usage per tumor") +
  theme(strip.text = element_text(size = 8)) 
dev.off()


### ### ### ### ### ###
# Correlation heat map
### ### ### ### ### ###
library(magrittr) 
library(corrplot)
library(Hmisc)

myeloid_md_trim <- myeloid_md %>% 
  dplyr::select(CellID, myeloid_state_collapsed)

# Difference in objects will be removal of Unresolved
md_all_plus_malignant <- md %>% 
  filter(CellType_final!="Unresolved") %>% 
  left_join(md_malignant, by=c("CellID", "SampleID", "case_barcode", "idh_codel_subtype", "care_id", "patient_id", "timepoint")) %>% 
  mutate(MalignantState = recode(State, `MP_AC1_MUT` = "AC-like",
                                 `MP_OPC_MUT` = "OPC-like",
                                 `MP_NPC_MUT` = "NPC-like",
                                 `MP_MES_MUT` = "MES-like",
                                 `MP_AC2_MUT` = "AC-like",
                                 "Undifferentiated" = "Undifferentiated")) %>% 
  left_join(myeloid_md_trim, by=c("CellID")) %>% 
  mutate(myeloid_state = recode(myeloid_state_collapsed, `Unresolved` = "TAM Unresolved",
                                `Inflammatory` = "TAM Inflammatory",
                                `Macrophage` = "TAM Macrophage",
                                `Microglia` = "TAM Microglia"))

# Redfine CellType variable so that it accounts for Malignant and Myeloid sub-compartments
md_all_plus_malignant$CellTypeRefined <- md_all_plus_malignant$CellType_final
md_all_plus_malignant$CellTypeRefined[md_all_plus_malignant$CellTypeRefined=="Myeloid"] <- md_all_plus_malignant$myeloid_state[md_all_plus_malignant$CellTypeRefined=="Myeloid"]
md_all_plus_malignant$CellTypeRefined[md_all_plus_malignant$CellTypeRefined=="Malignant"] <- md_all_plus_malignant$MalignantState[md_all_plus_malignant$CellTypeRefined=="Malignant"]

# 7 malignant cells from a sample that was removed from downstream malignant cell state assignment analyses
md_all_plus_malignant %>% 
  filter(is.na(CellTypeRefined)) %>% 
  dplyr::select(care_id)

# Summarize the frequency of each CellType across samples
care_all_celltypes_summary <- md_all_plus_malignant %>% 
  # Remove sample with very few malignant cells
  filter(care_id!="P99T3") %>% 
  group_by(SampleID, CellTypeRefined) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  dplyr::select(SampleID, CellTypeRefined, counts, freq) %>% 
  complete(SampleID, CellTypeRefined,
           fill = list(counts = 0, freq = 0)) %>%
  distinct() %>% 
  inner_join(caremut_md_case, by="SampleID") %>% 
  inner_join(patient_md, by=c("patient_id", "idh_codel_subtype")) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."))


state_freq <- care_all_celltypes_summary %>% 
  # Remove sample with no TME and time point 3 samples
  filter(timepoint!="T3", care_id!="P91T2") %>% 
  dplyr::select(care_id, CellTypeRefined, freq) %>% 
  pivot_wider(names_from = CellTypeRefined, values_from = freq)
# 69 total samples in this analysis.
n_distinct(state_freq$care_id)

# Performing adjustment for multiple hypothesis testing
state_freq_cor <- cor(as.matrix(state_freq[ ,2:17]), method = "pearson")
corRes = cor.mtest(as.matrix(state_freq_cor), conf.level = 0.95)
corRespAdj <- p.adjust(c(corRes[[1]]), method = "fdr")
resAdj <- matrix(corRespAdj, ncol = dim(corRes[[1]])[1])
colnames(resAdj) <- colnames(corRes$p)
rownames(resAdj) <- rownames(corRes$p)


pdf(paste0(fig_dir, "edf10_tme_malignant_state_correlation_pvalue_adjusted_heatmap.pdf"), width = 3.25, height = 2.75, useDingbats = FALSE, bg = "transparent")
corrplot(state_freq_cor, p.mat = resAdj, method = 'circle', type = 'lower',
         tl.col = "black",
         order = 'FPC', # 'FPC' for the first principal component order.
         insig='blank', number.cex = 0.5, cl.cex = 0.5, tl.cex = 0.5, diag=FALSE, col=rev(COL2("RdBu"))) 
dev.off()

### ### ### ### ### ### ###
# Delta malignant and delta TME cell states
### ### ### ### ### ### ###
# Separate TME and malignant to adjust for potential confounding from purity
care_md_summary_mal <- md_all_plus_malignant %>% 
  # Only consider T1 and T2 samples
  filter(timepoint!="T3") %>% 
  # Restrict to malignant cells
  filter(CellTypeRefined%in%c("MES-like", "NPC-like", "AC-like", "OPC-like", "Undifferentiated")) %>% 
  group_by(SampleID, CellTypeRefined) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  dplyr::select(SampleID, CellTypeRefined, counts, freq) %>% 
  complete(SampleID, CellTypeRefined,
           fill = list(counts = 0, freq = 0)) %>%
  distinct() %>% 
  inner_join(caremut_md_case, by="SampleID") %>% 
  inner_join(patient_md, by=c("patient_id", "idh_codel_subtype")) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."))

# Should be 70 samples and 35 patients
table(care_md_summary_mal$CellTypeRefined)

care_md_summary_tme <- md_all_plus_malignant %>% 
  # Only consider T1 and T2 samples
  filter(timepoint!="T3") %>% 
  # Restrict to Non-malignant cells
  filter(!CellTypeRefined%in%c("MES-like", "NPC-like", "AC-like", "OPC-like", "Undifferentiated")) %>% 
  group_by(SampleID, CellTypeRefined) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  dplyr::select(SampleID, CellTypeRefined, counts, freq) %>% 
  complete(SampleID, CellTypeRefined,
           fill = list(counts = 0, freq = 0)) %>%
  distinct() %>% 
  inner_join(caremut_md_case, by="SampleID") %>% 
  inner_join(patient_md, by=c("patient_id", "idh_codel_subtype")) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."))
# This results in 69 samples because one sample does not have any TME cells
table(care_md_summary_tme$CellTypeRefined)
# P91T2 is the same with missing non-malignant cells. We'll remove this patient from analysis since the lack of a microenvironment will yield unusual delta (T2-T1) values
unique(care_md_summary_mal$care_id)[!unique(care_md_summary_mal$care_id)%in%unique(care_md_summary_tme$care_id)]

# Should be 68 rows - one for sample (34 patients)
state_freq_adj <- care_md_summary_mal %>% 
  bind_rows(care_md_summary_tme) %>% 
  filter(patient_id !="P91") %>% 
  dplyr::select(care_id, CellTypeRefined, freq) %>% 
  pivot_wider(names_from = CellTypeRefined, values_from = freq)

# Create variables to be able to pivot. These represent the delta (T2-T1) values per cell type
state_freq_adj_t1t2 <- state_freq_adj %>% 
  mutate(patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2))) %>% 
  dplyr::select(-care_id) %>% 
  pivot_longer(cols= c(`AC-like`:`TAM Unresolved`),
               names_to = "cell_type",
               values_to = "freq") %>% 
  pivot_wider(names_from = timepoint, values_from = freq) %>% 
  mutate(delta_freq = T2-T1) %>% 
  dplyr::select(-c(T1, T2)) %>% 
  pivot_wider(names_from = cell_type, values_from = delta_freq) 

# Pearson correlation coefficient 0.57 and P = 0.0004
ggplot(state_freq_adj_t1t2, aes(x = `TAM Macrophage`, y= `MES-like`)) +
  geom_point() +
  stat_cor(method="pearson")

# Compute the Pearson correlation across all cell states
state_freq_adj_t1t2_cor <- cor(as.matrix(state_freq_adj_t1t2[ ,2:17]))

# We want to plot the delta malignant cell abundance against the microenvironmental cell abundance so restrict to those features.
state_freq_adj_t1t2_cor_input <- as.data.frame(state_freq_adj_t1t2_cor[1:5, 6:16])
# All rows are the malignant states
state_freq_adj_t1t2_cor_input$malignant_state <- rownames(state_freq_adj_t1t2_cor_input)

# Converting to long format so that we can merge downstream.
state_freq_adj_t1t2_cor_input_long <- state_freq_adj_t1t2_cor_input %>% 
  pivot_longer(cols = c(Astrocyte:`TAM Unresolved`), names_to = "tme_state", values_to = "pearson_cor") %>% 
  mutate(comparison = paste0(malignant_state, "_", tme_state))

# Calculate the p-values for each correlation
state_freq_adj_cor_t1t2_pval <- rcorr(as.matrix(state_freq_adj_t1t2[ ,2:17]))
state_freq_adj_cor_t1t2_pval_df <-state_freq_adj_cor_t1t2_pval$P

# Restrict to malignant vs TME
state_freq_adj_t1t2_pvalue_input <- as.data.frame(state_freq_adj_cor_t1t2_pval_df[1:5, 6:16])
state_freq_adj_t1t2_pvalue_input$malignant_state <- rownames(state_freq_adj_t1t2_pvalue_input)

state_freq_adj_t1t2_pvalue_input_long <- state_freq_adj_t1t2_pvalue_input %>% 
  pivot_longer(cols = c(Astrocyte:`TAM Unresolved`), names_to = "tme_state", values_to = "pearson_pval") %>% 
  mutate(comparison = paste0(malignant_state, "_", tme_state)) %>% 
  mutate(adj_pval = p.adjust(pearson_pval,  method = "fdr")) %>% 
  inner_join(state_freq_adj_t1t2_cor_input_long, by=c("comparison", "malignant_state", "tme_state"))

mal_order <- c("AC-like",
               "MES-like",
               "Undifferentiated",
               "OPC-like",
               "NPC-like")

tme_order <- c("Astrocyte", 
               "ExcNeuron",
               "InhNeuron",
               "Oligodendrocyte",
               "Endothelial",
               "Mural",
               "Lymphocyte",
               "TAM Macrophage",
               "TAM Inflammatory",
               "TAM Microglia",
               "TAM Unresolved")

state_freq_adj_t1t2_pvalue_input_long$malignant_state <- factor(state_freq_adj_t1t2_pvalue_input_long$malignant_state, levels=mal_order)
state_freq_adj_t1t2_pvalue_input_long$tme_state <- factor(state_freq_adj_t1t2_pvalue_input_long$tme_state, levels=tme_order)

pdf(paste0(fig_dir, "tme_malignant_delta_state_corr_heatmap.pdf"), height = 3, width = 3.25, useDingbats = FALSE, bg = "transparent")
ggplot(data = state_freq_adj_t1t2_pvalue_input_long, aes(x=tme_state, y=malignant_state, fill = pearson_cor)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "dodgerblue", high = "#FF0000", mid = "white", 
                       midpoint = 0, limit = c(-0.8,0.8), oob = scales::squish, space = "Lab",
                       name="Pearson's corr. coeff.") +
  guides(fill = guide_colourbar(barheight = 1, barwidth = 3.25)) +
  labs(x= "Tumor microenvironment cell type", y = "Malignant state") +
  plot_theme +
  theme(legend.position = "top",
        axis.text.x = element_text(angle=45, hjust=1),
        plot.title = element_text(size = 8)) +
  geom_text(data = filter(state_freq_adj_t1t2_pvalue_input_long, adj_pval<0.05), aes(x = 8, y = 2), label = "*")
dev.off()


### ### ### ### ### ### ### ###
# Organoid myeloid cell program usage
### ### ### ### ### ### ### ###
organoid_myeloid_program <- read.table(paste0(proj_dir, "data/snrna/organoid_myeloid_miller_cell_classification.tsv"), sep = "\t",  header = TRUE)
organoid_myeloid_program$grade <- gsub("Grade ", "G", organoid_myeloid_program$Grade) 

organoid_myeloid_program$organoid_id <- paste0(organoid_myeloid_program$GBO_ID, " (", organoid_myeloid_program$grade, ")")
organoid_myeloid_program$treatment <- sapply(strsplit(organoid_myeloid_program$SampleID, "_"), "[[", 2)

organoid_myeloid_activity <- organoid_myeloid_program %>% 
  dplyr::select(CellID, SampleID, organoid_id, treatment, Microglia:Monocyte) %>% 
  pivot_longer(cols = c(Microglia:Monocyte), names_to = "program", values_to = "activity") %>% 
  mutate(treatment = recode(treatment, `ctl` = "control",
                            `10Gy` = "irradiation")) %>% 
  filter(program%in%c("Microglia", "Macrophage"))

organoid_myeloid_activity$program <- factor(organoid_myeloid_activity$program, levels = c("Microglia", "Macrophage"))

pdf(paste0(fig_dir, "edf10_organoid_myeloid_miller.pdf"), width = 2.75, height = 2.75, useDingbats = FALSE, bg = "transparent")
ggplot(organoid_myeloid_activity, aes(x=treatment, y=activity)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1, seed = 42), 
             size = 0.5, alpha = 0.5) +
  facet_grid(organoid_id~program, scales = "free", space="free") +

  labs(x = "Treatment", y = "% program usage per cell (Miller et al)") +
  plot_theme +
  theme(strip.text = element_text(size = 8)) +
  stat_compare_means(method="wilcox", label="p.format", size = 2.25) +
  stat_n_text(size = 2.25) 
dev.off()


### ### ### ### ### ### ###
# Waterfall plot
### ### ### ### ### ### ###
library(grid)

res <- myeloid_care_md_summary_wide %>% 
  dplyr::select(patient_id, myeloid_state_collapsed, dT) %>% 
  pivot_wider(names_from = myeloid_state_collapsed, values_from = dT) %>% 
  inner_join(patient_md_trim, by = "patient_id") %>% 
  mutate(received_rt_t1t2 = recode(received_rt_t1t2, `1` = "RT", 
                                   `0` = "No RT")) 

ggplot(res, aes(x=factor(received_rt_t1t2), y=Macrophage)) +
  geom_boxplot() +
  stat_compare_means(method="wilcox") +
  plot_theme +
  stat_n_text()

ggplot(res, aes(x = Microglia, y = Macrophage)) +
  geom_point() +
  stat_cor(method="pearson") +
  plot_theme 

res_ordered <- res[order(res$Macrophage, decreasing=TRUE),]
res$patient_id <- factor(res$patient_id, levels= res_ordered$patient_id)
res$fill <- res$Macrophage > 0

gg_top <- ggplot(res, aes(x=patient_id, y = Macrophage*100, fill = factor(fill))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("dodgerblue","#FF0000")) +
  guides(fill = guide_colourbar(barheight = 1, barwidth = 3.25)) +
  labs(y = "Longitudinal macrophage state\nabundance change (%)") +
  plot_theme +
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


gg_middle_micro <- ggplot(res, aes(x=patient_id, y=1, fill=Microglia*100)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "dodgerblue", high = "#FF0000", mid = "white", 
                       midpoint = 0, space = "Lab",
                       guide = guide_colorbar(direction = "horizontal", barwidth = 2., barheight = 0.75)) +
  labs(y = "", fill="Microglia\nchange (%)") + 
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
        axis.ticks.length = unit(0, "pt"))


gg_bottom <- ggplot(res, aes(x=patient_id, y = 1, fill = received_rt_t1t2)) +
  geom_tile() +
  labs(y="", fill="Radiotherapy") +
  scale_fill_manual(values=c("white","black"),na.value="#E5E5E5") +
  plot_theme +
  theme(plot.title = element_text(size=8, hjust = 0.5),
        axis.ticks.x = element_blank(),
        text = element_text(size = 8), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        # no ticks / labels / titles for this strip
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks   = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        legend.position = "right")

p <- egg::ggarrange(gg_top, gg_middle_micro, gg_bottom, nrow = 3, heights = c(0.6, 0.2, 0.2))


pdf(paste0(fig_dir, "myeloid_longitudinal_abundance_waterfall.pdf"), width = 3.75, height = 2.75, useDingbats = FALSE)
grid.newpage()
grid.draw(p)
dev.off()
