##################################
# Examine distributions of the CARE IDH-mutant cell states across genetic alterations while controlling for grade
# Author: Kevin Johnson
# Date Updated: 2026.04.03
##################################

# Reproduce the following figures: 4a

library(tidyverse)
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

mal_care_md_state_summary_long <- care_state_freq %>% 
  bind_rows(care_state_freq_cc) %>% 
  inner_join(sample_md, by=c("care_id")) %>% 
  mutate(grade = paste0("G", grade_num),
         # Indicating whether the sample was from the first surgery or not
         primary_v_recurrence = ifelse(surgery_number==1, "Primary", "Recurrence"))


### Gene-level CNAs
# Sanity check for Grade 4, astrocytomas
sample_g4_md <- sample_md %>% 
  mutate(grade = paste0("G", grade_num)) %>% 
  dplyr::select(care_id, grade, cdkn2a_del, pdgfra_amp) %>% 
  distinct() %>% 
  filter(grade=="G4")

# Not all samples had high-quality copy number calls. Deletions being harder to confidently detect in WXS data than amplifications.
# 19/24 G4 astrocytomas with CDKN2A status
sum(is.na(sample_g4_md$cdkn2a_del))
table(sample_g4_md$cdkn2a_del)

# 23/24 G4 astrocytomas with PDGFRA amp status
sum(is.na(sample_g4_md$pdgfra_amp))
table(sample_g4_md$pdgfra_amp)

# SN05-1 (G4) doesn't have any accompanying genetic data
mal_care_md_state_summary_long %>% 
  filter(is.na(pdgfra_amp), grade=="G4")

# PDGFRA amplification - G4, Astro. 
ggplot(mal_care_md_state_summary_long %>% 
         filter(!is.na(pdgfra_amp), idh_codel_subtype.x=="Astro.", grade=="G4"), aes(x = as.factor(pdgfra_amp), y = freq*100)) + 
  geom_boxplot(aes(fill=cell_state), outlier.shape = NA) +
  geom_point() +
  plot_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  stat_compare_means(method = "wilcox", paired = FALSE, size =2, label="p.format") +
  facet_grid(grade~cell_state, scales="free") +
  labs(x="Grade 4 Astro.", y="Malignant cell abundance (%)") +
  stat_n_text()

tmp1 <- mal_care_md_state_summary_long %>% 
  filter(!is.na(pdgfra_amp), idh_codel_subtype.x=="Astro.", grade=="G4") 
n_distinct(tmp1$care_id)

# CDKN2A homozygous deletion - G4, Astro. 
ggplot(mal_care_md_state_summary_long %>% 
         mutate(cdkn2a_homdel = ifelse(cdkn2a_del==-2, "homdel", "WT")) %>% 
         filter(!is.na(cdkn2a_del), idh_codel_subtype.x=="Astro.", grade=="G4"), aes(x = as.factor(cdkn2a_homdel), y = freq*100)) + 
  geom_boxplot(aes(fill=cell_state), outlier.shape = NA) +
  geom_point() +
  plot_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  stat_compare_means(method = "wilcox", paired = FALSE, size =2, label="p.format") +
  facet_grid(grade~cell_state, scales="free") +
  labs(x="Grade 4 Astro.", y="Malignant cell abundance (%)") +
  stat_n_text()


# Treatment-associated hypermutation: G3, Oligo. (needs to be a recurrence for treatment-associated hypermutation) 
hypermutation_md <- sample_md %>% 
  filter(!is.na(hypermutation), grade_num=="3", idh_codel_subtype=="Oligo.", surgery_number!=1)
table(hypermutation_md$hypermutation)

hypermutant_to_plot <- mal_care_md_state_summary_long %>% 
  filter(primary_v_recurrence=="Recurrence",  !is.na(hypermutation), grade=="G3", idh_codel_subtype.x=="Oligo.") %>% 
  dplyr::select(care_id, grade, primary_v_recurrence, hypermutation) %>% 
  distinct()

# The one sample that is missing from the cell state proportions is the sample with insufficient malignant cells.
hypermutation_md$care_id[!hypermutation_md$care_id%in%hypermutant_to_plot$care_id]

ggplot(mal_care_md_state_summary_long %>% 
         filter(primary_v_recurrence=="Recurrence",  !is.na(hypermutation), grade=="G3", idh_codel_subtype.x=="Oligo."), aes(x = as.factor(hypermutation), y = freq)) + 
  geom_boxplot(aes(fill=cell_state), outlier.shape = NA) +
  geom_point() +
  plot_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6", 
                             "AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  stat_compare_means(method = "wilcox", paired = FALSE, size =2.5, label="p.format") +
  facet_grid(.~cell_state, scales="free") +
  stat_n_text() +
  labs(x="", y="Relative malignant cell\nproportion (%)")


# Reformat so that all data is presented in one image restricted to AC-like, Undifferentiated, and Cycling cells
olig_hm <- mal_care_md_state_summary_long %>% 
  filter(primary_v_recurrence=="Recurrence",  !is.na(hypermutation), grade=="G3", idh_codel_subtype.x=="Oligo.",
         cell_state%in%c("AC-like", "Undifferentiated", "Cycling")) %>% 
  mutate(genetic_alt = recode(hypermutation, `1` = "Alk-assoc. HM",
                              `0` = "Non-HM"),
         tumor_type = recode(idh_codel_subtype.x, `Oligo.` = "Oligo., G3",
                             `Astro.` = "Astro., G4")) %>% 
  dplyr::select(care_id, cell_state, freq, tumor_type, genetic_alt)

astro_g4_cdkn2a <- mal_care_md_state_summary_long %>% 
  mutate(genetic_alt = ifelse(cdkn2a_del==-2, "CDKN2A homdel", "CDKN2A WT")) %>% 
  filter(!is.na(genetic_alt), idh_codel_subtype.x=="Astro.", grade=="G4",
         cell_state%in%c("AC-like", "Undifferentiated", "Cycling")) %>% 
  mutate(tumor_type = recode(idh_codel_subtype.x, `Oligo.` = "Oligo., G3",
                             `Astro.` = "Astro., G4")) %>% 
  dplyr::select(care_id, cell_state, freq, tumor_type, genetic_alt)

astro_g4_pdgfra <- mal_care_md_state_summary_long %>% 
  filter(!is.na(pdgfra_amp), idh_codel_subtype.x=="Astro.", grade=="G4",
         cell_state%in%c("AC-like", "Undifferentiated", "Cycling")) %>% 
  mutate(genetic_alt = recode(pdgfra_amp, `1` = "PDGFRA amp",
                              `0` = "PDGFRA WT"),
         tumor_type = recode(idh_codel_subtype.x, `Oligo.` = "Oligo., G3",
                             `Astro.` = "Astro., G4")) %>% 
  dplyr::select(care_id, cell_state, freq, tumor_type, genetic_alt)

all_genetic_alt <- bind_rows(olig_hm, astro_g4_cdkn2a, astro_g4_pdgfra)
genetic_alt_order <- c("CDKN2A WT",
                       "CDKN2A homdel",
                       "PDGFRA WT",
                       "PDGFRA amp",
                       "Non-HM",
                       "Alk-assoc. HM")
jitter <- position_jitter(width = .15)
all_genetic_alt$genetic_alt <- factor(all_genetic_alt$genetic_alt, levels = genetic_alt_order)
all_genetic_alt$tumor_type <- factor(all_genetic_alt$tumor_type, levels = c("Oligo., G3", "Astro., G4"))
all_genetic_alt$cell_state <- factor(all_genetic_alt$cell_state, levels = c("AC-like", "Undifferentiated", "Cycling"))


pval_annotations <- data.frame(
  cell_state  = c("AC-like", "AC-like", "AC-like",
                  "Undifferentiated", "Undifferentiated", "Undifferentiated",
                  "Cycling", "Cycling", "Cycling"),
  tumor_type  = c("Astro., G4", "Astro., G4", "Oligo., G3",
                  "Astro., G4", "Astro., G4", "Oligo., G3",
                  "Astro., G4", "Astro., G4", "Oligo., G3"),
  genetic_alt = c("CDKN2A homdel", "PDGFRA amp", "Alk-assoc. HM",
                  "CDKN2A homdel", "PDGFRA amp", "Alk-assoc. HM",
                  "CDKN2A homdel", "PDGFRA amp", "Alk-assoc. HM"),
  y           = c(70, 70, 70, 90, 90, 90, 40, 40, 40),  
  label       = c("p=0.18", "p=0.07", "p=0.048",          
                  "p=0.09", "p=0.001", "p=1.6e-04",
                  "p=0.003", "p=3.3e-05", "p=0.15")
) 

pval_annotations$genetic_alt <- factor(pval_annotations$genetic_alt, levels = genetic_alt_order)
pval_annotations$tumor_type <- factor(pval_annotations$tumor_type, levels = c("Oligo., G3", "Astro., G4"))
pval_annotations$cell_state <- factor(pval_annotations$cell_state, levels = c("AC-like", "Undifferentiated", "Cycling"))

pdf(paste0(fig_dir, "fig4a_genetic_alteration_malignant_abundance.pdf"), width = 4.5, height = 4, useDingbats = FALSE, bg = "transparent")
ggplot(all_genetic_alt, aes(x = genetic_alt, y = freq*100, fill=genetic_alt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = jitter, alpha = 0.25, size=0.75) +
  labs(x = "",
       y = "Malignant cell abundance (%)",
       fill = "Gene alteration") +
  facet_grid(tumor_type~cell_state, scales="free", space="free") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=c("Non-HM" = "white", 
                             "Alk-assoc. HM" ="#33A02C",
                             "PDGFRA WT" ="white",
                             "PDGFRA amp" = "#E31A1C", 
                             "CDKN2A WT" ="white",
                             "CDKN2A homdel" = "#1F78B4")) +
  theme(strip.text.x = element_text(margin = margin(t = 15, unit = "pt"))) +
  guides(fill=FALSE) +
  stat_n_text(size=2.25) +
  geom_text(data = pval_annotations, 
              aes(x = genetic_alt, y = y, label = label),
              inherit.aes = FALSE,
              size = 2.5,
            nudge_x = -0.3)
dev.off()

pdf(paste0(fig_dir, "fig4a_genetic_alteration_malignant_abundance_legend.pdf"), width = 5, height = 5, useDingbats = FALSE, bg = "transparent")
ggplot(all_genetic_alt, aes(x = genetic_alt, y = freq*100, fill=genetic_alt)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = jitter, alpha = 0.25, size=0.75) +
  labs(x = "",
       y = "Malignant cell abundance (%)",
       fill = "Gene alteration") +
  facet_grid(tumor_type~cell_state, scales="free", space="free") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=c("Non-HM" = "white", 
                             "Alk-assoc. HM" ="#33A02C",
                             "PDGFRA WT" ="white",
                             "PDGFRA amp" = "#E31A1C", 
                             "CDKN2A WT" ="white",
                             "CDKN2A homdel" = "#1F78B4")) +
  theme(strip.text.x = element_text(margin = margin(t = 15, unit = "pt")),
        legend.position = "top") +
  stat_n_text(size=2.25) +
  geom_text(data = pval_annotations, 
            aes(x = genetic_alt, y = y, label = label),
            inherit.aes = FALSE,
            size = 2.5,
            nudge_x = -0.3)
dev.off()
