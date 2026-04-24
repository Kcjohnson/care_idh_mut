##################################
# Malignant cell state proportion associations with overall survival and post-recurrence survival
# Author: Kevin Johnson
# Date Updated: 2026.04.04
##################################

library(tidyverse) # 2.0.0
library(ggpubr) # 0.6.0
library(EnvStats) # 2.8.0
library(survival) # 3.5-7
library(survminer) # 0.4.9
library(forestmodel) # 0.6.2 

# Specify directories:
proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
source(paste0(proj_dir, "/scripts/utils/plot_theme.R"))
setwd(proj_dir)

### ### ### ### ### ### ### ### ###
## Clinical data set-up    ##
### ### ### ### ### ### ### ### ###

# Read in clinical genetic information for the samples present in CARE
sample_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240608.txt", sep="\t", header = TRUE)
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt", sep="\t", header = TRUE)

# Not all CARE IDH-mutant samples are from initial surgery, load in the information from all surgeries.
# The surgical_interval_mo values represents the *cumulative* time lapsed in months between surgeries.
care_sugeries_all <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/clinical/glass_clinical_surgeries.txt", sep="\t", header = TRUE)
cases <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/clinical/glass_clinical_cases.txt", sep="\t", header = TRUE)

# Select the time point 2 samples to extract the surgical interval of time lapsed until time point 2.
# Note that there is some information missing in some cases that we only have surgical interval information for T1 to T2. 
# Thus, we don't always have cumulative surgical interval.
sample_md_t2 <- sample_md %>% 
  filter(timepoint=="T2") %>% 
  inner_join(care_sugeries_all, by="sample_barcode") %>% 
  mutate(case_barcode = case_barcode.x, 
         time_lapsed_to_t2 = surgical_interval_mo.y) %>% 
  dplyr::select(case_barcode, time_lapsed_to_t2)


### ### ### ### ### ### ### ### ###
## snRNA proportions set-up     ##
### ### ### ### ### ### ### ### ###
# Malignant cell state assignment
care_state_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/scoring/caremut_malignant_cell_state_assignment.txt", sep="\t", header = TRUE)

# Relabel some of the features
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

# Determine the cell count and frequency for each of the malignant states.
care_state_freq <- care_state_md %>% 
  group_by(care_id, cell_state) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts/sum(counts)) %>% 
  ungroup() %>% 
  complete(care_id, cell_state,
           fill = list(counts = 0, freq = 0)) %>%
  inner_join(md_trim, by="care_id") 

care_state_freq_cc <- care_state_md %>% 
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

mal_care_md_state_summary <- care_state_freq %>% 
  bind_rows(care_state_freq_cc) %>% 
  inner_join(patient_md, by=c("case_barcode", "patient_id")) 

mal_care_md_state_summary$cell_state <- factor(mal_care_md_state_summary$cell_state, levels=c("Undifferentiated", "NPC-like", "OPC-like", "MES-like", "AC-like", "Cycling"))
mal_care_md_state_summary$idh_codel_subtype.x <- factor(mal_care_md_state_summary$idh_codel_subtype.x, levels=c("Oligo.", "Astro."))

# Confirm that all proportions sum to 100
ggplot(mal_care_md_state_summary %>% 
         filter(cell_state!="Cycling"), aes(x=care_id, y=freq*100, fill = cell_state)) +
  geom_bar(stat="identity")

mal_care_md_state_summary %>%
  # Ignore cycling, which is calculated separate from the core cell states.
  filter(cell_state!="Cycling") %>% 
  group_by(SampleID) %>%
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


mal_care_md_state_summary_t1t2 <- mal_care_md_state_summary %>% 
  # Keep it so that there are only 2 time points per patient
  filter(timepoint!="T3") %>% 
  inner_join(cases, by="case_barcode") %>% 
  mutate(patient_vital =  ifelse(case_vital_status=="alive", 0, 1),
         subtype = idh_codel_subtype.x) %>% 
  inner_join(sample_md_t2, by="case_barcode") %>% 
  mutate(post_t2_survival = case_overall_survival_mo-time_lapsed_to_t2)


### ### ### ### ### ### ### ### ###
## Create Kaplan-Meier and Forest plots for MES-like abundance at time point 2 (recurrence time point)
### ### ### ### ### ### ### ### ##
care_median_levels_t1t2 <- mal_care_md_state_summary_t1t2 %>% 
  group_by(subtype, cell_state, timepoint) %>% 
  mutate(median_levels = paste0("median_", ntile(freq, 2)),
         median_levels = recode(median_levels, `median_1` = "low",
                                `median_2` = "high"))  %>% 
  mutate(interval_event = 1) 
care_median_levels_t1t2$median_levels <- factor(care_median_levels_t1t2$median_levels, levels=c("low", "high"))

### MES-like ###
care_t2_mes <- care_median_levels_t1t2 %>% 
  filter(cell_state=="MES-like", timepoint=="T2") %>% 
  mutate(MES_T2_levels = recode(median_levels, `low` = "MES-like low",
                                `high` = "MES-like high"))

care_t2_mes$median_levels <- factor(care_t2_mes$median_levels, levels=c("MES-like low", "MES-like high"))

# Separate into two tumor types for independent testing
care_t2_mes_astro <- care_t2_mes %>% 
  filter(subtype=="Astro.") 

# Get a sense for the frequency values across these groups
care_t2_mes_astro %>% 
  group_by(MES_T2_levels) %>% 
  summarise(max_levels = max(freq),
            min_levels = min(freq),
            mean_levels = mean(freq))

care_t2_mes_oligo <- care_t2_mes %>% 
  filter(subtype=="Oligo.")
care_t2_mes_oligo %>% 
  group_by(MES_T2_levels) %>% 
  summarise(max_levels = max(freq),
            min_levels = min(freq),
            mean_levels = mean(freq))

dim(care_t2_mes_astro)
table(care_t2_mes_astro$MES_T2_levels)

fit_mes_t2_astro <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ MES_T2_levels,
                            data = care_t2_mes_astro)  

mes_t2_astro <- ggsurvplot(fit_mes_t2_astro, data = care_t2_mes_astro, risk.table = FALSE, pval= TRUE, 
                           surv.median.line = "hv",
                           palette = c("royalblue4", "tomato3"),
                           ylab = "Overall survival probability\nAstro. (n = 22)", xlab = "Time (months)",
                           font.x = 8,
                           font.y = 8,
                           font.tickslab = 8,
                           font.legend = 8,
                           pval.size = 2.8, 
                           font.title = 8) 

p <- mes_t2_astro
pdf(paste0(fig_dir, "fig5f_astro_recurrence_mes_median_os_km.pdf"), width = 3.25, height = 2, useDingbats = FALSE)
print(p$plot)
dev.off()

post_recur_fit_mes_astro <- survfit(Surv(post_t2_survival, patient_vital) ~ MES_T2_levels,
                                    data = care_t2_mes_astro)  

post_recur_mes_t2 <- ggsurvplot(post_recur_fit_mes_astro, data = care_t2_mes_astro, risk.table = FALSE, pval= TRUE, 
                                surv.median.line = "hv",
                                palette = c("royalblue4", "tomato3"),
                                ylab = "Post-recurrence survival\nprobability Astro. (n = 22)", xlab = "Time (months)")

p <- post_recur_mes_t2
pdf(paste0(fig_dir, "fig5f_astro_recurrence_mes_median_post_recurrence_km.pdf"), width = 5, height = 4, useDingbats = FALSE)
print(p$plot)
dev.off()

## Oligodendroglioma
dim(care_t2_mes_oligo)
table(care_t2_mes_oligo$MES_T2_levels)

fit_mes_t2_oligo <- survfit(Surv(case_overall_survival_mo, patient_vital) ~ MES_T2_levels,
                            data = care_t2_mes_oligo)  

mes_t2_oligo <- ggsurvplot(fit_mes_t2_oligo, data = care_t2_mes_oligo, risk.table = FALSE, pval= TRUE, 
                           surv.median.line = "hv",
                           palette = c("royalblue4", "tomato3"),
                           ylab = "Overall survival probability\nOligo. (n = 13)", xlab = "Time (months)",
                           font.x = 8,
                           font.y = 8,
                           font.tickslab = 8,
                           font.legend = 8,
                           pval.size = 2.8, 
                           font.title = 8)

p_oligo <- mes_t2_oligo
pdf(paste0(fig_dir, "fig5f_oligo_recurrence_mes_median_os_km.pdf"), width = 3.25, height = 2, useDingbats = FALSE, bg = "transparent")
print(p_oligo$plot)
dev.off()

post_recur_fit_mes_oligo <- survfit(Surv(post_t2_survival, patient_vital) ~ MES_T2_levels,
                                    data = care_t2_mes_oligo)  

post_recur_mes_t2_oligo <- ggsurvplot(post_recur_fit_mes_oligo, data = care_t2_mes_oligo, risk.table = FALSE, pval= TRUE, 
                                      surv.median.line = "hv",
                                      palette = c("royalblue4", "tomato3"),
                                      ylab = "Post-recurrence survival\nprobability Astro. (n = 22)", xlab = "Time (months)",
                                      font.x = 8,
                                      font.y = 8,
                                      font.tickslab = 8,
                                      font.legend = 8,
                                      pval.size = 2.8, 
                                      font.title = 8)

p <- post_recur_mes_t2_oligo
pdf(paste0(fig_dir, "fig5f_oligo_recurrence_mes_median_post_recurrence_km.pdf"), width = 5, height = 4, useDingbats = FALSE)
print(p$plot)
dev.off()


### ### ### ### ###
# Forest plot
### ### ### ### ###
all_t2_mes_filt <- care_t2_mes %>% 
  ungroup() %>% 
  # Grouping the grades this way as an alternative since Oligo. and Astro. grades are not necessarily the same.
  mutate(grade_num_t2 = paste0("G", T2)) %>% 
  dplyr::select(case_barcode, case_overall_survival_mo, post_t2_survival, treatment_t1t2, `Tumor type` = subtype, patient_vital, cell_state, `MES-like fraction` = freq,`Recur. MES level` = MES_T2_levels, 
                `Age at diagnosis` = case_age_diagnosis_years, `Recur. grade` = grade_num_t2) 


care_mes_grade_num_os <- forestmodel::forest_model(coxph(Surv(case_overall_survival_mo, patient_vital) 
                                                         ~ `Recur. MES level` + `Age at diagnosis` + `Tumor type` + `Recur. grade`,
                                                         data = all_t2_mes_filt))
care_mes_grade_num_os

pdf(paste0(fig_dir, "forest_plot_os_coxph_who_grades.pdf"), width = 9, height = 5, useDingbats = FALSE)
care_mes_grade_num_os
dev.off()


# Categorizing G2 and G3 together since there is only one G2 at recurrence.
table(all_t2_mes_filt$`Tumor type`, all_t2_mes_filt$`Recur. grade`)
care_t2_mes_astro <- all_t2_mes_filt %>% 
  filter(`Tumor type`=="Astro.") %>% 
  mutate(`Recur. grade` = ifelse(`Recur. grade`%in%c("G2","G3"), "G2/G3", "G4"))

# Post-recurrence survival
care_mes_group_idh_a_post_recur <- forestmodel::forest_model(coxph(Surv(post_t2_survival, patient_vital) 
                                                        ~ `Recur. MES level` + `Age at diagnosis` + `Recur. grade`,
                                                        data = care_t2_mes_astro))
care_mes_group_idh_a_post_recur


# Overall survival
care_mes_group_idh_a_os <- forestmodel::forest_model(coxph(Surv(case_overall_survival_mo, patient_vital) 
                                                           ~ `Recur. MES level` + `Age at diagnosis` + `Recur. grade`,
                                                           data = care_t2_mes_astro))

pdf(paste0(fig_dir, "astro_recur_mes_median_forest_plot_os_coxph.pdf"), width = 9, height = 5, useDingbats = FALSE)
care_mes_group_idh_a_os
dev.off()

# Attempt to perform coxph model with Oligo. tumors.
care_t2_mes_oligo <- all_t2_mes_filt %>% 
  filter(`Tumor type`=="Oligo.") %>% 
  dplyr::select(`Age at diagnosis`, `Recur. MES level`, case_overall_survival_mo, patient_vital)

# Not enough observations though the trend is there
care_mes_group_idh_o_os <- forestmodel::forest_model(coxph(Surv(case_overall_survival_mo, patient_vital) 
                                                           ~ `Recur. MES level` + `Age at diagnosis`,
                                                           data = care_t2_mes_oligo))

care_mes_group_idh_o_os


### END ###