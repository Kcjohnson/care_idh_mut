##################################
# Plot the genomic data availability across the cohort
# Author: Kevin Johnson
# Date Updated: 2026.03.30
##################################

library(tidyverse)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
source(paste0(proj_dir, "/scripts/utils/plot_theme.R"))
setwd(proj_dir)

# Curated data availability for CARE IDH-mutant cohort.
all_samples <- read.table("data/misc/genomic_data_availability_all.csv", sep = ",", header = TRUE)



# Re-level and reformat the data.
all_samples$data_type <- factor(all_samples$data_type, levels = rev(c("10x snRNA", "snATAC via 10x Multiome", "SmartSeq2 snRNA","Whole exome DNA", "Whole genome DNA","Tumor-only DNA", "Bulk RNAseq", "Bulk DNA methylation")))
all_samples$single_bulk <- ifelse(all_samples$data_type%in%c("10x snRNA", "snATAC via 10x Multiome", "SmartSeq2 snRNA"), "Single nucleus", "Bulk")
all_samples$single_bulk <- factor(all_samples$single_bulk, levels = c("Single nucleus", "Bulk"))
samples_annotated_plot <- all_samples %>% 
  mutate(available_plot = case_when(
    available=="yes" & single_bulk=="Single nucleus"~ "Single nucleus data avail.",
    available=="yes" & single_bulk=="Bulk"~ "Bulk data avail.",
    available=="no" ~ "Not avail.",
    TRUE ~ NA_character_
  ))
samples_annotated_plot$available_plot <- factor(samples_annotated_plot$available_plot, levels = c("Not avail.", "Single nucleus data avail.", "Bulk data avail."))


# Landscape of all samples
pdf(file = paste0(fig_dir, "/edf1a_sample_annotation.pdf"), height = 5, width = 9, useDingbats = FALSE)
ggplot(samples_annotated_plot, aes(x=care_id, y=data_type)) +
  geom_tile(aes(fill = factor(available_plot))) +
  scale_fill_manual(values = c("Not avail." = "gray90", "Single nucleus data avail." = "#8da0cb",
                               "Bulk data avail." = "#e78ac3",
                               "NA" = "gray70")) +
  plot_theme + 
  theme(
    panel.spacing = unit(0.15, "lines"),
    strip.text.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    legend.position = "bottom") +
  labs(y="Genomic profile", x="CARE IDHmut (n=75 samples)", fill="Data availability") +
  facet_grid(single_bulk~patient_id, scales="free", space="free") 
dev.off()

# Produce a plot for longitudinal data generation.
longitudinal_counts <- all_samples %>%
  filter(available == "yes") %>%                                
  mutate(timepoint = case_when(
    grepl("T1", care_id) ~ "S1",
    grepl("T2", care_id) ~ "S2",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(timepoint)) %>%                                 
  distinct(patient_id, data_type, timepoint) %>%               
  group_by(patient_id, data_type) %>%
  summarize(n_timepoints = n(), .groups = "drop") %>%         
  filter(n_timepoints == 2) %>%                                 
  count(data_type, name = "n_patients_longitudinal")   

# Manually add the number of patients with longitudinal DNA
new_row <- data.frame(
  data_type = "Whole exome or genome DNA",
  n_patients_longitudinal = 29)

longitudinal_counts <- rbind(longitudinal_counts, new_row)

longitudinal_counts$single_bulk <- ifelse(longitudinal_counts$data_type%in%c("10x snRNA", "snATAC via 10x Multiome", "SmartSeq2 snRNA"), "Single nucleus", "Bulk")
longitudinal_counts$single_bulk <- factor(longitudinal_counts$single_bulk, levels = c("Single nucleus", "Bulk"))

pdf(file = paste0(fig_dir, "/edf1c_longitudinal_patient_data_availability.pdf"), height = 4, width = 4, useDingbats = FALSE)
ggplot(longitudinal_counts, aes(x =reorder(data_type, -n_patients_longitudinal), 
                                y = n_patients_longitudinal)) +
  geom_bar(stat = "identity", aes(fill=single_bulk)) +
  scale_fill_manual(values = c("Single nucleus" = "#8da0cb",
                               "Bulk" = "#e78ac3")) +
  labs(x = "Genomic profile", 
       y = "Number of patients with longitudinal\ndata generated") +
  plot_theme +
  facet_grid(.~single_bulk, scales="free", space="free") +
  scale_y_continuous(breaks = seq(0, max(longitudinal_counts$n_patients_longitudinal), by = 5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  guides(fill=FALSE)
dev.off()


all_counts <- all_samples %>%
  filter(available == "yes") %>% 
  group_by(data_type) %>%
  summarize(counts = n(), .groups = "drop")

all_counts$single_bulk <- ifelse(all_counts$data_type%in%c("10x snRNA", "snATAC via 10x Multiome", "SmartSeq2 snRNA"), "Single nucleus", "Bulk")
all_counts$single_bulk <- factor(all_counts$single_bulk, levels = c("Single nucleus", "Bulk"))

pdf(file = paste0(fig_dir, "/edf1b_all_sample_data_availability.pdf"), height = 4, width = 4, useDingbats = FALSE)
ggplot(all_counts, aes(x =reorder(data_type, -counts), 
                       y = counts)) +
  geom_bar(stat = "identity", aes(fill=single_bulk)) +
  scale_fill_manual(values = c("Single nucleus" = "#8da0cb",
                               "Bulk" = "#e78ac3")) +
  labs(x = "Genomic profile", 
       y = "Number of samples with\ndata generated") +
  plot_theme +
  facet_grid(.~single_bulk, scales="free", space="free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  scale_y_continuous(breaks = seq(0, max(all_counts$counts), by = 15)) +
  guides(fill=FALSE)
dev.off()


### END ###