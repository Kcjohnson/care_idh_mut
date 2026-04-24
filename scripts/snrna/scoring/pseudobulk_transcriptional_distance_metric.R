##################################
# Calculate the pseudobulk transcriptional distance between T1 and T2 from snRNA data
# Author: Kevin Johnson
##################################


library(tidyverse)
library(Matrix)
library(Seurat)
library(ggpubr)
library(EnvStats)
library(scales)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
out_data_dir <- file.path(proj_dir, "processed_data/rna/")
source(paste0(proj_dir, "/scripts/utils/generate_matched_exp_profiles.R"))
source(paste0(proj_dir, "/scripts/utils/plot_theme.R"))
setwd(proj_dir)

# Load in the CAREmut UMI produced at the beginning of the project.
umi_data_all <- readRDS("data/snrna/care_mut_umi_data_all_20230729.RDS")

# Read in the cell state classification based on p-value assignment.
care_state_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/scoring/caremut_malignant_cell_state_assignment.txt", sep="\t", header = TRUE)

# Relabel some of the features.
mdata <- care_state_md %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."),
         patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2)), 
         cell_state = recode(State, `MP_AC1_MUT` = "AC-like",
                             `MP_OPC_MUT` = "OPC-like",
                             `MP_NPC_MUT` = "NPC-like",
                             `MP_MES_MUT` = "MES-like",
                             `MP_AC2_MUT` = "AC-like",
                             "Undifferentiated" = "Undifferentiated")) %>% 
  dplyr::select(CellID, Sample = SampleID, Patient = patient_id, Timepoint = timepoint, State = cell_state) %>% 
  # Only consider T1-T2 (Initial vs. Recurrence for this analysis)
  filter(Timepoint!="T3")

# Purpose is to determine a high-level metric to assess transcriptional differences or distance across time points while controlling for malignant state
patients <- unique(mdata$Patient)
states <- c("AC-like", "MES-like", "NPC-like", "OPC-like", "Undifferentiated")
min_cells <- 25
min_exp <- 4

#  Configure this to all genes or chose to remove the pseudogenes, anti-sense genes, etc. to make computation less noisy
junk_genes <- c(rownames(umi_data_all[[1]])[grep("\\.", rownames(umi_data_all[[1]]))], rownames(umi_data_all[[1]])[grep("-AS*", rownames(umi_data_all[[1]]))], rownames(umi_data_all[[1]])[grep("LINC", rownames(umi_data_all[[1]]))],  rownames(umi_data_all[[1]])[grep("^RP[S|L]", rownames(umi_data_all[[1]]))], rownames(umi_data_all[[1]])[grep("^MT-", rownames(umi_data_all[[1]]))])
valid_genes <- rownames(umi_data_all[[1]])[!rownames(umi_data_all[[1]]) %in% junk_genes]
length(valid_genes)

set.seed(123)
t2_vs_t1_matched_exp_profiles <- generate_matched_exp_profiles(meta_data = mdata, exp_data = umi_data_all, patients = patients, states = states, min_cells = min_cells, min_exp = min_exp, valid_genes = valid_genes)

observed_dist <- t2_vs_t1_matched_exp_profiles$profiles

observed_dist <- observed_dist %>%
  group_by(Patient) %>%
  summarise(Dist = euclidean_dist(T1, T2)) %>% 
  # Rescaled to be between 0 and 1
  mutate(Dist01 = rescale(Dist)) %>% 
  mutate(patient_id = Patient, transcriptional_dist = Dist, scaled_transcriptional_dist = Dist01) %>% 
  dplyr::select(patient_id, transcriptional_dist, scaled_transcriptional_dist)

# Output the table
write.table(observed_dist, paste0(out_data_dir, "pseudobulk_transcriptional_euclidean_distance.txt"), sep="\t", row.names=FALSE)

