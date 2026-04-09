##################################
# Visualize the CIBERSORTx estimated cell type/state abundance in the Glioma Longitudinal AnalySiS (GLASS) RNA sequencing dataset across grade, time points, and stratified by genetics.
# Author: Kevin Johnson
# Date Updated: 2026.04.07
##################################

# This script generates figures related to analysis of GLASS cohort, largely using CIBERSORTx results.

library(tidyverse)
library(ggpubr)
library(EnvStats)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
base_data_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/glass"
setwd(proj_dir)
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")

# Load all relevant glass database version (September release 2025), available on Synapse.
# Some text files have slightly different format than the verhaak lab-internal glass SQL database.
aliquots <- read.table(file.path(base_data_dir, "glass_biospecimen_aliquots.txt"), sep = "\t", header = TRUE)
cases <- read.table(file.path(base_data_dir, "glass_clinical_cases.txt"), sep = "\t", header = TRUE)
surgeries <- read_tsv("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/glass/glass_clinical_surgeries.txt")
samples <- read.table(file.path(base_data_dir, "glass_clinical_samples.txt"), sep = "\t", header = TRUE)
subtypes <- read.table(file.path(base_data_dir, "glass_clinical_subtypes.txt"), sep = "\t", header = TRUE)
rna_blocklist <- read.table(file.path(base_data_dir, "glass_analysis_rna_blocklist.txt"), sep = "\t", header = TRUE)
files <-read_tsv("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/glass/glass_analysis_files.txt")
rna_silver_set <- read.table(file.path(base_data_dir, "glass_analysis_rna_silver_set.txt"), sep = "\t", header = TRUE)
tumor_rna_clinical_comparison <- read_tsv(file.path(base_data_dir, "glass_analysis_tumor_rna_clinical_comparison.txt"),
                                            col_types = cols(
                                              received_tmz             = col_character(),
                                              received_rt    = col_character(),
                                              received_alk  = col_character(),
                                              received_treatment    = col_character(),
                                              received_pd1 = col_character(),
                                              received_bev = col_character()
                                            ))
top_transcriptional_subtype   <- read.table(file.path(base_data_dir, "glass_analysis_top_transcriptional_subtype.txt"), sep = "\t", header = TRUE)
sample_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240608.txt", sep="\t", header = TRUE)
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt", sep="\t", header = TRUE)
longitudinal_genetic_df <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/glass/case_longitudinal_acquired_genetic_event_annotation.txt", sep="\t", header = TRUE)
mp_gsea <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/results/glass/caremut_malignant_plus_myeloid_metaprograms_ssgsea_20260407.txt", sep="\t", row.names = 1, header = TRUE)
colnames(mp_gsea) <- gsub("\\.","-",colnames(mp_gsea))
mp_gsea <- as.data.frame(t(mp_gsea))
mp_gsea$aliquot_barcode <- row.names(mp_gsea)
mp_gsea$case_barcode <- substr(mp_gsea$aliquot_barcode, 1, 12)
mp_gsea$sample_barcode <- substr(mp_gsea$aliquot_barcode, 1, 15)

# Keep only those samples with clear with tumor type/subtype information.
rna_aliquots <- aliquots %>% 
  filter(aliquot_analyte_type=="R") %>% 
  mutate(case_barcode = substr(aliquot_barcode, 1, 12)) %>% 
  inner_join(subtypes, by="case_barcode")

# CARE IDH-mutant cases also have GLASS IDs.
care_cases <- patient_md

## snRNA proportions ##
# Load in the CAREmut metadata - updated 2026
snrna_md <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20260320.txt", sep = "\t", row.names = 1, header = TRUE)
celltype_index <- which(colnames(snrna_md)=="CellType_final")
colnames(snrna_md)[celltype_index] <- "CellType"
snrna_md <- snrna_md %>% 
  filter(CellType!="Unresolved")

care_subtypes <- snrna_md %>% 
  dplyr::select(sample_barcode, idh_codel_subtype) %>% 
  distinct() %>% 
  mutate(case_barcode = substr(sample_barcode, 1, 12)) %>% 
  dplyr::select(case_barcode, idh_codel_subtype)

# Malignant cell state
md_malignant <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/scoring/caremut_malignant_cell_state_assignment.txt", sep="\t", header = TRUE)

# Load in CIBERSORTx results, which were applied to all GLASS RNA samples. Need to remove IDH-wildtype cases from analysis as these are not relevant for our analyses.
csx_fractions <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/cibersortx/glass_all_cibersortx_celltype_state_fractions.txt", header = TRUE, sep = "\t")
csx_fractions$aliquot_barcode <- gsub("\\.", "-", csx_fractions$Mixture)
csx_fractions$sample_barcode <- substr(csx_fractions$aliquot_barcode, 1, 15)
csx_fractions$sample_barcode[which(duplicated(csx_fractions$sample_barcode))]

snrna_md_proportions_wide <- snrna_md %>% 
  mutate(CellType = recode(CellType, `ExcNeuron` = "Neuron",
                           `InhNeuron` = "Neuron",
                           `Mural` = "MuralEndothelial",
                           `Endothelial` = "MuralEndothelial")) %>% 
  group_by(sample_barcode, CellType) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  dplyr::select(sample_barcode, CellType, counts, freq) %>% 
  complete(sample_barcode, CellType,
           fill = list(counts = 0, freq = 0)) %>%
  distinct() %>% 
  dplyr::select(sample_barcode, cell_type = CellType, freq) %>% 
  pivot_wider(names_from = cell_type, values_from = freq) %>% 
  filter(sample_barcode%in%csx_fractions$sample_barcode) %>% 
  as.data.frame()

csx_fractions_wide <- csx_fractions %>% 
  filter(sample_barcode%in%snrna_md_proportions_wide$sample_barcode) %>% 
  mutate(Malignant = NPC.like+OPC.like+Undifferentiated+AC.like+MES.like) %>% 
  dplyr::select(colnames(snrna_md_proportions_wide))

cell_type <- colnames(csx_fractions_wide)[2:8]
num_columns <- length(csx_fractions_wide)
correlation_results <- data.frame()
all(snrna_md_proportions_wide$sample_barcode==csx_fractions_wide$sample_barcode)

# Examine the correlations and p-values for Pearson's correlation coefficient
for (i in 2:num_columns) {
  col1 <- snrna_md_proportions_wide[, i]
  col2 <- csx_fractions_wide[, i]
  correlation_coefficient <- cor(col1, col2, method = "pearson")
  correlation_pvalue <- cor.test(col1, col2, method = "pearson")$p.value
  correlation_results <- rbind(correlation_results, data.frame(snRNA_cell_type = colnames(snrna_md_proportions_wide)[i], bulk_cell_type = colnames(csx_fractions_wide)[i], Correlation = correlation_coefficient, Pvalue = correlation_pvalue))
}
correlation_results

snrna_md_proportions_long <- snrna_md %>% 
  mutate(CellType = recode(CellType, `ExcNeuron` = "Neuron",
                           `InhNeuron` = "Neuron",
                           `Mural` = "MuralEndothelial",
                           `Endothelial` = "MuralEndothelial")) %>% 
  group_by(sample_barcode, CellType) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  dplyr::select(sample_barcode, CellType, counts, freq) %>% 
  complete(sample_barcode, CellType,
           fill = list(counts = 0, freq = 0)) %>%
  dplyr::select(sample_barcode, cell_type = CellType, freq) %>% 
  mutate(data_type = "snrna") %>% 
  filter(sample_barcode%in%csx_fractions$sample_barcode)

csx_fractions_long <- csx_fractions %>% 
  filter(sample_barcode%in%snrna_md_proportions_long$sample_barcode) %>%
  mutate(Malignant = NPC.like+OPC.like+Undifferentiated+AC.like+MES.like) %>% 
  dplyr::select(sample_barcode, Oligodendrocyte, Neuron, Myeloid, Astrocyte, Lymphocyte, MuralEndothelial, Malignant) %>% 
  pivot_longer(cols= c(Oligodendrocyte, Neuron, Myeloid, Astrocyte, Lymphocyte, MuralEndothelial, Malignant),
               names_to = "cell_type",
               values_to = "freq") %>% 
  mutate(data_type = "bulk_rna")

# Confirm that all sample sum 1
snrna_md_proportions_long %>%
  group_by(sample_barcode) %>%
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

csx_fractions_long %>%
  group_by(sample_barcode) %>%
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


pairwise_fractions <- csx_fractions_long %>% 
  bind_rows(snrna_md_proportions_long) %>% 
  pivot_wider(names_from = data_type, values_from = freq)  %>% 
  distinct()
table(pairwise_fractions$sample_barcode)

pairwise_fractions$cell_type <- factor(pairwise_fractions$cell_type, levels=c("Malignant", "Oligodendrocyte", "Myeloid", "Neuron", "MuralEndothelial", "Astrocyte", "Lymphocyte"))

ggplot(pairwise_fractions, aes(x=snrna*100, y=bulk_rna*100)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  stat_cor(method="pearson") + 
  plot_theme +
  facet_wrap(.~cell_type, scales="free", nrow = 2) +
  labs(x="snRNA cell proportions (%)", y="Bulk CIBERSORTx proportions (%)",  title="IDHmut (n=38) samples with matched snRNA + bulk RNA")

pairwise_fractions$cell_type <- factor(pairwise_fractions$cell_type, levels=c("Malignant", "Oligodendrocyte", "Myeloid", "Neuron"))

pdf(paste0(fig_dir, "edf3d_matched_snrna_bulk_rna_n38_pearson.pdf"), width = 4.5, height = 4, useDingbats = FALSE)
ggplot(pairwise_fractions %>% 
         filter(cell_type%in%c("Malignant", "Oligodendrocyte", "Myeloid", "Neuron")), aes(x=snrna*100, y=bulk_rna*100, color=cell_type)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  scale_color_manual(values=c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3",  "Neuron" = "#BC80BD")) +
  stat_cor(method="pearson") + 
  plot_theme +
  theme(legend.position = "bottom") +
  labs(x="snRNA abundance (%)", y="Bulk CIBERSORTx abundance (%)",  
       color="") +
  guides(color=guide_legend(nrow=1,byrow=TRUE))
dev.off()


# Add malignant compartment annotation to inspect CIBERSORTx performance
md_all_plus_malignant <- snrna_md %>% 
  left_join(md_malignant, by=c("CellID", "SampleID", "case_barcode", "idh_codel_subtype", "care_id")) %>% 
  mutate(MalignantState = recode(State, `MP_AC1_MUT` = "AC-like",
                                 `MP_OPC_MUT` = "OPC-like",
                                 `MP_NPC_MUT` = "NPC-like",
                                 `MP_MES_MUT` = "MES-like",
                                 `MP_AC2_MUT` = "AC-like",
                                 "Undifferentiated" = "Undifferentiated")) %>% 
  mutate(patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2))) 


# Creating a new variable to further split Malignant cell type into more granular cell states
md_all_plus_malignant$CellTypeRefined <- md_all_plus_malignant$CellType
md_all_plus_malignant$CellTypeRefined[md_all_plus_malignant$CellTypeRefined=="Malignant"] <- md_all_plus_malignant$MalignantState[md_all_plus_malignant$CellTypeRefined=="Malignant"]
#  The few cells that were not included in scoring are dropped from analysis
sample_to_drop <- md_all_plus_malignant$SampleID[is.na(md_all_plus_malignant$CellTypeRefined)][1]

message("Dropping the following sample since it does not have malignant cells classified: ", sample_to_drop)

# Summarise the frequency of each CellType across samples
snrna_md_all_proportions_long <- md_all_plus_malignant %>% 
  filter(SampleID != sample_to_drop) %>% 
  mutate(CellTypeRefined = recode(CellTypeRefined, `ExcNeuron` = "Neuron",
                                  `InhNeuron` = "Neuron",
                                  `Mural` = "MuralEndothelial",
                                  `Endothelial` = "MuralEndothelial")) %>% 
  group_by(sample_barcode, CellTypeRefined) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  dplyr::select(sample_barcode, CellTypeRefined, counts, freq) %>% 
  complete(sample_barcode, CellTypeRefined,
           fill = list(counts = 0, freq = 0)) %>%
  dplyr::select(sample_barcode, cell_type = CellTypeRefined, freq) %>% 
  mutate(data_type = "snrna") %>% 
  filter(sample_barcode%in%csx_fractions$sample_barcode)

# Sanity checks.
table(snrna_md_all_proportions_long$sample_barcode)
snrna_md_all_proportions_long %>%
  group_by(sample_barcode) %>%
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

# This should be n = 37 samples due to dropping SJ02-3, which has too few cells to be considered reliable for estimating proportions.
csx_all_fractions_long <- csx_fractions %>% 
  filter(sample_barcode%in%snrna_md_all_proportions_long$sample_barcode) %>%
  dplyr::select(sample_barcode, AC.like:Lymphocyte) %>% 
  pivot_longer(cols= c(AC.like:Lymphocyte),
               names_to = "cell_type",
               values_to = "freq") %>% 
  mutate(data_type = "bulk_rna") %>% 
  mutate(cell_type = gsub("\\.", "-", cell_type))

pairwise_fractions <- csx_all_fractions_long %>% 
  bind_rows(snrna_md_all_proportions_long) %>% 
  pivot_wider(names_from = data_type, values_from = freq)  %>% 
  distinct()

ggplot(pairwise_fractions, aes(x=snrna, y=bulk_rna)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  stat_cor(method="pearson") + 
  plot_theme +
  facet_wrap(.~cell_type, scales="free", nrow = 3) +
  labs(x="snRNA cell proportions", y="Bulk CIBERSORTx proportions",  title="IDHmut (n=37) samples with matched snRNA + bulk RNA")

pairwise_fractions_select <- pairwise_fractions %>% 
  filter(cell_type%in%c("Undifferentiated","AC-like", "NPC-like", "OPC-like", "MES-like"))
pairwise_fractions_select$cell_type <- factor(pairwise_fractions_select$cell_type, levels=c("AC-like", "MES-like", "Undifferentiated", "OPC-like", "NPC-like"))
n_distinct(pairwise_fractions_select$sample_barcode)

pdf(paste0(fig_dir, "edf7b_matched_snrna_bulk_rna_malignant_correlations_n37.pdf"), width = 4, height = 4, useDingbats = FALSE)
ggplot(pairwise_fractions_select %>% 
         filter(cell_type%in%c("Undifferentiated","AC-like", "NPC-like", "OPC-like", "MES-like")), aes(x=snrna*100, y=bulk_rna*100, color=cell_type)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  stat_cor(method="pearson", size = 2.25) + 
  plot_theme +
  scale_color_manual(values=c("AC-like" = "#AA2756", 
                              "MES-like"="#F77D58",
                              "NPC-like" = "#7fbf7b",
                              "OPC-like"="#E8F5A3",
                              "Undifferentiated" = "gray70")) +
  guides(color=FALSE) +
  facet_wrap(.~cell_type, scales="free", nrow = 2) +
  labs(x="snRNA cell proportions (%)", y="Bulk CIBERSORTx proportions (%)") 
dev.off()

# Information on whether the patient received treatment or had a tumor grade change at recurrent timepoint
tumor_rna_clinical_comparison_filt <- tumor_rna_clinical_comparison %>% 
  dplyr::select(case_barcode, received_treatment, received_rt, grade_change)

surgeries_grade <- surgeries %>% 
  mutate(grade_num = case_when(grade=="II" ~ "G2",
                               grade=="III" ~ "G3",
                               grade=="IV" ~ "G4",
                               TRUE ~ NA_character_)) %>% 
  dplyr::select(sample_barcode, grade_num, surgery_number, surgical_interval_mo) 


rna_silver_samples <- rna_silver_set %>% 
  inner_join(subtypes, by="case_barcode") %>% 
  dplyr::select(case_barcode:tumor_barcode_b) %>% 
  pivot_longer(-case_barcode) %>% 
  mutate(timepoint =ifelse(name=="tumor_barcode_a", "T1", "T2")) %>% 
  dplyr::select(case_barcode, aliquot_barcode = value, timepoint) %>% 
  inner_join(top_transcriptional_subtype, by="aliquot_barcode") %>% 
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15)) %>% 
  inner_join(surgeries_grade, by="sample_barcode") %>% 
  inner_join(tumor_rna_clinical_comparison_filt, by="case_barcode") %>% 
  inner_join(subtypes, by="case_barcode") %>% 
  # Do not include IDH-wildtype or any CARE IDH-mutant cases where we are already measuring by snRNA
  filter(idh_codel_subtype!="IDHwt", !case_barcode%in%care_cases$case_barcode) 

mut_csx_fractions <- rna_silver_samples %>% 
  inner_join(csx_fractions, by="sample_barcode") %>% 
  # Multiple aliquots from the same timepoint (i.e., "02R").
  filter(aliquot_barcode.y!="GLSS-CU-P101-R1-02R-RNA-8AOGKJ") 

mut_csx_fractions_long <- mut_csx_fractions %>% 
  dplyr::select(-c(Mixture, P.value, Correlation, RMSE)) %>% 
  pivot_longer(cols= c(AC.like:Lymphocyte),
               names_to = "cell_type",
               values_to = "freq") %>%
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-noncodel` = "Astro.",
                                    `IDHmut-codel` = "Oligo."))
mut_csx_fractions_long$idh_codel_subtype <- factor(mut_csx_fractions_long$idh_codel_subtype, levels=c("Oligo.", "Astro."))

ggplot(mut_csx_fractions_long %>% 
         filter(cell_type%in%c("Neuron", "Myeloid", "Lymphocyte", "Oligodendrocyte")), aes(x=idh_codel_subtype, y=freq*100)) +
  geom_boxplot(aes(fill=cell_type)) +
  stat_compare_means(method="wilcox", label="p.format") + 
  plot_theme +
  scale_fill_manual(values=c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3",  "Neuron" = "#BC80BD")) +
  facet_grid(.~cell_type, scales="free") +
  stat_n_text() +
  theme(legend.position = "none") +
  labs(y="CIBERSORTx abundance (%)", x = "Subtype")


rna_silver_samples <- rna_silver_set %>% 
  dplyr::select(case_barcode:tumor_barcode_b) %>% 
  pivot_longer(-case_barcode) %>% 
  mutate(timepoint =ifelse(name=="tumor_barcode_a", "T1", "T2")) %>% 
  dplyr::select(case_barcode, aliquot_barcode = value, timepoint) %>% 
  inner_join(top_transcriptional_subtype, by="aliquot_barcode") %>% 
  mutate(sample_barcode = substr(aliquot_barcode, 1, 15)) %>% 
  inner_join(surgeries_grade, by="sample_barcode") %>% 
  inner_join(tumor_rna_clinical_comparison_filt, by="case_barcode")  %>% 
  inner_join(subtypes, by="case_barcode") %>% 
  # Do not include IDH-wildtype or any CARE IDH-mutant cases where we are already measuring by snRNA
  filter(idh_codel_subtype!="IDHwt", !case_barcode%in%care_cases$case_barcode)


all_mut_csx_fractions <- rna_silver_samples %>% 
  inner_join(csx_fractions, by="sample_barcode") %>% 
  # Remove the aliquot_barcode where the sample_barcode is duplicated. Here it's the second portion.
  filter(aliquot_barcode.y!="GLSS-CU-P101-R1-02R-RNA-8AOGKJ") 

all_mut_csx_fractions_long <- all_mut_csx_fractions %>% 
  dplyr::select(-c(Mixture, P.value, Correlation, RMSE)) %>% 
  pivot_longer(cols= c(AC.like:Lymphocyte),
               names_to = "cell_type",
               values_to = "freq") 

# Comparing all cell types. One key observation is that MuralEndothelial increases at Recurrence. This may be due to MES-like malignant signature. That is, picking up the same signal.
ggplot(all_mut_csx_fractions_long, aes(x = timepoint, y = freq)) + 
  geom_line(aes(group=case_barcode), color="gray70", linetype=2) +
  geom_boxplot(aes(fill=cell_type)) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  labs(x="Matched longitudinal analysis", y="Cell proportion") +
  theme(strip.background = element_blank()) +
  facet_grid(.~cell_type, scales="free") +
  stat_n_text(size = 2.25)

# Need to combine ALL malignant cell states in order to reaction the total percentage of malignant cells.
# We can then perform analyses that represent the relative fraction of malignant cells.
# This is important because we don't want sample purity to impact associations. A tumor with 50% MES-like in a tumor with 20% purity would be diluted out.
all_mut_csx_fractions_malignant <- all_mut_csx_fractions %>% 
  mutate(Malignant = OPC.like+Undifferentiated+AC.like+MES.like+NPC.like) %>% 
  dplyr::select(sample_barcode, Malignant)

# Calculating the relative percentage.
all_mut_csx_fractions_malignant_adj <- all_mut_csx_fractions %>%  
  dplyr::select(-c(P.value, Correlation, RMSE)) %>% 
  pivot_longer(cols= c(AC.like:Lymphocyte),
               names_to = "cell_type",
               values_to = "freq") %>% 
  inner_join(all_mut_csx_fractions_malignant, by="sample_barcode") %>% 
  filter(cell_type%in%c("OPC.like", "Undifferentiated", "AC.like", "MES.like", "NPC.like")) %>% 
  mutate(malignant_freq = freq/Malignant) %>% 
  mutate(cell_type = gsub("\\.", "-", cell_type)) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-noncodel` = "Astro.",
                                    "IDHmut-codel" = "Oligo.")) 

# Setting the order by which the cell types and tumor types appear.
all_mut_csx_fractions_malignant_adj$cell_type <- factor(all_mut_csx_fractions_malignant_adj$cell_type, levels=c("AC-like", "MES-like", "NPC-like", "OPC-like", "Undifferentiated"))
all_mut_csx_fractions_malignant_adj$idh_codel_subtype <- factor(all_mut_csx_fractions_malignant_adj$idh_codel_subtype, levels=c("Oligo.", "Astro."))

# Plot with all GLASS longitudinal RNA data. This EXCLUDES CARE samples.
unique(all_mut_csx_fractions_malignant_adj$case_barcode)

pdf(paste0(fig_dir, "edf8i_all_available_glass_rna_n65.pdf"), width = 3.5, height = 3, useDingbats = FALSE)
ggplot(all_mut_csx_fractions_malignant_adj %>% 
         filter(cell_type%in%c("AC-like", "OPC-like", "NPC-like", "MES-like", "Undifferentiated"), !case_barcode%in%care_cases$case_barcode), aes(x = timepoint, y = malignant_freq*100)) + 
  geom_line(aes(group=case_barcode), color="gray70", linetype=2) +
  geom_boxplot(aes(fill=cell_type), outlier.shape = NA) +
  geom_point(size = 0.6) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  labs(x="Matched longitudinal analysis", y="GLASS CIBERSORTx relative\nmalignant cell proportion (%)") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  theme(strip.background = element_blank()) +
  facet_grid(.~cell_type, scales="free") +
  stat_n_text(size = 2.25) +
  scale_x_discrete(labels = c("Initial" = "Init.",
                              "Recurrence" = "Rec."))
dev.off()

#### Tumor grade 
pdf(paste0(fig_dir, "glass_grade_rna_n130_wilcox.pdf"),  width = 4, height = 3, useDingbats = FALSE)
ggplot(all_mut_csx_fractions_malignant_adj %>% 
         filter(cell_type%in%c("AC-like", "OPC-like", "NPC-like", "MES-like", "Undifferentiated"), !case_barcode%in%care_cases$case_barcode), aes(x = grade_num, y = malignant_freq*100)) + 
  geom_boxplot(aes(fill=cell_type), outlier.shape = NA) +
  geom_point(size = 0.6) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", paired = FALSE, size = 2.25, label="p.format") +
  labs(x="GLASS IDHmut bulk RNA samples (n = 130)", y="GLASS CIBERSORTx relative\nmalignant cell proportion (%)") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  theme(strip.background = element_blank()) +
  facet_grid(idh_codel_subtype~cell_type, scales="free") +
  stat_n_text(size = 2.25) +
  ylim(-10, 100)
dev.off()

pdf(paste0(fig_dir, "edf7e_glass_grade_rna_n130_kruskal.pdf"),  width = 4, height = 3, useDingbats = FALSE)
ggplot(all_mut_csx_fractions_malignant_adj %>% 
         filter(cell_type%in%c("AC-like", "OPC-like", "NPC-like", "MES-like", "Undifferentiated"), !case_barcode%in%care_cases$case_barcode), aes(x = grade_num, y = malignant_freq*100)) + 
  geom_boxplot(aes(fill=cell_type), outlier.shape = NA) +
  geom_point(size = 0.6) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal", size = 2.25, label="p.format") +
  labs(x="GLASS IDHmut bulk RNA samples (n = 130)", y="GLASS CIBERSORTx relative\nmalignant cell proportion (%)") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  theme(strip.background = element_blank()) +
  facet_grid(idh_codel_subtype~cell_type, scales="free") +
  stat_n_text(size = 2.25) +
  ylim(-10, 100)
dev.off()


#### Longitudinal genetic changes
# Add whether the cases also had evidence for one of the recurrence-assoc. genetic alterations.
all_mut_csx_fractions_malignant_adj_genetics <- all_mut_csx_fractions %>%  
  dplyr::select(-c(P.value, Correlation, RMSE)) %>% 
  pivot_longer(cols= c(AC.like:Lymphocyte),
               names_to = "cell_type",
               values_to = "freq") %>% 
  inner_join(all_mut_csx_fractions_malignant, by="sample_barcode") %>% 
  filter(cell_type%in%c("OPC.like", "Undifferentiated", "AC.like", "MES.like", "NPC.like")) %>% 
  mutate(malignant_freq = freq/Malignant) %>% 
  mutate(cell_type = gsub("\\.", "-", cell_type)) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-noncodel` = "Astro.",
                                    "IDHmut-codel" = "Oligo.")) %>% 
  inner_join(longitudinal_genetic_df, by="case_barcode") %>% 
  mutate(acquired_genetic_alt_S1S2 = recode(acquired_genetic_alt_S1S2, 
                                            `genetic_alt` = "Acquired genetic alt.",
                                            `no_genetic_alt` = "Did not acquire alt."))

# 57 patients with overlapping DNA/RNA from the same tumor samples.
n_distinct(all_mut_csx_fractions_malignant_adj_genetics$case_barcode)

pdf(paste0(fig_dir, "fig3e_glass_longitudinal_acquired_genetic_alt.pdf"), width = 2.75, height = 3, useDingbats = FALSE)
ggplot(all_mut_csx_fractions_malignant_adj_genetics %>% 
         filter(cell_type%in%c("AC-like", "Undifferentiated"), !case_barcode%in%care_cases$case_barcode, !is.na(acquired_genetic_alt_S1S2)), aes(x = timepoint, y = malignant_freq*100)) + 
  geom_boxplot(aes(fill=cell_type), outlier.shape = NA) +
  geom_line(aes(group=case_barcode), color="gray70", linetype=2, size = 0.6) +
  geom_point(size = 0.6) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 2.25, label="p.format") +
  labs(x="Matched longitudinal analysis", y="GLASS CIBERSORTx relative\nmalignant cell abundance (%)") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  theme(strip.background = element_blank()) +
  facet_grid(acquired_genetic_alt_S1S2~cell_type, scales="free") +
  stat_n_text(size = 2.25)  +
  scale_x_discrete(labels = c("T1" = "Init.",
                              "T2" = "Recur."))
dev.off()


all_mut_csx_fractions_malignant_adj_wide <- all_mut_csx_fractions_malignant_adj_genetics %>% 
  dplyr::select(case_barcode, timepoint, cell_type, malignant_freq, acquired_genetic_alt_S1S2, idh_codel_subtype, received_treatment) %>% 
  pivot_wider(names_from = c(timepoint, cell_type), values_from = malignant_freq) %>% 
  mutate(`AC-like` = `T2_AC-like`-`T1_AC-like`, 
         `Undifferentiated` = `T2_Undifferentiated`-`T1_Undifferentiated`,
         `MES-like` = `T2_MES-like`-`T1_MES-like`) %>% 
  dplyr::select(case_barcode:received_treatment,  `Undifferentiated`, `AC-like`, `MES-like`) %>% 
  pivot_longer(cols = c(`Undifferentiated`, `AC-like`), names_to = "cell_type", values_to = "percentage_change") %>% 
  mutate(acquired_genetic_alt_S1S2 = recode(acquired_genetic_alt_S1S2, `genetic_alt` = "Yes",
                                            `no_genetic_alt` = "No"),
         received_treatment = recode(received_treatment, `1` = "Yes",
                                     `0` = "No")) 

genetic_alt_vs_treatment_df <- all_mut_csx_fractions_malignant_adj_wide %>% 
  dplyr::select(case_barcode, received_treatment, acquired_genetic_alt_S1S2) %>% 
  distinct()
# Statistically significant association between acquired genetic alt. and treatment. P = 0.006
fisher.test(table(genetic_alt_vs_treatment_df$acquired_genetic_alt_S1S2, genetic_alt_vs_treatment_df$received_treatment))

pdf(paste0(fig_dir, "fig3g_glass_longitudinal_differences_by_genetics.pdf"), width = 3, height = 3, useDingbats = FALSE)
ggplot(all_mut_csx_fractions_malignant_adj_wide, aes(x=acquired_genetic_alt_S1S2, y=percentage_change*100)) +
  geom_boxplot(aes(fill=cell_type), outlier.shape = NA) +
  geom_point(size = 0.6) +
  facet_grid(.~cell_type, scales="free") +
  stat_compare_means(method = "wilcox", size = 2.25, label="p.format") +
  stat_n_text(size = 2.25) +
  plot_theme +
  theme(legend.position = "none") +
  labs(x="Acquired recurrence-associated genetic change?", y="GLASS CIBERSORTx longitudinal\nmalignant cell change (%)") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  theme(strip.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = 2)
dev.off()

# Tidy up the GSEA scores.
mp_gsea_filt <- mp_gsea %>% 
  dplyr::select(aliquot_barcode, MP_CC_malignant)

glass_long_ssgsea <- rna_silver_set %>% 
  inner_join(subtypes, by="case_barcode") %>% 
  # Do not include IDH-wildtype or any CARE IDH-mutant cases where we are already measuring by snRNA
  filter(idh_codel_subtype!="IDHwt", !case_barcode%in%care_cases$case_barcode) %>% 
  inner_join(mp_gsea_filt, by=c("tumor_barcode_a"="aliquot_barcode")) %>% 
  inner_join(mp_gsea_filt, by=c("tumor_barcode_b"="aliquot_barcode")) %>% 
  dplyr::select(case_barcode, MP_CC_malignant.x, MP_CC_malignant.y) %>% 
  pivot_longer(cols=c(MP_CC_malignant.x, MP_CC_malignant.y), names_to = "timepoint", values_to = "ssGSEA") %>% 
  mutate(timepoint = recode(timepoint, `MP_CC_malignant.x` = "Init.",
                                `MP_CC_malignant.y` = "Recur.")) %>% 
    mutate(paired = rep(1:(n()/2), each=2)) 

glass_long_ssgsea_genetics <- glass_long_ssgsea %>% 
  inner_join(longitudinal_genetic_df, by = "case_barcode") %>% 
  mutate(acquired_genetic_alt_S1S2 = recode(acquired_genetic_alt_S1S2, 
                                            `genetic_alt` = "Acquired genetic alt.",
                                            `no_genetic_alt` = "Did not acquire alt."))
         
pdf(paste0(fig_dir, "fig3f_glass_cycling_differences_by_genetics.pdf"), width = 1.75, height = 3, useDingbats = FALSE)
ggplot(glass_long_ssgsea_genetics, aes(x = timepoint, y =ssGSEA)) + 
  geom_boxplot(outlier.shape = NA, fill="#6BAED6") +
  geom_line(aes(group=paired), color="gray70", linetype=2, size = 0.6) +
  geom_point(size = 0.6) +
  scale_linetype_manual(values="dashed") +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.format", method = "wilcox", paired = TRUE, size = 2.25) +
  labs(x="Matched longitudinal analysis", y="GLASS RNA ssGSEA - cell cycle score") +
  theme(strip.background = element_blank()) +
  stat_n_text(size = 2.25) +
  facet_grid(acquired_genetic_alt_S1S2~., scales="free") 
dev.off()


### END ###