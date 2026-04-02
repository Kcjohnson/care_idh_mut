##################################
# Visualize the cell type abundance frequency across tumor types and time points
# Author: Kevin Johnson
# Date Updated: 2026.03.30
##################################

library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(EnvStats)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "results/figures/rna/proportions/")
setwd(proj_dir)

source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")

# Load in the CAREmut cell type assignment metadata
caremut_md <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20260320.txt", sep = "\t", row.names = 1, header = TRUE)

# Renaming variable for clarity. CellType_final represents the integrated expression and infercnv results
celltype_index <- which(colnames(caremut_md)=="CellType_final")
colnames(caremut_md)[celltype_index] <- "CellType"

# **Remove unresolved cells from further analysis** 
caremut_md <- caremut_md %>% 
  filter(CellType!="Unresolved")

caremut_case_linker <- caremut_md %>% 
  mutate(patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2))) %>% 
  dplyr::select(case_barcode, patient_id) %>% 
  distinct()
rownames(caremut_case_linker) <- NULL
caremut_md_case <- caremut_md %>% 
  mutate(patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2))) %>% 
  dplyr::select(SampleID, lab, sample_barcode, idh_codel_subtype, care_id, patient_id, timepoint) %>% 
  distinct()
rownames(caremut_md_case) <- NULL

# Clinical and genetic summary table for samples and patients
sample_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240608.txt",  sep="\t", header = TRUE)
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt",  sep="\t", header = TRUE)


# Summarize the frequency of each CellType across samples
care_md_summary <- caremut_md %>% 
  group_by(SampleID, CellType) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  dplyr::select(SampleID, CellType, counts, freq) %>% 
  complete(SampleID, CellType,
           fill = list(counts = 0, freq = 0)) %>%
  distinct() %>% 
  inner_join(caremut_md_case, by="SampleID") %>% 
  inner_join(patient_md, by=c("patient_id", "idh_codel_subtype")) %>% 
  left_join(sample_md, by=c("sample_barcode", "case_barcode", "idh_codel_subtype", "patient_id", "timepoint", "care_id")) %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."))


# Confirm that all sample sum 1 and that all cell types are measured.
care_md_summary %>%
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

table(care_md_summary$CellType)==n_distinct(care_md_summary$SampleID)

# Inspect the relative frequency across tumor type and time point.
care_md_summary_type <- care_md_summary %>% 
  dplyr::select(idh_codel_subtype, timepoint, CellType, freq) %>% 
  group_by(CellType, timepoint, idh_codel_subtype) %>% 
  summarise(avg_freq = mean(freq))  

# Determine plotting order
cell_state_order <- c("Malignant",
                      "Oligodendrocyte",
                      "Astrocyte",
                      "ExcNeuron",
                      "InhNeuron",
                      "Lymphocyte" ,
                      "Myeloid",
                      "Endothelial",
                      "Mural")
care_md_summary$CellType <- factor(care_md_summary$CellType, levels=rev(cell_state_order))
care_md_summary$idh_codel_subtype <- factor(care_md_summary$idh_codel_subtype, levels=c("Oligo.", "Astro."))


# Stacked bar plot across all samples
care_md_summary %>% 
  ggplot(aes(x=care_id, fill = factor(CellType), y=freq*100)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values=c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                             "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
                             "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7")) +
  labs(y = "Cell proportion", fill="Cell Type", x = "") +
  theme_bw() +
  theme(panel.grid.major = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(.~idh_codel_subtype, scales = "free_x", space = "free_x") 

# Stacked bar plot for the longitudinal pairs
care_md_summary %>% 
  filter(timepoint!="T3") %>% 
  ggplot(aes(x=patient_id, fill = factor(CellType), y=freq*100)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values=c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                             "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
                             "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7")) +
  labs(y = "Cell proportion", fill="Cell Type", x = "") +
  theme(panel.grid.major = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(timepoint~idh_codel_subtype, scales = "free", space = "free") +
  plot_theme



ggplot(care_md_summary %>% 
         filter(timepoint!="T3") %>% 
         mutate(timepoint = recode(timepoint, "T1" = "I",
                                   "T2" = "R")), aes(x = timepoint, y = freq*100)) + 
  geom_line(aes(group=patient_id), color="gray70", linetype=2) +
  geom_boxplot(aes(fill=CellType)) +
  scale_linetype_manual(values="dashed") +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 4, label="p.format") +
  labs(x="Longitudinal IDHmut pairs", y="Cell abundance (%)") +
  scale_fill_manual(values=c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                             "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
                             "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7")) +
  theme(strip.background = element_blank()) +
  facet_grid(.~CellType, scales="free") +
  stat_n_text() +
  plot_theme 


pdf(paste0(fig_dir, "caremut_cell_composition_paired_t1t2.pdf"), width = 6, height = 4, useDingbats = FALSE)
ggplot(care_md_summary %>% 
         filter(timepoint!="T3") %>% 
         mutate(timepoint = recode(timepoint, "T1" = "I",
                                   "T2" = "R")), aes(x = timepoint, y = freq*100)) + 
  geom_line(aes(group=patient_id), color="gray70", linetype=2) +
  geom_boxplot(aes(fill=CellType)) +
  scale_linetype_manual(values="dashed") +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", paired = TRUE, size = 4, label="p.format") +
  labs(x="Longitudinal IDHmut pairs", y="Cell abundance (%)") +
  scale_fill_manual(values=c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                             "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
                             "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7")) +
  theme(strip.background = element_blank()) +
  facet_grid(idh_codel_subtype~CellType, scales="free") +
  stat_n_text() +
  plot_theme 
dev.off()


# Stacked bar plot
care_md_summary_avg <- care_md_summary %>%
  # Exclude time point 3 values from initial vs. recurrence comparisons
  filter(timepoint!="T3") %>%
  group_by(timepoint, idh_codel_subtype, CellType) %>% 
  summarise(avg_freq = mean(freq)) %>% 
  mutate(timepoint = recode(timepoint, "T1" = "I",
                            "T2" = "R"),
         subtype = recode(idh_codel_subtype, "IDH-O" = "Oligo.",
                          "IDH-A" = "Astro."))

care_md_summary_avg$subtype <- factor(care_md_summary_avg$subtype, levels=rev(c("Astro.","Oligo.")))

pdf(paste0(fig_dir, "caremut_celltype_stacked_barplot.pdf"), width = 5, height = 4, useDingbats = FALSE, bg="transparent")
ggplot(care_md_summary_avg, aes(x=timepoint, fill = factor(CellType), y=avg_freq*100)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values=c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                             "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
                             "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7")) +
  facet_grid(.~subtype, scales="free") +
  labs(x = "Time point", y="Mean cellular abundance (%)", fill="Cell type") +
  plot_theme
dev.off()

# Compare T1+T2 cell type abundance so that each patient only has two samples.
# Statistically significant cell types include: Mural (higher in oligodendrogliomas) and Myeloid (higher in astrocytomas)
# Borderline statistical significance includes: Oligodendrocytes (P = 0.054) and Endothelial (P = 0.085) both higher in oligodendrogliomas
ggplot(care_md_summary %>% 
         filter(timepoint!="T3"), aes(x = idh_codel_subtype, y = freq)) + 
  geom_boxplot(aes(fill=CellType)) +
  theme_bw() +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", size =2) +
  labs(x="IDHmut subtype comparisons", y="Cell proportion") +
  scale_fill_manual(values=c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                             "Astrocyte" = "#BFBADA", "ExcNeuron" = "#BC80BD", 
                             "InhNeuron" = "#FFED6F", "Endothelial" = "#FCCDE5", "Mural" = "#FFFFB3", "Lymphocyte" = "#8DD3C7")) +
  theme(strip.background = element_blank()) +
  facet_grid(.~CellType, scales="free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_n_text() + 
  plot_theme


### END ###