##################################
# Analyze malignant cell abundance across tumor types and grade
# Author: Kevin Johnson
##################################

library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(openxlsx)
library(EnvStats)
library(ggrepel)
library(cowplot)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
source(paste0(proj_dir, "/scripts/utils/plot_theme.R"))
setwd(proj_dir)


# Selected genetic events - IDH1, CDK4, PDGFRA, CDKN2A plus relevant genetic events for Oligo. drivers were also added
cbio_genetics <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/public/cbioportal_tcga_driver_alterations_across_samples.txt", header = TRUE)
cbio_cna <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/public/tcga_lgggbm_2016_select_cna.txt", header = TRUE)

# Restricting to the genetics for IDH mutation.  We're only looking at IDH1 mutations here. For actual analyses, we are using the clinical metadata for tumor type.
cbio_genetics_idh <- cbio_genetics %>% 
  filter(IDH1..MUT!="no alteration") %>% 
  mutate(`CDKN2A del` = ifelse(CDKN2A..HOMDEL=="HOMDEL (driver)", "homdel", "no alt."),
         `PDGFRA amp` = ifelse(PDGFRA..AMP=="AMP (driver)", "amp", "no alt."),
         `CDK4 amp` = ifelse(CDK4..AMP=="AMP (driver)", "amp", "no alt."),
         `CIC mut` = ifelse(CIC..MUT=="no alteration", "no alt.", "mut"),
         `FUBP1 mut` = ifelse(FUBP1..MUT=="no alteration", "no alt.", "mut"),
         `NOTCH1 mut` = ifelse(NOTCH1..MUT=="no alteration", "no alt.", "mut"),
         `PIK3CA mut` = ifelse(PIK3CA..MUT=="no alteration", "no alt.", "mut"),
         `PIK3R1 mut` = ifelse(PIK3R1..MUT=="no alteration", "no alt.", "mut"),
         `NF1 mut` = ifelse(NF1..MUT=="no alteration", "no alt.", "mut"))%>% 
  dplyr::select(Sample.ID, idh1 = IDH1..MUT, `CDKN2A del`:`NF1 mut`)
cbio_genetics_idh$idh1 <- gsub(" \\(driver\\)", "", cbio_genetics_idh$idh1)

cbio_cna <- cbio_cna %>% 
  mutate(cdkn2a = ifelse(CDKN2A=="-2", "homdel", "no alt."),
         pdgfra = ifelse(PDGFRA=="2", "amp", "no alt."),
         cdk4 = ifelse(CDK4=="2", "amp", "no alt."))
table(cbio_cna$cdkn2a, cbio_cna$pdgfra)
table(cbio_genetics_idh$`PDGFRA amp`, cbio_genetics_idh$`CDKN2A del`)

#############################
## TCGA
#############################
tcga_clinical <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/public/tcga_lgggbm_2016_clinical_data.txt", sep = "\t", header = TRUE)
table(tcga_clinical$IDH.codel.subtype)
table(tcga_clinical$IDH.status)

# I applied Qhianghu's RNA classifier to the TCGA expression data
rna_class <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/public/tcga_lgggbm_2016_bulk_rna_classification_top_subtype.txt", sep = "\t", header = TRUE)
rna_class$Sample.ID <- substr(rna_class$aliquot_barcode, 1, 15)

tcga_clinical_idhmut <- tcga_clinical %>% 
  # left_join since not all samples may have been included.
  left_join(rna_class, by="Sample.ID") %>% 
  dplyr::select(Patient.ID, Sample.ID, IDH.codel.subtype, Grade = Neoplasm.Histologic.Grade, Age = Diagnosis.Age, overall_survival_mo = Overall.Survival..Months., vital_status = Overall.Survival.Status, signature_name:p_value) %>% 
  mutate(vital_status = as.numeric(sapply(strsplit(vital_status, ":"), "[[", 1))) %>% 
  distinct() %>% 
  filter(IDH.codel.subtype%in%c("IDHmut-codel", "IDHmut-non-codel")) %>% 
  mutate(IDH.codel.subtype = recode(IDH.codel.subtype, `IDHmut-non-codel` = "Astro.",
                                    `IDHmut-codel` = "Oligo.")) %>% 
  left_join(cbio_genetics_idh, by=c("Sample.ID"="Sample.ID")) %>% 
  mutate(revised_grade = ifelse(Grade%in%c("G2", "G3") & `CDKN2A del`=="homdel"&IDH.codel.subtype=="Astro.", "G4", Grade)) %>% 
  mutate(revised_grade = ifelse(is.na(revised_grade) & `CDKN2A del`=="homdel"&IDH.codel.subtype=="Astro.", "G4", revised_grade)) 

# CIBERSORTx estimated cellular fractions based on input CARE-MUT signatures. We found these to perform fairly well.
tcga_fractions <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/public/tcga_lgggbm_cibersortx_celltype_state_fractions.txt", header = TRUE, sep = "\t")
tcga_fractions$aliquot_barcode <- gsub("\\.", "-", tcga_fractions$Mixture)
tcga_fractions$Sample.ID <- substr(tcga_fractions$aliquot_barcode, 1, 15)

tcga_fractions_subtype <- tcga_fractions %>% 
  inner_join(tcga_clinical_idhmut, by="Sample.ID") %>% 
  dplyr::select(-c(Mixture, P.value, Correlation, RMSE)) %>% 
  pivot_longer(cols= c(AC.like:Lymphocyte),
               names_to = "cell_type",
               values_to = "freq") 

tcga_fractions_subtype$IDH.codel.subtype <- factor(tcga_fractions_subtype$IDH.codel.subtype, levels=c("Oligo.", "Astro."))

pdf(paste0(fig_dir, "tcga_celltype_tumor_type_comparisons.pdf"), width = 2.75, height = 3, useDingbats = FALSE, bg="transparent")
ggplot(tcga_fractions_subtype %>% 
         filter(cell_type%in%c("MuralEndothelial","Oligodendrocyte", "Myeloid", "Neuron")), aes(x = IDH.codel.subtype, y = freq*100)) + 
  geom_boxplot(aes(fill=cell_type), outlier.shape = NA) +
  geom_point(size = 0.5) +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", size = 2.25, label = "p.format") +
  labs(x="Tumor type", y="TCGA CIBERSORTx\ncell abundance (%)") +
  scale_fill_manual(values=c("Malignant" = "#FB8072", "Oligodendrocyte" = "#B3DE69", "Myeloid" = "#80B1D3", 
                             "Astrocyte" = "#BFBADA", "Neuron" = "#BC80BD", 
                             "InhNeuron" = "#FFED6F", "MuralEndothelial" = "#FCCDE5", "MuralEndothelial" = "#FFFFB3", "Lymphocyte" = "#8DD3C7")) +
  theme(strip.background = element_blank()) +
  facet_grid(.~cell_type, scales="free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_n_text(size = 2.25) + 
  plot_theme
dev.off()

# Extract exact p-value for Myeloid cells
tcga_myeloid <- tcga_fractions_subtype %>% 
  filter(cell_type%in%c("Myeloid"))
wilcox.test(tcga_myeloid$freq~tcga_myeloid$IDH.codel.subtype)$p.value # 1.47e-39

### ### ### ### ### ### ### ### ### ###
### Malignant relative proportions  ###
### ### ### ### ### ### ### ### ### ###
# Determine the malignant compartment fraction
tcga_fraction_malignant <- tcga_fractions %>% 
  inner_join(tcga_clinical_idhmut, by="Sample.ID") %>% 
  mutate(Malignant = NPC.like+OPC.like+Undifferentiated+AC.like+MES.like) %>% 
  dplyr::select(Sample.ID, Malignant)

tcga_fractions_subtype_adj <- tcga_fractions %>%  
  inner_join(tcga_clinical_idhmut, by="Sample.ID") %>% 
  dplyr::select(-c(Mixture, P.value, Correlation, RMSE)) %>% 
  pivot_longer(cols= c(AC.like:Lymphocyte),
               names_to = "cell_type",
               values_to = "freq") %>% 
  inner_join(tcga_fraction_malignant, by="Sample.ID") %>% 
  filter(cell_type%in%c("NPC.like","OPC.like", "Undifferentiated", "AC.like", "MES.like")) %>% 
  mutate(malignant_freq = freq/Malignant) %>% 
  mutate(cell_type = gsub("\\.", "-", cell_type))


tcga_fractions_subtype_adj$IDH.codel.subtype <- factor(tcga_fractions_subtype_adj$IDH.codel.subtype, levels=c("Oligo.", "Astro."))

pdf(paste0(fig_dir, "tcga_types_malignant_states.pdf"), width = 3.75, height = 2.75, useDingbats = FALSE)
ggplot(tcga_fractions_subtype_adj, aes(x = IDH.codel.subtype, y = malignant_freq*100)) + 
  geom_boxplot(aes(fill=cell_type), outlier.shape = NA) +
  geom_point(size = 0.6) +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", size = 2.25, label="p.format") +
  labs(x="Tumor type", y="TCGA CIBERSORTx relative\nmalignant cell proportion (%)") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  theme(strip.background = element_blank()) +
  facet_grid(.~cell_type, scales="free") +
  stat_n_text(size = 2.25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Extract exact p-value for Myeloid cells
tcga_opc <- tcga_fractions_subtype_adj %>% 
  filter(cell_type%in%c("OPC-like"))
wilcox.test(tcga_opc$freq~tcga_opc$IDH.codel.subtype)$p.value # 6.603952e-18

pdf(paste0(fig_dir, "tcga_tumor_types_molecular_grade_kruskal.pdf"), width = 4, height = 2.75, useDingbats = FALSE)
ggplot(tcga_fractions_subtype_adj %>% 
         filter(!is.na(revised_grade)), aes(x = revised_grade, y = malignant_freq*100)) + 
  geom_boxplot(aes(fill=cell_type), outlier.shape = NA) +
  geom_point(size = 0.6) +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(method = "kruskal", size = 2.25, label="p.format") +
  labs(x="Tumor grade", y="TCGA CIBERSORTx relative\nmalignant cell proportion (%)") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  theme(strip.background = element_blank()) +
  facet_grid(IDH.codel.subtype~cell_type, scales="free") +
  stat_n_text(size = 2.25)
dev.off()

pdf(paste0(fig_dir, "tcga_tumor_types_malignant_states_molecular_grade_wilcox.pdf"), width = 4, height = 2.75, useDingbats = FALSE)
ggplot(tcga_fractions_subtype_adj %>% 
         filter(!is.na(revised_grade)), aes(x = revised_grade, y = malignant_freq*100)) + 
  geom_boxplot(aes(fill=cell_type), outlier.shape = NA) +
  geom_point(size = 0.6) +
  plot_theme +
  theme(legend.position = "none") +
  stat_compare_means(method = "wilcox", size = 2.25, label="p.format") +
  labs(x="Tumor grade", y="TCGA CIBERSORTx relative\nmalignant cell proportion (%)") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  theme(strip.background = element_blank()) +
  facet_grid(IDH.codel.subtype~cell_type, scales="free") +
  stat_n_text(size = 2.25)
dev.off()