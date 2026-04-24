##################################
# Analyze malignant cell abundance association with common genetic events in IDH-mutant glioma TCGA dataset.
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


# Selected genetic events - IDH1, CDK4, PDGFRA, CDKN2A plus relevant genetic events for IDH-O drivers were also added
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
  mutate(IDH.codel.subtype = recode(IDH.codel.subtype, `IDHmut-non-codel` = "IDH-A",
                                    `IDHmut-codel` = "IDH-O")) %>% 
  left_join(cbio_genetics_idh, by=c("Sample.ID"="Sample.ID")) %>% 
  mutate(revised_grade = ifelse(Grade%in%c("G2", "G3") & `CDKN2A del`=="homdel"&IDH.codel.subtype=="IDH-A", "G4", Grade)) %>% 
  mutate(revised_grade = ifelse(is.na(revised_grade) & `CDKN2A del`=="homdel"&IDH.codel.subtype=="IDH-A", "G4", revised_grade)) 

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


tcga_oligo_long <- tcga_fractions_subtype_adj %>%
  filter(!is.na(`NOTCH1 mut`)) %>% 
  pivot_longer(cols = c("CDKN2A del", "NOTCH1 mut", "CIC mut", "FUBP1 mut", "PIK3CA mut", "PIK3R1 mut", "NF1 mut"),
               names_to = "alteration", values_to = "status") %>%
  mutate(has_alt = ifelse(status != "no alt." & !is.na(status), 1, 0))

tcga_oligo_results <- tcga_oligo_long %>%
  filter(IDH.codel.subtype=="IDH-O", cell_type=="AC-like") %>% 
  group_by(alteration) %>% 
  summarise(
    n_alt = sum(has_alt == 1, na.rm = TRUE),   # number of altered samples
    n_noalt = sum(has_alt == 0, na.rm = TRUE),
    diff = mean(malignant_freq[has_alt == 1], na.rm = TRUE) -
      mean(malignant_freq[has_alt == 0], na.rm = TRUE),
    pval = tryCatch(
      wilcox.test(malignant_freq ~ has_alt)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(padj = p.adjust(pval, method = "fdr"),
         log10padj = -log10(padj),
         log10pval = -log10(pval),
         alt_type = case_when(alteration == "CDK4 amp" ~ "amplification",
                              alteration == "PDGFRA amp" ~ "amplification",
                              alteration == "CDKN2A del" ~ "deletion",
                              grepl("mut", alteration) ~ "mutation"))

oligo_genetic_res <- ggplot(tcga_oligo_results, aes(x = diff*100, y = log10padj, color = alt_type)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_text_repel(
    aes(label = paste0(alteration, " (n=", n_alt, ")")),
    size = 2.25, show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  labs(
    x = "Mean % AC-like abundance diff.\n(alt vs no alt) TCGA oligo.",
    y = "-log10 adj. p-value",
    color = "Gene event"
  )  + 
  scale_color_manual(values =c("amplification" = "red",
                               "deletion" = "dodgerblue",
                               "mutation" = "springgreen4")) +
  plot_theme + 
  xlim(-25, 25) +
  ylim(0, 4)  +
  guides(color = FALSE)


## Astrocytoma
tcga_astro_long <- tcga_fractions_subtype_adj %>%
  filter(!is.na(`CDKN2A del`)) %>% 
  pivot_longer(cols = c("CDKN2A del", "PDGFRA amp", "CDK4 amp", "PIK3CA mut", "PIK3R1 mut", "NF1 mut", "NOTCH1 mut"),
               names_to = "alteration", values_to = "status") %>%
  mutate(has_alt = ifelse(status != "no alt." & !is.na(status), 1, 0))


tcga_astro_results <- tcga_astro_long %>%
  filter(IDH.codel.subtype=="IDH-A", cell_type=="AC-like") %>% 
  group_by(alteration) %>% 
  summarise(
    n_alt = sum(has_alt == 1, na.rm = TRUE),   # number of altered samples
    n_noalt = sum(has_alt == 0, na.rm = TRUE),
    diff = mean(malignant_freq[has_alt == 1], na.rm = TRUE) -
      mean(malignant_freq[has_alt == 0], na.rm = TRUE),
    pval = tryCatch(
      wilcox.test(malignant_freq ~ has_alt)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(padj = p.adjust(pval, method = "fdr"),
         log10padj = -log10(padj),
         log10pval = -log10(pval),
         alt_type = case_when(alteration == "CDK4 amp" ~ "amplification",
                              alteration == "PDGFRA amp" ~ "amplification",
                              alteration == "CDKN2A del" ~ "deletion",
                              grepl("mut", alteration) ~ "mutation"))

astrocytoma_genetic_results <- ggplot(tcga_astro_results, aes(x = diff*100, y = log10padj, color = alt_type)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_text_repel(
    aes(label = paste0(alteration, " (n=", n_alt, ")")),
    size = 2.25, show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  labs(
    x = "Mean % AC-like abundance diff.\n(alt vs no alt) TCGA astro.",
    y = "-log10 adj. p-value",
    color = "Genetic event"
  )  + 
  scale_color_manual(values =c("amplification" = "red",
                               "deletion" = "dodgerblue",
                               "mutation" = "springgreen4")) +
  plot_theme + 
  xlim(-25, 25) +
  ylim(0, 4)  +
  guides(color = FALSE)

pdf("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/figures/edf9_tcga_genetics_ac_like_abundance.pdf", width = 3.75, height = 2.75, useDingbats = FALSE, bg = "transparent")
plot_grid(oligo_genetic_res, astrocytoma_genetic_results, ncol = 2)
dev.off()

### END ###