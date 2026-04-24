##################################
# Plot the assigned malignant state abundance across different conditions in PDGFRAi experiments
# Author: Kevin Johnson
##################################

library(tidyverse)
library(ggpubr)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
source(paste0(proj_dir, "/scripts/utils/plot_theme.R"))
setwd(proj_dir)


# 10x cell state assignment
pdgfrai_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/perturbation/pdgfrai/pdgfrai_caremut_mp_select_state_assignment_20260406.txt", sep = "\t", row.names = 1, header = TRUE)

# Create an easy linker file to add back after calculating cell state frequencies.
linker_info <- pdgfrai_md %>% 
  dplyr::select(SampleID, exp_group:dose) %>% 
  distinct()

# Recode and collapse malignant states from AC1 and AC2.
malignant_pval_freq <- pdgfrai_md %>% 
  mutate(State = recode(State, `MP_OPC_MUT` = "OPC-like",
                        `MP_NPC_MUT` = "NPC-like",
                        `MP_AC1_MUT` = "AC-like", 
                        `MP_AC2_MUT` = "AC-like",
                        `MP_MES_MUT` = "MES-like",
                        `Undifferentiated` = "Undifferentiated")) %>% 
  group_by(SampleID, State) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(SampleID, State,
           fill = list(counts = 0, freq = 0)) %>%
  inner_join(linker_info, by="SampleID")

# Confirm that all sample sum 1 and that all cell types are measured.
malignant_pval_freq %>%
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
table(malignant_pval_freq$State)

# Independently calculate the cell cycle.
pdgfrai_state_freq_cc <- pdgfrai_md %>% 
  group_by(SampleID, isCC) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(SampleID, isCC,
           fill = list(counts = 0, freq = 0)) %>%
  filter(isCC=="TRUE") %>% 
  dplyr::select(-isCC) %>% 
  mutate(State = "Cycling") %>% 
  inner_join(linker_info, by="SampleID") 

malignant_pval_freq <- malignant_pval_freq %>% 
  bind_rows(pdgfrai_state_freq_cc)

malignant_pval_freq$State <- factor(malignant_pval_freq$State, levels = c("MES-like", "AC-like", "OPC-like", "NPC-like", "Undifferentiated", "Cycling"))

malignant_pval_freq %>% 
  # Omitted the cycling cells so that all frequencies sum to 100%. As a reminder the cycling cells are separately enumerated.
  filter(State!="Cycling") %>% 
  ggplot(aes(x=exp_group, fill = factor(State), y=freq*100)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "Relative cell abundance (%)", fill="Cell State", x = "") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  facet_grid(.~cell_line, scales = "free_x", space = "free_x") +
  plot_theme +
  theme(panel.grid.major = element_line(),
        axis.text.x = element_text(angle = 45, hjust=1)) 


# There is not a clear reduction in proliferation due to PDGFRAi, perhaps due to the FACS enrichment.
# We observed a viability reduction in experimental assays
pdgfrai_state_freq_cc %>% 
  ggplot(aes(x=exp_group, fill = factor(State), y=freq*100)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "Relative cell abundance (%)", fill="Is Cycling?", x = "") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6")) +
  facet_grid(.~cell_line, scales = "free_x", space = "free_x") +
  plot_theme +
  theme(panel.grid.major = element_line(),
        axis.text.x = element_text(angle = 45, hjust=1)) 

### ### ### ### ### ### ### ###
### 2D graph for directionality of changes
### ### ### ### ### ### ### ###
malignant_pval_freq_wide <- malignant_pval_freq %>%
  mutate(freq = freq*100) %>% 
  dplyr::select(SampleID, State, freq) %>% 
  pivot_wider(names_from = State,  values_from = freq) %>% 
  mutate(Undiff_Stem = Undifferentiated+`NPC-like`,
         AC_MES = `AC-like`+`MES-like`) %>% 
  inner_join(linker_info, by="SampleID") 


arrow_df <- malignant_pval_freq_wide %>%
  group_by(cell_line) %>%
  # Identify control row per cell line
  mutate(ctrl_x = AC_MES[which(treatment == "DMSO")],
         ctrl_y = Undiff_Stem[which(treatment == "DMSO")]) %>%
  ungroup() %>%
  filter(treatment != "DMSO") %>%
  mutate(x = ctrl_x, y = ctrl_y, xend = AC_MES, yend = Undiff_Stem)

pdf(file = paste0(fig_dir, "fig4_pdgfra_inhibition_acmes_stem.pdf"), width = 4, height = 3.75,  useDingbats = FALSE)
ggplot(malignant_pval_freq_wide %>% 
         mutate(text_to_add = paste(treatment, " ", dose)), aes(x = AC_MES, y = Undiff_Stem)) +
  geom_point(aes(color = treatment), size = 2) +
  geom_segment(data = arrow_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 25),
               linetype = "dashed",
               size = 0.25,
               color = "black") +
  geom_text(
    aes(label = text_to_add),
    nudge_x = 2,    
    nudge_y = 0.75,    
    size = 2.25
  ) +
  facet_grid(.~cell_line, space = "free") +
  scale_color_manual(values=c("DMSO" = "#00BFC4",
                              "Dasatinib" = "#9590FF",
                              "CP-673451" = "#FF62BC")) +
  labs(x = "Astrocyte lineage (AC/MES-like %)", y = "Stem-like (Undiff./NPC-like %)", shape = "Group", color = "Condition") +
  plot_theme +
  guides(color = FALSE)
dev.off()

arrow_df <- malignant_pval_freq_wide %>%
  group_by(cell_line) %>%
  # Identify control row per cell line
  mutate(ctrl_x = AC_MES[which(treatment == "DMSO")],
         ctrl_y = Undifferentiated[which(treatment == "DMSO")]) %>%
  ungroup() %>%
  filter(treatment != "DMSO") %>%
  mutate(x = ctrl_x, y = ctrl_y, xend = AC_MES, yend = Undifferentiated)

pdf(file = paste0(fig_dir, "fig4_pdgfra_inhibition_acmes_undifferentiated.pdf"), width = 5, height = 3,  useDingbats = FALSE)
ggplot(malignant_pval_freq_wide %>% 
         mutate(text_to_add = paste(treatment, " ", dose)), aes(x = AC_MES, y = Undifferentiated)) +
  geom_point(aes(color = treatment), size = 2) +
  geom_segment(data = arrow_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 25),
               linetype = "dashed",
               size = 0.25,
               color = "black") +
  geom_text(
    aes(label = text_to_add),
    nudge_x = 2,    
    nudge_y = 0.75,    
    size = 2.5
  ) +
  facet_wrap(~cell_line, scales = "fixed", nrow = 1) +
  scale_color_manual(values=c("DMSO" = "#00BFC4",
                              "Dasatinib" = "#9590FF",
                              "CP-673451" = "#FF62BC")) +
  labs(x = "Astrocyte lineage abundance (AC/MES-like %)", y = "Undifferentiated abundance %", shape = "Group", color = "Condition") +
  plot_theme +
  guides(color = FALSE)
dev.off()

