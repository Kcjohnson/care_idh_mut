##################################
# Visualize MES-like abundance change in monoculture, co-culture experiments following irradiation
# Author: Kevin Johnson
##################################

library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(EnvStats)

fig_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/figures/"
out_data_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/perturbation/coculture/"
setwd("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/")
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")

# Malignant metaprogram assignment to human MGG152 cells that were either monoculture/co-culture and irradiated/control
mdata_out <- read.delim(file = paste0(out_data_dir, "mgg152_coculture_caremut_mp_select_state_assignment.txt"), sep = "\t", header = TRUE)
linker_info <- mdata_out %>% 
  dplyr::select(SampleID, species, exp_group, condition, technical_batch) %>% 
  distinct()
rownames(linker_info) <- NULL

malignant_pval_freq <- mdata_out %>% 
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

malignant_pval_freq$State <- factor(malignant_pval_freq$State, levels = rev(c("MES-like", "AC-like", "OPC-like", "NPC-like", "Undifferentiated")))

malignant_pval_freq <- malignant_pval_freq %>% 
  mutate(plot_id = case_when(exp_group == "malignant_only_control" ~ "Malignant\nmonoculture",
                             exp_group == "malignant+macrophage_control" ~ "Mal.+Mac.\nco-culture",
                             exp_group == "malignant_only_irradiation" ~ "Mal. monoculture\nirradiation",
                             exp_group == "malignant+macrophage_irradiation" ~ "Mal.+Mac. co-culture\nirradiation",
                             TRUE ~ NA_character_
  ))
malignant_pval_freq$plot_id <- factor(malignant_pval_freq$plot_id, levels=c("Malignant\nmonoculture", "Mal.+Mac.\nco-culture", "Mal. monoculture\nirradiation", "Mal.+Mac. co-culture\nirradiation"))

mono_vs_co_plot <- malignant_pval_freq %>% 
  filter(State=="MES-like", exp_group%in%c("malignant_only_control", "malignant+macrophage_control")) %>% 
  ggplot(aes(x=plot_id, y=freq*100)) +
  geom_boxplot(aes(fill = factor(State))) +
  geom_point(size = 0.5) +
  stat_compare_means(method="wilcox", label="p.format", size = 2.25) +
  labs(y = "MES-like cell abundance (%) - no irradiation ctrl", x = "Condition", fill="Cell State") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(fill=FALSE) +
  stat_n_text(size = 2.25)


# Define a batch controlled pairwise difference in MES-like abundance.
mes_abundance_diff <-  malignant_pval_freq %>% 
  filter(State=="MES-like")  %>%
  select(State, freq, species, condition, technical_batch) %>%
  pivot_wider(names_from = condition, values_from = freq) %>%
  mutate(mes_irradiated_change = irradiation - control,
         State = "MES-like") %>% 
  mutate(plot_id = case_when(species == "malignant_only" ~ "Malignant\nmonoculture",
                             species == "malignant+macrophage" ~ "Mal.+Mac.\nco-culture",
                             TRUE ~ NA_character_
  ))
mes_abundance_diff$plot_id <- factor(mes_abundance_diff$plot_id, levels=c("Malignant\nmonoculture", "Mal.+Mac.\nco-culture"))

# We reach statistical significance if we examine the pairwise difference.
irradiated_control_plot <- ggplot(mes_abundance_diff, aes(x=plot_id, y=mes_irradiated_change*100)) +
  geom_boxplot(aes(fill = factor(State))) +
  scale_fill_manual(values=c("MES-like"="#F77D58")) +
  geom_point(size = 0.5) +
  stat_compare_means(method = "wilcox", label="p.format", size = 2.25) + 
  plot_theme +
  labs(x = "Condition", y = "MES-like change (%) after irradiation") +
  guides(fill = FALSE)  +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  stat_n_text(size = 2.25)



pdf(paste0(fig_dir, "fig5e_monoculture_coculture_irradiation.pdf"), width = 2.5, height = 2.75, bg = "transparent")
plot_grid(mono_vs_co_plot, irradiated_control_plot, ncol = 2)
dev.off()

### END ###