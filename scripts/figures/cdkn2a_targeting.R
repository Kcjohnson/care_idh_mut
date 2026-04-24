##################################
# Plot the assigned malignant state abundance across different conditions in CDKN2A-/- experiments
# Author: Kevin Johnson
##################################

library(tidyverse)
library(ggpubr)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
source(paste0(proj_dir, "/scripts/utils/plot_theme.R"))
setwd(proj_dir)


cdkn2a_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/perturbation/cdkn2a/cdkn2a_caremut_mp_select_state_assignment.txt", sep = "\t", row.names = 1, header = TRUE)

# Recode and collapse malignant states from AC1 and AC2.
malignant_pval_freq <- cdkn2a_md %>% 
  mutate(State = recode(State, `MP_OPC_MUT` = "OPC-like",
                        `MP_NPC_MUT` = "NPC-like",
                        `MP_AC1_MUT` = "AC-like", 
                        `MP_AC2_MUT` = "AC-like",
                        `MP_MES_MUT` = "MES-like",
                        `Undifferentiated` = "Undifferentiated")) %>% 
  group_by(exp_id, State) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(exp_id, State,
           fill = list(counts = 0, freq = 0)) 


# Independently calculate the cell cycle.
cdkn2a_state_freq_cc <- cdkn2a_md %>% 
  group_by(exp_id, isCC) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(exp_id, isCC,
           fill = list(counts = 0, freq = 0)) %>%
  filter(isCC=="TRUE") %>% 
  dplyr::select(-isCC) %>% 
  mutate(State="Cycling") 

malignant_pval_freq <- malignant_pval_freq %>% 
  bind_rows(cdkn2a_state_freq_cc)

malignant_pval_freq$State <- factor(malignant_pval_freq$State, levels = rev(c("MES-like", "AC-like", "OPC-like", "NPC-like", "Undifferentiated", "Cycling")))
malignant_pval_freq$exp_group <- ifelse(malignant_pval_freq$exp_id%in%c("sgNTC","Parental"), "Control", "sgRNA")

malignant_pval_freq_wide <- malignant_pval_freq %>%
  mutate(freq = freq*100) %>% 
  dplyr::select(exp_id, State, freq) %>% 
  pivot_wider(names_from = State,  values_from = freq) %>% 
  mutate(Undiff_Stem = Undifferentiated+`NPC-like`,
         AC_MES = `MES-like`+`AC-like`) %>% 
  mutate(exp_group = ifelse(exp_id%in%c("sgNTC","Parental"), "Control", "CDKN2A-/-"))


malignant_pval_freq_wide <- malignant_pval_freq_wide %>% 
  mutate(exp_type = recode(exp_id, `sgRNA1` = "CDKN2A-/-",
                           `sgRNA3` = "CDKN2A-/-",
                           `sgNTC` = "NTC",
                           `Parental` = "Parental")) 


parent <- malignant_pval_freq_wide[malignant_pval_freq_wide$exp_id == "Parental", ]

arrows_df <- malignant_pval_freq_wide %>%
  filter(exp_id != "Parental") %>%
  mutate(y = parent$Undiff_Stem,
         x = parent$AC_MES,
         yend = Undiff_Stem,
         xend = AC_MES)

pdf(file = paste0(fig_dir, "edf9_cdkn2a_undiffstem.pdf"), width = 2.5, height = 2.5,  useDingbats = FALSE)
ggplot(malignant_pval_freq_wide, aes(x = AC_MES, y = Undiff_Stem)) +
  geom_point(aes(fill = exp_group), shape = 21, size = 2, color = "black") +
  scale_fill_manual(values=c("Control" = "white", 
                             "CDKN2A-/-" = "#1F78B4")) +
  geom_segment(data = arrows_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 25),
               linetype = "dashed",
               size = 0.25,
               color = "black") +
  geom_text(
    aes(label = exp_type),
    nudge_x = -2.25,    
    nudge_y = 0,    
    size = 2.5
  ) +
  labs(x = "Astrocyte lineage (AC/MES-like %)", y = "Stem-like (Undiff./NPC-like %)") +
  plot_theme +
  #ylim(85, 93) +
  #xlim(0, 15) +
  guides(fill = FALSE)
dev.off()


