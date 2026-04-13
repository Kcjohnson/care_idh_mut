##################################
# Plot the proportion of malignant cells cycling per malignant state
# Author: Kevin Johnson
# Date Updated: 2026.04.10
##################################

library(tidyverse)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
setwd(proj_dir)
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")


# Read in the cell state classification based on p-value assignment.
care_state_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/scoring/caremut_malignant_cell_state_assignment.txt", sep="\t", header = TRUE)

# Relabel some of the features.
care_state_md <- care_state_md %>% 
  mutate(tumor_type = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."),
         isCC = recode(as.character(isCC), `TRUE` = "Cycling",
                       `FALSE` = "Non-cycling"),
         cell_state = recode(State, `MP_AC1_MUT` = "AC-like",
                             `MP_OPC_MUT` = "OPC-like",
                             `MP_NPC_MUT` = "NPC-like",
                             `MP_MES_MUT` = "MES-like",
                             `MP_AC2_MUT` = "AC-like",
                             "Undifferentiated" = "Undiff.")) 

care_state_freq_cc <- care_state_md %>% 
  group_by(tumor_type, cell_state, isCC) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts/sum(counts)) %>% 
  ungroup()

care_state_freq_cc$tumor_type <- factor(care_state_freq_cc$tumor_type, levels=c("Oligo.", "Astro."))
care_state_freq_cc$cell_state <- factor(care_state_freq_cc$cell_state, levels=c("Undiff.", "NPC-like", "OPC-like", "MES-like", "AC-like"))
care_state_freq_cc$isCC <- factor(care_state_freq_cc$isCC, levels=c("Non-cycling", "Cycling"))


pdf(paste0(fig_dir, "cycling_percentage_by_malignant_state.pdf"), width=2.75, height=2.25, bg = "transparent", useDingbats = FALSE)
ggplot(care_state_freq_cc %>% 
         filter(isCC == "Cycling"), aes(x = cell_state, y = freq * 100, fill = isCC)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = paste0(round(freq * 100, 1), "%")), vjust = -0.5, size = 2) +  # Add labels above bars
  scale_fill_manual(values = c("#6BAED6"),
                    labels = c("Non-cycling", "Cycling")) +
  labs(fill = "% Cycling") + 
  xlab("") + 
  plot_theme + 
  facet_wrap(.~tumor_type) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(y = "Cycling cell abundance (%)") +
  guides(fill = FALSE)
dev.off()


### END ###