##############################
## Visualize experimental PDGFRAi analyses
## Author: Kevin Johnson
##############################

# Listed as 637451 in the data but should be 673451

library(tidyverse)
library(ggpubr)
library(EnvStats)
library(openxlsx)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "results/figures/perturbation/")
setwd(proj_dir)

source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")

# Experimental data that accompanied the PDGFRAi experiments
# Trypan-based viability
viability_df <- readWorkbook("data/perturbation/pdgfrai/norlux_cell_viability_pdgfrai.xlsx", sheet = 1, startRow = 1)

# FACS cell cycle phase
cellcycle_df <- readWorkbook("data/perturbation/pdgfrai/norlux_cellcycle_phase_pdgfrai.xlsx", sheet = 1, startRow = 1)

# Reformatting to extract replicates, treatment conditions, and SampleID
viability_df$replicate <- sapply(strsplit(viability_df$exp_replicate, " R"), "[[", 2)
viability_df$treatment <- sapply(strsplit(viability_df$exp_replicate, " R"), "[[", 1)
viability_df$treatment_group <- sapply(strsplit(viability_df$exp_replicate, " "), "[[", 1)
viability_df <- viability_df %>% 
  mutate(SampleID = paste0(cell_line, " ", treatment))
viability_df$SampleID <- factor(viability_df$SampleID, levels=c("T394NS DMSO ctr", "T394NS Dasatinib 1uM", "T394NS Dasatinib 5uM", "T394NS CP-637451 1uM", "T394NS CP-637451 5uM",
                                                                "T407NS DMSO ctr", "T407NS Dasatinib 1uM", "T407NS Dasatinib 5uM", "T407NS CP-637451 1uM", "T407NS CP-637451 5uM"))

viability_df$treatment <- factor(viability_df$treatment, levels=c("DMSO ctr", "Dasatinib 1uM", "Dasatinib 5uM", "CP-637451 1uM", "CP-637451 5uM"))

pdf(paste0(fig_dir, "pdgfrai_trypan_blue_relative_viability.pdf"), width = 2.75, height = 2.75, useDingbats = FALSE)
ggplot(viability_df, aes(x = treatment, y=relative2dmso_live)) +
  geom_boxplot(aes(fill = treatment_group), outlier.shape = NA) +
  geom_point(size = 0.5) + 
  plot_theme +
  facet_grid(.~cell_line, space = "free", scales = "free") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = c(
    "DMSO" = "#00BFC4",
    "Dasatinib" = "#9590FF",
    "CP-637451" = "#FF62BC"
  )) +
  labs(y = "% relative viable cell count",
       x = "Treatment") +
  ylim(0, 105) +  
  stat_compare_means(method="kruskal", label="p.format", size = 2.25) +
  stat_n_text(size = 2.25)
dev.off()

# Reformat the cell cycle analyses
cellcycle_df$replicate <- sapply(strsplit(cellcycle_df$exp_replicate, " R"), "[[", 2)
cellcycle_df$treatment <- sapply(strsplit(cellcycle_df$exp_replicate, " R"), "[[", 1)

cellcycle_df_long <- cellcycle_df %>% 
  filter(!is.na(G1|S|G2)) %>% 
  pivot_longer(cols = c("G1", "S", "G2"), names_to = "phase", values_to = "percent") %>% 
  mutate(SampleID = paste0(cell_line, " ", treatment))

# The FACS fitting algorithm/software does not assign a small proportion of events to any phase. 
# Rescaling G0/G1 + S + G2/M are to a total of 100% after these events have been excluded.
cellcycle_df_long <- cellcycle_df_long %>%
  group_by(SampleID, exp_replicate) %>%
  mutate(rescaled_percent = (percent / sum(percent)) * 100) %>%
  ungroup()


cellcycle_df_long$SampleID <- factor(cellcycle_df_long$SampleID, levels=c("T394NS DMSO ctr", "T394NS Dasatinib 1uM", "T394NS Dasatinib 5uM", "T394NS CP-637451 1uM", "T394NS CP-637451 5uM",
                                                                          "T407NS DMSO ctr", "T407NS Dasatinib 1uM", "T407NS Dasatinib 5uM", "T407NS CP-637451 1uM", "T407NS CP-637451 5uM"))

cellcycle_df_long$treatment <- factor(cellcycle_df_long$treatment, levels = c("DMSO ctr", "Dasatinib 1uM", "Dasatinib 5uM", "CP-637451 1uM", "CP-637451 5uM"))

pdf(paste0(fig_dir, "pdgfrai_facs_cell_cycle_phases_boxplot.pdf"), width = 2.25, height = 3, useDingbats = FALSE)
ggplot(cellcycle_df_long, aes(x = treatment, y = percent, fill = phase)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.6) + 
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c(
    "G1" = "#3b4fe5",
    "S"  = "#c9c940",
    "G2" = "#2e6f1e"
  )) +
  facet_grid(phase~cell_line) +
  labs(x="Treatment", y="Cell cycle phase (%)")
dev.off()

# Visually inspect the breakdown for these two replicates' averages.
cellcycle_summary <- cellcycle_df_long %>% 
  group_by(treatment, phase, cell_line) %>%
  summarise(
    mean = mean(rescaled_percent),
    sem  = sd(rescaled_percent) / sqrt(n()),
    .groups = "drop"
  )


cellcycle_summary$phase <- factor(cellcycle_summary$phase, levels = rev(c("G1", "S", "G2")))
cellcycle_summary$treatment <- factor(cellcycle_summary$treatment, levels = c("DMSO ctr", "Dasatinib 1uM", "Dasatinib 5uM", "CP-637451 1uM", "CP-637451 5uM"))


ggplot(cellcycle_summary,
       aes(x = treatment, y = mean, fill = phase)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c(
    "G1" = "#3b4fe5",
    "S"  = "#c9c940",
    "G2" = "#2e6f1e"
  )) +
  labs(y = "Mean cell cycle phase (%)",
       x = "Treatment",
       fill = "Phase") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  facet_grid(.~cell_line) +
  geom_hline(yintercept = 100, linetype = "dashed")


### END ###