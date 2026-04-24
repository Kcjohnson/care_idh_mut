##################################
# Plot malignant metaprogram scores across a heatmap stratified by assigned cell state
# Author: Kevin Johnson
# Date Updated: 2026.04.10
##################################

library(tidyverse)
library(RColorBrewer)
library(viridis)
library(openxlsx)
library(scales)
library(reshape2)
library(ggpubr)
library(patchwork)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
out_data_dir <- file.path(proj_dir, "processed_data/rna/")
script_dir  <- file.path(proj_dir, "scripts")

setwd(proj_dir)

# Helper scripts
source(file.path(script_dir, "utils", "plot_theme.R"))

# Read in the cell state classification based on p-value assignment.
care_state_md <- read.delim(file = paste0(proj_dir, "/results/scoring/caremut_malignant_cell_state_assignment.txt"), sep="\t", header = TRUE)

# Relabel some of the features.
care_state_md <- care_state_md %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."),
         patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2)), 
         MalignantState = recode(State, `MP_AC1_MUT` = "AC-like",
                             `MP_OPC_MUT` = "OPC-like",
                             `MP_NPC_MUT` = "NPC-like",
                             `MP_MES_MUT` = "MES-like",
                             `MP_AC2_MUT` = "AC-like",
                             "Undifferentiated" = "Undiff.")) 

set.seed(34)
sampled_cells <- care_state_md %>% 
  sample_n(10000)

sampled_cells <- sampled_cells %>%
  dplyr::select(CellID, SampleID, patient_id, timepoint, starts_with("MP"), isCC, MalignantState)

# The following code excludes the MPs we don't want to show
sampled_cells <- sampled_cells %>%
  select(-"MP_AC2_MUT")

mp_order <- c("MP_NPC",
              "MP_OPC",
              "MP_AC1",
              "MP_AC2",
              "MP_MES",
              "MP_CC")

scores <- sampled_cells %>%
  dplyr::select(MP_AC = MP_AC1_MUT, MP_MES = MP_MES_MUT, MP_OPC = MP_OPC_MUT, MP_NPC = MP_NPC_MUT, MP_CC = MP_CC_MUT)
scores <- as.matrix(scores)
rownames(scores) <- sampled_cells$CellID
scores <- t(scores)

dm <- melt(scores) %>%
  as_tibble()
colnames(dm) <- c("MP", "CellID", "Score")

dm <- dm %>%
  left_join(sampled_cells %>%
              select(CellID, MalignantState), by = "CellID")

# Set the order of facets
dm$MalignantState <- factor(dm$MalignantState, c("NPC-like",
                                                 "OPC-like",
                                                 "Undiff.",
                                                 "MES-like",
                                                 "AC-like"))

# Set the order of rows (y-axis)
dm$MP <- gsub("MP_", "", dm$MP)
dm$MP <- factor(dm$MP, rev(c("NPC",
                             "OPC",
                             "MES",
                             "AC",
                             "CC")))

# I ordered the cells according to the order of the cell cycle score
dm$CellID <- factor(dm$CellID, colnames(scores)[order(scores["MP_CC", ])])


p1 <- dm %>%
  filter(MP != "CC") %>%
  ggplot(aes(x = CellID, y = MP, fill = Score)) +
  facet_grid(cols = vars(MalignantState), scales = "free", space = "free_x") +
  geom_tile() +
  scale_fill_gradient2(
    name = "MP Score",
    low = "dodgerblue", mid = "white", high = "red",
    labels = c("-1", "0", "1"),
    breaks = c(-1, 0, 1),
    limits = c(-1, 1), oob = squish
  ) +
  xlab("") +
  ylab("Malignant metaprogram (MP)") +
  guides(fill = guide_colourbar(
    ticks.colour = "black",
    frame.colour = "black",
    barwidth  = unit(0.1, "in"),
    barheight = unit(0.6, "in")
  )) +
  theme_pvsr(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    strip.text = element_text(size = 7, angle = 0, hjust = 0.5),
    legend.position = "right"
  ) +
  theme(
    strip.background = element_blank(),
    strip.clip = "off",
    plot.margin = margin(t = 10, r = 5, b = -5, l = 5, unit = "pt")
  )

p2 <- dm %>%
  filter(MP == "CC") %>%
  ggplot(aes(x = CellID, y = MP, fill = Score)) +
  facet_grid(cols = vars(MalignantState), scales = "free", space = "free_x") +
  geom_tile() +
  scale_fill_gradient2(
    name = "MP Score",
    low = "dodgerblue", mid = "white", high = "red",
    labels = c("-1", "0", "1"),
    breaks = c(-1, 0, 1),
    limits = c(-1, 1), oob = squish
  ) +
  xlab("Cells") +
  ylab("") +
  guides(fill = guide_colourbar(ticks.colour = "black", frame.colour = "black")) +
  theme_pvsr(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    strip.text = element_blank(),
    legend.position = "none"
  ) +
  theme(
    strip.background = element_blank(),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt")
  )

pdf(paste0(fig_dir, "malignant_state_heatmap.pdf"), width = 4.17, height = 2.6, bg = "transparent", useDingbats = FALSE)
p1 / p2 + plot_layout(heights = c(4, 1))
dev.off()


p1 <- dm %>%
  filter(MP != "CC") %>%
  ggplot(aes(x = CellID, y = MP, fill = Score)) +
  facet_grid(cols = vars(MalignantState), scales = "free", space = "free_x") +
  geom_tile() +
  scale_fill_gradient2(
    name = "MP Score",
    low = "dodgerblue", mid = "white", high = "red",
    labels = c("-1", "0", "1"),
    breaks = c(-1, 0, 1),
    limits = c(-1, 1), oob = squish
  ) +
  xlab("") +
  ylab("Malignant metaprogram (MP)") +
  guides(fill = "none") +  # hide legend in p1
  theme_pvsr(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 7),
    strip.text = element_text(size = 7, angle = 0, hjust = 0.5),
    legend.position = "none"
  ) +
  theme(
    strip.background = element_blank(),
    strip.clip = "off",
    plot.margin = margin(t = 10, r = 5, b = -5, l = 5, unit = "pt")
  )

p2 <- dm %>%
  filter(MP == "CC") %>%
  ggplot(aes(x = CellID, y = MP, fill = Score)) +
  facet_grid(cols = vars(MalignantState), scales = "free", space = "free_x") +
  geom_tile() +
  scale_fill_gradient2(
    name = "MP Score",
    low = "dodgerblue", mid = "white", high = "red",
    labels = c("-1", "0", "1"),
    breaks = c(-1, 0, 1),
    limits = c(-1, 1), oob = squish
  ) +
  xlab("Cells") +
  ylab("") +
  guides(fill = guide_colourbar(
    ticks.colour = "black",
    frame.colour = "black",
    barwidth  = unit(0.6, "in"),
    barheight = unit(0.08, "in"),
    title.position = "top",
    title.hjust = 0.5
  )) +
  theme_pvsr(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    strip.text = element_blank(),
    legend.position = "bottom"
  ) +
  theme(
    strip.background = element_blank(),
    plot.margin = margin(t = 0, r = 5, b = 5, l = 5, unit = "pt")
  )

pdf(paste0(fig_dir, "malignant_state_heatmap_legend.pdf"), width = 4.17,  height = 2.6, bg = "transparent", useDingbats = FALSE)
p1 / p2 + plot_layout(heights = c(4, 1))
dev.off()
