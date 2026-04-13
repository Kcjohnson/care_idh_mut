##################################
# Plot malignant states across IDH-mutant hierarchy previously described by Tirosh et al.
# Author: Kevin Johnson
# Date Updated: 2026.04.10
##################################

library(tidyverse)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(EnvStats)
library(scales)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
out_data_dir <- file.path(proj_dir, "processed_data/rna/")
script_dir  <- file.path(proj_dir, "scripts")

setwd(proj_dir)

# Helper scripts
source(file.path(script_dir, "utils", "plot_theme.R"))
source(file.path(script_dir, "utils", "caremut_utils.R"))

# IDH-mutant hierarchy scores
mp_hierarchy_scores <- read.delim(file = paste0(proj_dir, "/results/scoring/malignant_malignant_hierarchy_signature_scores.txt"), sep = "\t", row.names = 1, header = TRUE)

# Read in the cell state classification based on p-value assignment.
care_state_md <- read.delim(file = paste0(proj_dir, "/results/scoring/caremut_malignant_cell_state_assignment.txt"), sep="\t", header = TRUE)

# Relabel some of the features.
care_state_md <- care_state_md %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."),
         patient_id = sapply(strsplit(care_id, "T"), "[[", 1),
         timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2)), 
         cell_state = recode(State, `MP_AC1_MUT` = "AC-like",
                             `MP_OPC_MUT` = "OPC-like",
                             `MP_NPC_MUT` = "NPC-like",
                             `MP_MES_MUT` = "MES-like",
                             `MP_AC2_MUT` = "AC-like",
                             "Undifferentiated" = "Undiff.")) 

mp_hierarchy_scores_annot <- mp_hierarchy_scores %>% 
  inner_join(care_state_md, by=c("SampleID", "CellID"))

# Inspect how the malignant states are distributed across the previously defined malignant state hierarchy
mp_hierarchy_scores_annot_trim <- mp_hierarchy_scores_annot %>% 
  dplyr::select(SampleID, idh_codel_subtype) %>% 
  distinct()

state_hierarchy_median_scores <- mp_hierarchy_scores_annot %>% 
  group_by(SampleID, cell_state) %>% 
  summarise(median_lineage = median(LineagePlot), median_stemness = median(Stemness),
            counts = n())  %>% 
  ungroup() %>% 
  left_join(mp_hierarchy_scores_annot_trim, by="SampleID") %>% 
  filter(counts > 20)


state_hierarchy_median_scores$cell_state <- factor(state_hierarchy_median_scores$cell_state, levels=c("NPC-like", "Undiff.", "OPC-like", "MES-like", "AC-like"))
state_hierarchy_median_scores$idh_codel_subtype <- factor(state_hierarchy_median_scores$idh_codel_subtype, levels=c("Oligo.", "Astro."))

my_comparisons <- list(c("AC-like", "MES-like"),c("AC-like", "OPC-like"), c("OPC-like", "MES-like"), c("OPC-like", "Undiff."),
                       c("OPC-like", "NPC-like"), c("NPC-like", "Undiff."))

pdf(paste0(fig_dir, "hierarchy_median_stemness_score_comparisons.pdf"), width = 5, height = 5, bg = "transparent", useDingbats = FALSE)
ggplot(state_hierarchy_median_scores, aes(x=cell_state, y=median_stemness)) +
  geom_boxplot(outlier.shape = NA, aes(fill=cell_state)) +
  geom_point() +
  stat_compare_means(comparisons=my_comparisons, method="wilcox", label = "p.format", size=2) +
  plot_theme +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "OPC-like"="#E8F5A3",
                             "NPC-like"="#6BAED6",
                             "Undiff." = "gray90")) +
  facet_grid(idh_codel_subtype~., scales = "fixed", space = "free") +
  labs(y="Median stemness score\nper state per sample (min. 20 cells)", x="State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  guides(fill=FALSE) + 
  stat_n_text()
dev.off()


mp_hierarchy_scores_annot$cell_state <- factor(mp_hierarchy_scores_annot$cell_state, levels=c("NPC-like", "Undiff.", "OPC-like", "MES-like", "AC-like"))
mp_hierarchy_scores_annot$idh_codel_subtype <- factor(mp_hierarchy_scores_annot$idh_codel_subtype, levels=c("Oligo.", "Astro."))

# Adding some background guides that give context for the IDH-mutant malignant state hierarchy
# A resizing factor
k <- 0.9

hierarchy_plot <- ggplot(mp_hierarchy_scores_annot, aes(x = LineagePlot, y = Stemness, color=idh_codel_subtype)) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 2.25 * k), color = "gray90", size = 3, lineend = "round", linejoin = "round") +  # Vertical line
  geom_segment(aes(x = 0, xend = 3 * cos(pi / 4) * k, y = 0, yend = -3.5 * sin(pi / 4) * k), color = "gray90", size = 3, lineend = "round", linejoin = "round") +  # Diagonal line 1
  geom_segment(aes(x = 0, xend = -3 * cos(pi / 4) * k, y = 0, yend = -3.5 * sin(pi / 4) * k), color = "gray90", size = 3, lineend = "round", linejoin = "round") +  # Diagonal line 2
  geom_density_2d() +
  scale_color_manual(values = c("Astro." = "#800074", "Oligo." = "#298C8C")) +
  plot_theme +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "\nLineage Score", y = "\nStemness Score") +
  xlim(-2.25, 2.25) +
  ylim(-2.25, 2.25) +
  facet_grid(idh_codel_subtype ~ cell_state, scales = "fixed", space = "free") +
  guides(color = FALSE) 

pdf(paste0(fig_dir, "hierarchy_guides_and_density.pdf"),  width = 4.17, height = 3.5, bg = "transparent", useDingbats = FALSE)
hierarchy_plot
dev.off()

png(paste0(fig_dir, "hierarchy_guides_and_density.png"), width = 4.17, height = 3.5, units = "in", res = 300)
hierarchy_plot
dev.off()


# Due to how the PDFs are loaded with the guides, I need to export the background guid as a png and the density as a PDF.
# Create a shared theme for these two images so that the png can be placed behind the pdf
shared_theme <- list(
  plot_theme,
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    strip.text = element_text(size = 7)),
  xlim(-2.25, 2.25),
  ylim(-2.25, 2.25),
  facet_grid(idh_codel_subtype ~ cell_state, scales = "fixed", space = "fixed"),
  labs(x = "Lineage score", y = "Stemness score"),
  guides(color = "none")
)

# Then update both export calls:
png(paste0(fig_dir, "hierarchy_guides.png"), 
    width = 4, height = 3.25,units = "in", res = 300)
ggplot(mp_hierarchy_scores_annot, aes(x = LineagePlot, y = Stemness)) +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = 2.25 * k), 
               color = "gray90", size = 3, lineend = "round") +
  geom_segment(aes(x = 0, xend = 3 * cos(pi/4) * k, y = 0, yend = -3.5 * sin(pi/4) * k), 
               color = "gray90", size = 3, lineend = "round") +
  geom_segment(aes(x = 0, xend = -3 * cos(pi/4) * k, y = 0, yend = -3.5 * sin(pi/4) * k), 
               color = "gray90", size = 3, lineend = "round") +
  shared_theme
dev.off()


pdf(paste0(fig_dir, "hierarchy_density.pdf"), 
    width = 4, height = 3.25, bg = "transparent", useDingbats = FALSE)
ggplot(mp_hierarchy_scores_annot, aes(x = LineagePlot, y = Stemness, color = idh_codel_subtype)) +
  geom_density_2d() +
  scale_color_manual(values = c("Astro." = "#800074", "Oligo." = "#298C8C")) +
  shared_theme
dev.off()



### ### ### ### ### ### ###
### Hierarchy cycling cell proximity 
### ### ### ### ### ### ###
# Oligodendrogliomas
set.seed(123)
mp_hierarchy_assigned_oligo <- mp_hierarchy_scores_annot %>% 
  filter(idh_codel_subtype=="Oligo.") %>% 
  sample_n(10000, replace = FALSE) 

# Function to calculate the fraction of "isCC" within 0.3 distance (selected based on Tirosh) from the input coordinate
get_fraction_within_distance <- function(x_coord, y_coord, df) {
  distance <- sqrt((df$LineagePlot - x_coord)^2 + (df$Stemness - y_coord)^2)
  points_within_distance <- df[distance <= 0.3, ]
  fraction_within_distance_and_CC <- mean(points_within_distance$isCC)
  return(fraction_within_distance_and_CC)
}

fractions_within_distance_oligo <- mapply(
  get_fraction_within_distance,
  mp_hierarchy_assigned_oligo$LineagePlot,
  mp_hierarchy_assigned_oligo$Stemness,
  MoreArgs = list(df = mp_hierarchy_assigned_oligo)
)

# Add the fractions to the original data frame
mp_hierarchy_assigned_oligo$fraction_within_0.3_distance <- fractions_within_distance_oligo

oligo_cycling_plot <- ggplot(mp_hierarchy_assigned_oligo %>% 
                               # A few outliers
                               filter(fraction_within_0.3_distance<0.4), aes(x = LineagePlot, y = Stemness, color=fraction_within_0.3_distance*100)) +
  geom_point() +
  plot_theme +
  scale_color_viridis(name = "Neighboring cycling\ncells %") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 7)) +
  labs(x="Lineage score", y="Stemness score") 

pdf(paste0(fig_dir, "oligo_neighboring_cycling_cells.pdf"), width = 4.5, height = 4.25, bg = "transparent", useDingbats = FALSE)
oligo_cycling_plot
dev.off()

png(paste0(fig_dir, "oligo_neighboring_cycling_cells.png"), width = 4.5, height = 4.25, res = 300, units = "in")
oligo_cycling_plot
dev.off()


# Repeat for astrocytomas
set.seed(21)
mp_hierarchy_assigned_astro <- mp_hierarchy_scores_annot %>% 
  filter(idh_codel_subtype=="Astro.") %>% 
  sample_n(10000, replace = FALSE) 

fractions_within_distance_astro <- mapply(
  get_fraction_within_distance,
  mp_hierarchy_assigned_astro$LineagePlot,
  mp_hierarchy_assigned_astro$Stemness,
  MoreArgs = list(df = mp_hierarchy_assigned_astro)
)

# Add the fractions to the original data frame
mp_hierarchy_assigned_astro$fraction_within_0.3_distance <- fractions_within_distance_astro

ggplot(mp_hierarchy_assigned_astro, aes(x = LineagePlot, y = Stemness, color=fraction_within_0.3_distance*100)) +
  geom_point() +
  plot_theme +
  scale_color_viridis(name = "Neighboring cycling\ncells %") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Lineage score", y="Stemness score") 

hp_astro <- ggplot(mp_hierarchy_assigned_astro %>% 
                        # A few outliers
                        filter(fraction_within_0.3_distance<0.4), aes(x = LineagePlot, y = Stemness, color=fraction_within_0.3_distance*100)) +
  geom_point() +
  plot_theme +
  scale_color_viridis(name = "Neighboring cycling\ncells %") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Lineage score", y="Stemness score") 

ggsave(paste0(fig_dir, "astro_neighboring_cycling_cells.pdf"), hp_astro,  width = 4.5, height = 4.25, dpi = 300)
ggsave(paste0(fig_dir, "astro_neighboring_cycling_cells.png"), hp_astro,  width = 4.5, height = 4.25, dpi = 300)

set.seed(42)
mp_hierarchy_assigned_all_ds <- mp_hierarchy_assigned_astro %>% 
  bind_rows(mp_hierarchy_assigned_oligo) %>% 
  group_by(idh_codel_subtype) %>% 
  sample_n(5000, replace = FALSE) 


mp_hierarchy_assigned_all_ds$idh_codel_subtype <- factor(mp_hierarchy_assigned_all_ds$idh_codel_subtype, levels=c("Oligo.", "Astro."))

hp_neighboring_cc <- ggplot(mp_hierarchy_assigned_all_ds, aes(x = LineagePlot, y = Stemness, color=fraction_within_0.3_distance*100)) +
  geom_point() +
  plot_theme +
  scale_color_viridis(name = "Neighboring cycling\ncells %", limits = c(0, 30), oob = squish) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "bottom") +
  labs(x="Lineage score", y="Stemness score") +
  facet_grid(idh_codel_subtype~., space = "free") 


