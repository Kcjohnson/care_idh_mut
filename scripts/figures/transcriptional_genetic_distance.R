##################################
# Plot the difference in genetic distance (private mutations) versus pseudobulk state-controlled transcriptional distance
# Author: Kevin Johnson
##################################

library(tidyverse)
library(ggpubr)
library(EnvStats)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/"
fig_dir     <- file.path(proj_dir, "figures/")
setwd(proj_dir)

source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")

# Mutation distance represent 1-fraction_shared_mutation. In the event that there were multiple tumor pairs (e.g., WGS and WXS data), averages were taken.
mutation_distance <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/genetics/mutation_genetic_distance.tsv", header = TRUE, sep = "\t")
# State controlled Euclidean distance in pseudobulk gene expression
transcriptional_distance <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/rna/pseudobulk_transcriptional_euclidean_distance.txt", header = TRUE, sep = "\t")
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt", sep="\t", header = TRUE)

genomic_distance <- mutation_distance %>% 
  inner_join(transcriptional_distance, by="patient_id") %>% 
  inner_join(patient_md, by="patient_id") %>% 
  mutate(HM = ifelse(is.na(acquired_hypermutation_t1t2), "No", "Yes"))

pdf(file = paste0(fig_dir, "edf8e_scaled_transcriptomic_dist_mutation_genetic_dist.pdf"), height = 2.5, width = 3, useDingbats = FALSE)
ggplot(genomic_distance, aes(x=scaled_snv_distance, y=scaled_transcriptional_dist)) +
  geom_point(aes(color=HM), size = 0.8) +
  scale_color_manual(values = c("Yes" = "red",
                                "No"="black")) +
  labs(x = "Scaled genetic distance (I-R Mutations)", 
       y="Scaled transcriptional distance\n(I-R Euclidean distance)",
       color="Hypermutation") +
  stat_cor(method="spearman", cor.coef.name = "rho", size = 2.25) +
  plot_theme +
  theme(legend.position = "bottom")
dev.off()
