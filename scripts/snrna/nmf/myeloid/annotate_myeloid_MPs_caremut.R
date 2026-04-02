##################################
# Derive and annotate Myeloid MPs for CAREmut dataset
# Author: Kevin Johnson
# Date: 2026.03.30
##################################


proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures")
out_data_dir <- file.path(proj_dir, "processed_data")
script_dir  <- file.path(proj_dir, "scripts")

setwd(proj_dir)

# Load the NMF helper functions (implementation of the MP algorithm) 
source(file.path(script_dir, "utils", "PvsR-NMF-caremut.R"))

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Derive NMF MPs
#
# In this part we load the NMF results that were computed on HPC. We expect an NMF file for each sample/cell type and
# generate the MPs for each cell type independently
#
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


# Path to NMF results directory, e.g.
nmf_res_dir <- "/vast/palmer/pi/verhaak/kcj28/care_mut/results/nmf_res_2026/myeloid_n73/"

# List files in the directory
nmf_files <- list.files(nmf_res_dir)

# Generate the NMF scores matrix across all samples and candidate programs
Genes_nmf_w_basis <- lapply(nmf_files, function(x) {
  f <- readRDS(file = paste0(nmf_res_dir, x))
  
  res <- get_nmf_programs(f$fit, x)
  
  rm(f)
  gc()
  
  return(res)
})
names(Genes_nmf_w_basis) <- nmf_files
length(Genes_nmf_w_basis)

# Call the NMF algorithm
myeloid_nmf_metaprograms <- derive_NMF_metaprograms(Genes_nmf_w_basis = Genes_nmf_w_basis, 
                                                      save_path = "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/metaprograms/myeloid/", 
                                                      verbose = T, 
                                                      n_genes = 50)

# The myeloid_nmf_metaprograms object is a list containing the clusters (i.e. the NMF programs from which each MP was derived)
# and the MPs in tabular and list form

# saveRDS(myeloid_nmf_metaprograms, file= "/vast/palmer/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/myeloid/myeloid_nmf_metaprograms_out.RDS")
myeloid_nmf_metaprograms <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/myeloid/myeloid_nmf_metaprograms_out.RDS")

# Inspect the out and determine whether there is representation from the IDH-O tumors in some of the MPs
# meta from care to separate based on codel status.

# Visualize metaprograms:
library(tidyverse)
meta_care <- read.table("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/caremut_sample_identifier_linker.txt", sep = "\t", row.names = 1, header = TRUE)
mp_res_df <- data.frame(
  Cluster = rep(names(myeloid_nmf_metaprograms$clusters), sapply(myeloid_nmf_metaprograms$clusters, length)),
  Sample = unlist(myeloid_nmf_metaprograms$clusters)
)
mp_res_df$SampleID <- sapply(strsplit(mp_res_df$Sample, "_"), "[[", 1)
mp_res_df_annot <- mp_res_df %>% 
  inner_join(meta_care, "SampleID")

# What's the distribution of clusters across subtype and timepoint - 
table(mp_res_df_annot$Cluster, mp_res_df_annot$idh_codel_subtype)
table(mp_res_df_annot$Cluster, mp_res_df_annot$lab)
table(mp_res_df_annot$Cluster, mp_res_df_annot$timepoint)

mp_res_df_annot_distinct <- mp_res_df_annot %>% 
  dplyr::select(Cluster, care_id:timepoint, sample_barcode:idh_codel_subtype) %>% 
  distinct()
mp_by_subtype <- mp_res_df_annot_distinct %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "IDH-O",
                                    `IDHmut-noncodel` = "IDH-A")) %>% 
  group_by(idh_codel_subtype, Cluster) %>% 
  summarise(counts = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = idh_codel_subtype, values_from = counts)
mp_by_subtype[is.na(mp_by_subtype)] <- 0
mp_by_subtype$IDH_A_freq <- mp_by_subtype$`IDH-A` / 45
mp_by_subtype$IDH_O_freq <- mp_by_subtype$`IDH-O` / 30

mp_by_subtype_long <- mp_by_subtype %>% 
  dplyr::select(Cluster, IDH_A_freq, IDH_O_freq) %>% 
  mutate(IDH_O_freq = ifelse(is.na(IDH_O_freq), 0, IDH_O_freq)) %>% 
  pivot_longer(cols=c(IDH_A_freq, IDH_O_freq), names_to = "subtype", values_to = "freq") %>% 
  mutate(Cluster = gsub("Cluster_", "MP_", Cluster))

mp_by_subtype_long$Cluster <- factor(mp_by_subtype_long$Cluster, levels=c("MP_1", "MP_2", "MP_3", "MP_4", "MP_5", "MP_6",
                                                                          "MP_7", "MP_8", "MP_9", "MP_10", "MP_11", "MP_12",
                                                                          "MP_13", "MP_14", "MP_15"))
mp_by_subtype_long$subtype <- factor(mp_by_subtype_long$subtype, levels=c("IDH_O_freq", "IDH_A_freq"))

source("/vast/palmer/pi/verhaak/kcj28/care_mut/scripts/misc/plot_theme.R")
fig_dir <-"/vast/palmer/pi/verhaak/kcj28/care_mut/results/figures/classification/tme/"

pdf(paste0(fig_dir, "caremut_myeloid_metaprogram_sample_contribution_by_subtype.pdf"), width = 5, height = 4, useDingbats = FALSE)
ggplot(mp_by_subtype_long, aes(x=Cluster, y=freq*100, fill=subtype)) +
  geom_bar(position="dodge", stat="identity") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  geom_hline(yintercept = 25) +
  scale_fill_manual(values=c("IDH_A_freq" = "#67A9CF", 
                             "IDH_O_freq"="#EF8A62"),
                    labels=c("IDH_A_freq" = "IDH-A",
                             "IDH_O_freq" = "IDH-O")) +
  labs(x="Metaprogram", y="Percent contributing samples\n(proportional to subtype)")
dev.off()

mp_by_lab <- mp_res_df_annot_distinct %>% 
  filter(idh_codel_subtype=="IDHmut-noncodel") %>% 
  group_by(lab, Cluster) %>% 
  summarise(counts = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = lab, values_from = counts)
mp_by_lab$Iavarone_lab_freq <- mp_by_lab$`Iavarone lab` / 8
mp_by_lab$Suva_lab_freq <- mp_by_lab$`Suva lab` / 16
mp_by_lab$Verhaak_lab_freq <- mp_by_lab$`Verhaak lab` / 21

mp_by_lab_long <- mp_by_lab %>% 
  dplyr::select(Cluster, Iavarone_lab_freq:Verhaak_lab_freq) %>% 
  #mutate(IDH_O_freq = ifelse(is.na(IDH_O_freq), 0, IDH_O_freq)) %>% 
  pivot_longer(cols=c( Iavarone_lab_freq:Verhaak_lab_freq), names_to = "lab", values_to = "freq") %>% 
  mutate(Cluster = gsub("Cluster_", "MP_", Cluster))

mp_by_lab_long$Cluster <- factor(mp_by_lab_long$Cluster, levels=c("MP_1", "MP_2", "MP_3", "MP_4", "MP_5", "MP_6",
                                                                  "MP_7", "MP_8", "MP_9", "MP_10", "MP_11", "MP_12",
                                                                  "MP_13", "MP_14", "MP_15"))

pdf(paste0(fig_dir, "caremut_myeloid_metaprogram_noncodel_sample_contribution_by_lab.pdf"), width = 5, height = 4, useDingbats = FALSE)
ggplot(mp_by_lab_long, aes(x=Cluster, y=freq*100, fill=lab)) +
  geom_bar(position="dodge", stat="identity") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  scale_fill_manual(values =c( "Iavarone_lab_freq"= "#7570b3", "Suva_lab_freq"="#d95f02", "Verhaak_lab_freq"="#1b9e77"),
                    labels =c("Iavarone_lab_freq"="Iavarone(n=8)", "Verhaak_lab_freq"="Verhaak(n=21)", "Suva_lab_freq"="Suva(n=16)"),
                    name="Lab") +
  labs(x="Metaprogram", y="Percent contributing samples\n(proportional to lab)")
dev.off()


table(mp_res_df_annot_distinct$Cluster, mp_res_df_annot_distinct$lab)
table(mp_res_df_annot_distinct$Cluster, mp_res_df_annot_distinct$timepoint)


### END ###