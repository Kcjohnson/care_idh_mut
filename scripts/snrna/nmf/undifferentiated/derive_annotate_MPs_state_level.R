##################################
# Derive and annotate MPs for Undifferentiated malignant cells in the CAREmut dataset
# Author: Kevin Johnson
# Date: 2024.04.15
##################################
proj_dir <- "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/"
fig_dir <- "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/figures/"
out_data_dir <- "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/processed_data/"
setwd(proj_dir)

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

# Load the NMF helper functions (implementation of the MP algorithm)
source("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/scripts/nmf/PvsR-NMF-caremut.R")

# Path to NMF results directory
nmf_res_dir <- "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/nmf_res/malignant_undifferentiated/"

# List files in the directory - there's 73/75 that had sufficient undifferentiated cells (>10) on which to run NMF based on input (SJ12-1 not included)
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
malignant_nmf_metaprograms <- derive_NMF_metaprograms(Genes_nmf_w_basis = Genes_nmf_w_basis, 
                                                      save_path = "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/undifferentiated/", 
                                                      population = "undifferentiated",
                                                      verbose = T, 
                                                      n_genes = 50)

# Save output as RDS file
saveRDS(malignant_nmf_metaprograms, file= "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/undifferentiated/undifferentiated_nmf_metaprograms_out.RDS")
# malignant_nmf_metaprograms <- readRDS("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/undifferentiated/undifferentiated_nmf_metaprograms_out.RDS")


# Visualize metaprograms:
library(tidyverse)
meta_care <- read.table("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/data/metadata/caremut_sample_identifier_linker.txt", sep = "\t", row.names = 1, header = TRUE)
mp_res_df <- data.frame(
  Cluster = rep(names(malignant_nmf_metaprograms$clusters), sapply(malignant_nmf_metaprograms$clusters, length)),
  Sample = unlist(malignant_nmf_metaprograms$clusters)
)
mp_res_df$SampleID <- sapply(strsplit(mp_res_df$Sample, "_"), "[[", 1)

mp_res_df_annot <- mp_res_df %>% 
  inner_join(meta_care, "SampleID")

# What's the distribution of clusters across subtype and timepoint 
table(mp_res_df_annot$Cluster, mp_res_df_annot$idh_codel_subtype)  # Cluster 4 (NPC-like) picked up mostly for IDH-A
table(mp_res_df_annot$Cluster, mp_res_df_annot$lab)
table(mp_res_df_annot$Cluster, mp_res_df_annot$timepoint) # Cluster 4 (NPC-like) not dominated by T2 

mp_res_df_annot_distinct <- mp_res_df_annot %>% 
  dplyr::select(Cluster, care_id:timepoint, sample_barcode:idh_codel_subtype) %>% 
  distinct()
sample_subtype_uniq <- mp_res_df_annot %>% 
  dplyr::select(sample_barcode:idh_codel_subtype) %>% 
  distinct()
mp_by_subtype <- mp_res_df_annot_distinct %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "IDH-O",
                                    `IDHmut-noncodel` = "IDH-A")) %>% 
  group_by(idh_codel_subtype, Cluster) %>% 
  summarise(counts = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = idh_codel_subtype, values_from = counts)
mp_by_subtype$IDH_A_freq <- mp_by_subtype$`IDH-A` / 45
mp_by_subtype$IDH_O_freq <- mp_by_subtype$`IDH-O` / 28

mp_by_subtype_long <- mp_by_subtype %>% 
  dplyr::select(Cluster, IDH_A_freq, IDH_O_freq) %>% 
  mutate(IDH_O_freq = ifelse(is.na(IDH_O_freq), 0, IDH_O_freq)) %>% 
  pivot_longer(cols=c(IDH_A_freq, IDH_O_freq), names_to = "subtype", values_to = "freq") %>% 
  mutate(Cluster = gsub("Cluster_", "MP_", Cluster))

mp_by_subtype_long$Cluster <- factor(mp_by_subtype_long$Cluster, levels=c("MP_1", "MP_2", "MP_3", "MP_4", "MP_5", "MP_6",
                                                                          "MP_7", "MP_8", "MP_9", "MP_10"))
mp_by_subtype_long$subtype <- factor(mp_by_subtype_long$subtype, levels=c("IDH_O_freq", "IDH_A_freq"))

source("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/scripts/misc/plot_theme.R")
fig_dir <-"/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/figures/classification/"

pdf(paste0(fig_dir, "undifferentiated_metaprogram_sample_contribution_by_subtype.pdf"), width = 5, height = 4, useDingbats = FALSE)
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
  labs(x="Undifferentiated malignant metaprogram", y="Percent contributing samples\n(proportional to subtype)", fill="Subtype")
dev.off()



### ### ### ### ### ### ### ### ### ### ### ### ###
## Repeat for OPC-like malignant cells
### ### ### ### ### ### ### ### ### ### ### ### ### 
nmf_res_dir <- "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/nmf_res/malignant_opc/"

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
opc_nmf_metaprograms <- derive_NMF_metaprograms(Genes_nmf_w_basis = Genes_nmf_w_basis, 
                                                      save_path = "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/opc/", 
                                                      population = "opc",
                                                      verbose = T, 
                                                      n_genes = 50)

# Save output as RDS file
saveRDS(opc_nmf_metaprograms, file= "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/opc/opc_nmf_metaprograms_out.RDS")
#opc_nmf_metaprograms <- readRDS("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/opc/opc_nmf_metaprograms_out.RDS")

### END ###