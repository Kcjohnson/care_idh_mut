##################################
# Test for enrichment of published gene sets among the 50 MP malignant genes
# Author: Kevin Johnson
# Date Updated: 2026.04.01
##################################

library(tidyverse)
library(openxlsx)
library(msigdbr)
library(scales)
library(viridis)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
script_dir  <- file.path(proj_dir, "scripts")

setwd(proj_dir)

source(paste0(script_dir, "/utils/metaprograms_enrichment.R"))

# Load in metaprograms and relevant gene sets from prior publications:

# Load final pan-IDH-mutant metaprogram gene set (snRNA)
care_mut_mps <- read.table(paste0(proj_dir, "/results/metaprograms/malignant_downsampled/meta_programs_generated_for_all_mutants_n74_min_groupsize5_2026-04-01.csv"), sep = ",", row.names = 1, header = T)
mut_mp_list <- lapply(names(care_mut_mps), function(col_name) care_mut_mps[[col_name]])
names(mut_mp_list) <-colnames(care_mut_mps)

# Load the CGGA dataset that contains the largest cohort of all IDH-mutant scRNA data
cgga_mut_mps <- read.table("data/misc/cgga_scrna_n46_malignant_metaprograms.csv", sep = ",", row.names = 1, header = T)
cgga_mut_mp_list <- lapply(names(cgga_mut_mps), function(col_name) cgga_mut_mps[[col_name]])
names(cgga_mut_mp_list) <-colnames(cgga_mut_mps)

# Load the CARE IDH-wildtype data for reference (Nomura et al Nature Genetics) processed by the same laboratories for snRNA.
care_wt_mps <- read.delim("data/misc/carewt_malignant_metaprograms.txt", sep = "\t", header = T)
colnames(care_wt_mps) <- gsub("MP_", "CARE_IDHwt_MP", colnames(care_wt_mps))

# Venteicher Astrocytoma
venteicher <- readWorkbook("data/misc/venteicher_table_s3.xlsx", startRow = 5, colNames = TRUE)
venteicher_signatures <- venteicher %>% 
  dplyr::select("Venteicher2017_OC"=`Oligo-program.(Fig..2C)`, "Venteicher2017_AC"=`Astro-program.(Fig..2C)`,
                "Venteicher2017_Stem"=`Stemness.program.(Fig..3C)`)
venteicher_signatures$Venteicher2017_OC <- trimws(venteicher_signatures$Venteicher2017_OC)
venteicher_signatures$Venteicher2017_AC <- trimws(venteicher_signatures$Venteicher2017_AC)
venteicher_signatures$Venteicher2017_Stem <- trimws(venteicher_signatures$Venteicher2017_Stem)
venteicher_sig_list <- as.list(venteicher_signatures[,1:3])
venteicher_sig_list <- lapply(venteicher_sig_list, function(x) x[!is.na(x)])

# Tirosh Oligodendroglioma
tirosh <- readWorkbook("data/misc/tirosh_nature_2016_supplementary_table_1.xlsx", startRow = 9, colNames = TRUE)
tirosh_signatures <- tirosh %>% 
  dplyr::select("Tirosh2016_OC"=`OC.(PCA-only)`,
                "Tirosh2016_AC"=`AC.(PCA-only)`,
                "Tirosh2016_Stem"=`stemness`)
tirosh_signatures$Tirosh2016_OC <- trimws(tirosh_signatures$Tirosh2016_OC)
tirosh_signatures$Tirosh2016_AC <- trimws(tirosh_signatures$Tirosh2016_AC)
tirosh_signatures$Tirosh2016_Stem <- trimws(tirosh_signatures$Tirosh2016_Stem)
tirosh_signatures_list <- as.list(tirosh_signatures[ ,1:3])
tirosh_signatures_list <- lapply(tirosh_signatures_list, function(x) x[!is.na(x)])

# Read in markers from Liu et al Cell 2023 - developing brain cells
liu_markers <- read.csv("data/misc/hNSPC_marker_genes.csv", header = T, stringsAsFactors = F) 

liu_markers <- apply(liu_markers, 2, function(x) {
  sapply(strsplit(x, split = "_"), function(y) y[[2]][1]) %>%
    head(100)
})
liu_markers_list <- setNames(lapply(1:ncol(liu_markers), function(i) liu_markers[, i]), colnames(liu_markers))
names(liu_markers_list) <- paste0("Liu2023_", names(liu_markers_list))

# Malignantn metaprograms from the pan-cancer analysis presented in Gavish et al Nature 2023
gavish_signatures <- readWorkbook("data/misc/gavish_nature_2023_signatures.xlsx", sheet = 1)
colnames(gavish_signatures) <- paste0("Gavish2023_", colnames(gavish_signatures))

# A set of previously curated brain signatures, including a reactive astrocyte signature.
load("data/misc/brain_signatures.RData")
names(brain_signatures)[c(1,4)]
astrocyte_list <- brain_signatures[c(1, 4)]
names(astrocyte_list) <- c("Sadick2022_protective_astro.", "Sadick2022_reactive_astro.")
astrocyte_list <- lapply(astrocyte_list, function(x) x[!is.na(x)])

# Load the MSigDB genesets that would be relevant to test for enrichment among
m_df_h <- msigdbr(species = "Homo sapiens", category = "H")
m_df_c5 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
m_df_c8 <- msigdbr(species = "Homo sapiens", category = "C8")

# Format the results
m_df_long <- m_df_c8 %>% 
  bind_rows(m_df_c5) %>% 
  bind_rows(m_df_h) %>% 
  dplyr::select(gs_name, human_gene_symbol) %>% 
  distinct() %>% 
  mutate(gs_name = paste0("MSigDB.", gs_name))
msigdb_genesets <- split(m_df_long$human_gene_symbol, m_df_long$gs_name)

# Combine all signatures to be tested.
all_signatures <- c(gavish_signatures, astrocyte_list, msigdb_genesets, liu_markers_list, tirosh_signatures_list, venteicher_sig_list, care_wt_mps)

# The background set from which genes could be selected for metaprogram enrichment.
# total_gene_set <- readRDS("data/misc/snrna_genes_in_seurat_obj.RDS")
total_gene_set <- readRDS("data/snrna/unique_genes_across_malignant_nmf_input_matrices.RDS")

res <- metaprograms_enrichment(mut_mp_list, pathways = all_signatures)

# Need to calculate a p-value to assign the level of enrichment in the overlap between the two gene sets
# Using a hypergeometric test to do this:
# Using format: phyper(Overlap-1, Group2, Total-Group2, Group1, lower.tail= FALSE)
res$hyper_value <- phyper(res$MP_Int-1, res$N_p, length(total_gene_set)-res$N_p, res$N_mp, lower.tail= FALSE)
res$p.adj <- p.adjust(res$hyper_value, method = "fdr")

# Select the top 5 five pathways to visualize in the heatmap
top_5_p_values <- res %>%
  group_by(MP) %>%
  top_n(5, -`p.adj`) %>% 
  filter(p.adj<0.05)

# Check a few interesting overlap sets
res$int_genes[res$MP=="MP_8"&res$Pathway=="CARE_IDHwt_MP8_GPC"]
res$int_genes[res$MP=="MP_5"&res$Pathway=="CARE_IDHwt_MP6_MES"]

terms_to_plot <- unique(top_5_p_values$Pathway)
res_plot <- res %>%
  filter(Pathway%in%terms_to_plot)

# Specifying the order in which the metaprograms should appear along the x-axis
mp_order <- c("MP_1",
              "MP_2",
              "MP_3",
              "MP_4",
              "MP_5",
              "MP_6",
              "MP_7",
              "MP_8",
              "MP_9",
              "MP_10",
              "MP_11",
              "MP_12")

res_plot$MP <- factor(res_plot$MP, levels=mp_order)
res_plot$Pathway <- as.factor(res_plot$Pathway)
res_plot$Pathway <- factor(res_plot$Pathway, levels=terms_to_plot)

ggplot(res_plot, aes(x = MP, y = Pathway, fill = -log10(hyper_value))) +
  geom_tile() +
  scale_fill_gradient2(limits=c(0,20), low ="white", high = "#756bb1", midpoint = 10, oob=squish, name="Significance\n-log10(adj. P)") +
  labs(y = "Enriched Pathway", x = "Malignant MP", fill="Significance\n-log10(p-value)") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1), text=element_text(size=8))

y_labels_to_display <- c("Liu2023_Pre.OPC", "Tirosh2016_OC", "CARE_IDHwt_MP2_OPC", 
                         "Liu2023_Astrocyte", "CARE_IDHwt_MP4_AC",
                         "CARE_IDHwt_MP3_CC",
                         "Gavish2023_MP16.MES.(glioma)", "Sadick2022_reactive_astro.", "CARE_IDHwt_MP6_MES",
                         "Tirosh2016_AC", "CARE_IDHwt_MP8_GPC","Gavish2023_MP25_Astrocytes", 
                         "CARE_IDHwt_MP10_Stress1")



pdf(paste0(fig_dir, "edf4d_malignant_hypergeometric_gene_overlap_enrichment.pdf"), width = 6, height = 4.5)
ggplot(res_plot, aes(x = MP, y = Pathway, fill = -log10(hyper_value))) +
  geom_tile() +
  scale_fill_viridis(limits=c(0,20), option = "D", oob=squish, name="Significance\n-log10(p-value)") +
  labs(y = "Enriched Genesets", x = "CARE IDH-mutant malignant MPs", fill="Significance\n-log10(p-value)") + 
  theme_bw() +
  theme(axis.ticks.x=element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1)) +
  scale_y_discrete(breaks = y_labels_to_display) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1), text=element_text(size=8))
dev.off()


### END ###