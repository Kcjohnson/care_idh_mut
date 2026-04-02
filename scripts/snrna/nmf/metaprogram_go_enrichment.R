##################################
# Perform GO term enrichment for each MP (50 genes) for CAREmut malignant states
# Author: Kevin Johnson
# Date: 2026.04.01
##################################

library(tidyverse)
library(Seurat)
library(topGO)
library(org.Hs.eg.db)
library(dplyr)
library(pheatmap)
library(scales)
library(viridis)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
script_dir  <- file.path(proj_dir, "scripts")

setwd(proj_dir)

# Load CAREmut metaprogram information
malignant_metaprograms <- read.table(paste0(proj_dir, "/results/metaprograms/malignant_downsampled/meta_programs_generated_for_all_mutants_n74_min_groupsize5_2026-04-01.csv"), sep = ",", row.names = 1, header = T)
mp_list <- as.list(data.frame(malignant_metaprograms))

# The background gene set amongst which the malignant metaprograms genes would have been pulled from.
total_gene_set <- readRDS("data/snrna/unique_genes_across_malignant_nmf_input_matrices.RDS")

############################################
# For each MP run a TopGO analysis 
############################################
## Function for enrichment of significant NMF genes against covered background.
selFun = function(x) {
  ifelse(x==1, TRUE, FALSE)
}

goList <- list()
for(i in 1:length(mp_list)){
  cat("\r", i)
  mysig <- mp_list[[i]]
  mp_name <- names(mp_list)[i]
  all_genes <- total_gene_set
  
  gene_list_mp <- ifelse(all_genes%in%mysig, 1, 0)
  names(gene_list_mp) <- all_genes
  
  # Functional enrichment of metaprogram signature
  mpGOdata <- new("topGOdata",
                      ontology = "BP",
                      allGenes = gene_list_mp,
                      geneSel = selFun,
                      annot=annFUN.org,
                      mapping = 'org.Hs.eg.db', # The annotation package for the human genome.
                      ID = 'symbol', # We're using gene symbols.
                      nodeSize = 10)
  
  # Fishers test
  resultFisher <- runTest(mpGOdata, algorithm = "classic", statistic = "fisher")
  fishRes <- GenTable(mpGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
  fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")
  
  fishRes[,"MP"] <- mp_name
  goList[[i]] <- fishRes
}

goRes <- do.call(rbind, goList)
goRes$raw.p.value <- as.numeric(goRes$raw.p.value)

# Helpful to examine the larger number of gene sets. Some sets functionally repeat.
top_10_p_values <- goRes %>%
  group_by(MP) %>%
  top_n(10, -`raw.p.value`) %>% 
  filter(raw.p.value<0.05)

# For plotting purposes, it can be beneficial to simply plot the top hit either by p-value or q-value depending on strength of association.
top_1_q_values <- goRes %>%
  group_by(MP) %>%
  top_n(1, -q.value) %>% 
  filter(q.value<0.05)

top_5_p_values <- goRes %>%
  group_by(MP) %>%
  top_n(5, -`raw.p.value`) %>% 
  filter(raw.p.value<0.05)

terms_to_plot <- unique(top_10_p_values$Term)
goRes_plot <- goRes %>%
  filter(Term%in%terms_to_plot)

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

goRes_plot$MP <- factor(goRes_plot$MP, levels=mp_order)
goRes_plot$Term <- as.factor(goRes_plot$Term)
goRes_plot$Term <- factor(goRes_plot$Term, levels=terms_to_plot)

y_labels_to_display <- c("neuron recognition", "developmental process", 
                         "neurogenesis",
                         "protein metabolic process",
                         "cell cycle process", 
                         "response to wounding",
                         "neuron differentiation",
                         "chemical synaptic transmission",
                         "cell communication",
                         "response to stress")


pdf(paste0(fig_dir, "edf4c_malignant_mp_gene_ontology_enrichment.pdf"), width = 6, height = 4.5)
ggplot(goRes_plot, aes(x = MP, y = Term, fill = -log10(raw.p.value))) +
  geom_tile() +
  scale_fill_viridis(limits=c(0, 15), option = "D", oob=squish, name="Significance\n-log10(p-value)") +
  labs(y = "Enriched GO Term", x = "Malignant MP", fill="Significance\n-log10(p-value)") + 
  theme_bw() +
  theme(axis.ticks.x=element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1)) +
  scale_y_discrete(breaks = y_labels_to_display) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1), text=element_text(size=8))
dev.off()

### END ###