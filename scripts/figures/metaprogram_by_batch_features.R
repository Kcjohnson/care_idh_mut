##################################
# Visualize the metaprograms across different batch features (tumor type and lab)
# Author: Kevin Johnson
# Date: 2026.04.01
##################################

# Libraries, paths, and helper functions
library(tidyverse)


proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures")

source(file.path(proj_dir, "scripts/utils", "plot_theme.R"))

# Load in the the derive_NMF_metaprograms() function.
# malignant_nmf_metaprograms <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/malignant_nrun10/malignant_nmf_metaprograms_out.RDS")
# The malignant_nmf_metaprograms object is a list containing the clusters (i.e. the NMF programs from which each MP was derived)
# and the MPs in tabular and in list form

# Use the simplified linker
meta_care <- read.table("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/caremut_sample_identifier_linker.txt", sep = "\t", row.names = 1, header = TRUE)
mp_res_df <- data.frame(
  Cluster = rep(names(malignant_nmf_metaprograms$clusters), sapply(malignant_nmf_metaprograms$clusters, length)),
  Sample = unlist(malignant_nmf_metaprograms$clusters)
)
mp_res_df$SampleID <- sapply(strsplit(mp_res_df$Sample, "_"), "[[", 1)
mp_res_df_annot <- mp_res_df %>% 
  inner_join(meta_care, "SampleID")

# Examine the distribution of clusters across subtype and time point.
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
mp_by_subtype$IDH_A_freq <- mp_by_subtype$`IDH-A` / 44
mp_by_subtype$IDH_O_freq <- mp_by_subtype$`IDH-O` / 29

mp_by_subtype_long <- mp_by_subtype %>% 
  dplyr::select(Cluster, IDH_A_freq, IDH_O_freq) %>% 
  mutate(IDH_O_freq = ifelse(is.na(IDH_O_freq), 0, IDH_O_freq)) %>% 
  pivot_longer(cols=c(IDH_A_freq, IDH_O_freq), names_to = "subtype", values_to = "freq") %>% 
  mutate(Cluster = gsub("Cluster_", "MP_", Cluster))

mp_by_subtype_long$Cluster <- factor(mp_by_subtype_long$Cluster, levels=c("MP_1", "MP_2", "MP_3", "MP_4", "MP_5", "MP_6",
                                                                          "MP_7", "MP_8", "MP_9", "MP_10", "MP_11", "MP_12"))
mp_by_subtype_long$subtype <- factor(mp_by_subtype_long$subtype, levels=c("IDH_O_freq", "IDH_A_freq"))



pdf(paste0(fig_dir, "caremut_metaprogram_sample_contribution_by_tumor_type.pdf"), width = 5, height = 4, useDingbats = FALSE)
ggplot(mp_by_subtype_long, aes(x=Cluster, y=freq*100, fill=subtype)) +
  geom_bar(position="dodge", stat="identity") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  scale_fill_manual(values=c("IDH_A_freq" = "#800074", 
                             "IDH_O_freq"="#298C8C"),
                    labels=c("IDH_A_freq" = "Astro.",
                             "IDH_O_freq" = "Oligo.")) +
  labs(x="Metaprogram", y="Percent contributing samples\n(proportional to subtype)", fill="Subtype")
dev.off()

# Restrict to Astrocytoma samples.
mp_by_lab <- mp_res_df_annot_distinct %>% 
  filter(idh_codel_subtype=="IDHmut-noncodel") %>% 
  group_by(lab, Cluster) %>% 
  summarise(counts = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = lab, values_from = counts)
mp_by_lab$Iavarone_lab_freq <- mp_by_lab$`Iavarone lab` / 8
mp_by_lab$Suva_lab_freq <- mp_by_lab$`Suva lab` / 16
mp_by_lab$Verhaak_lab_freq <- mp_by_lab$`Verhaak lab` / 20

mp_by_lab_long <- mp_by_lab %>% 
  dplyr::select(Cluster, Iavarone_lab_freq:Verhaak_lab_freq) %>% 
  #mutate(IDH_O_freq = ifelse(is.na(IDH_O_freq), 0, IDH_O_freq)) %>% 
  pivot_longer(cols=c( Iavarone_lab_freq:Verhaak_lab_freq), names_to = "lab", values_to = "freq") %>% 
  mutate(Cluster = gsub("Cluster_", "MP_", Cluster))

mp_by_lab_long$Cluster <- factor(mp_by_lab_long$Cluster, levels=c("MP_1", "MP_2", "MP_3", "MP_4", "MP_5", "MP_6",
                                                                  "MP_7", "MP_8", "MP_9", "MP_10", "MP_11", "MP_12"))

pdf(paste0(fig_dir, "caremut_metaprogram_noncodel_sample_contribution_by_lab.pdf"), width = 5, height = 4, useDingbats = FALSE)
ggplot(mp_by_lab_long, aes(x=Cluster, y=freq*100, fill=lab)) +
  geom_bar(position="dodge", stat="identity") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  scale_fill_manual(values =c( "Iavarone_lab_freq"= "#7570b3", "Suva_lab_freq"="#d95f02", "Verhaak_lab_freq"="#1b9e77"),
                    labels =c("Iavarone_lab_freq"="Iavarone(n=8)", "Verhaak_lab_freq"="Verhaak(n=21)", "Suva_lab_freq"="Suva(n=16)"),
                    name="IDH-A (lab)") +
  labs(x="Metaprogram", y="Percent contributing samples\n(proportional to lab)") 
dev.off()


table(mp_res_df_annot_distinct$Cluster, mp_res_df_annot_distinct$lab)
table(mp_res_df_annot_distinct$Cluster, mp_res_df_annot_distinct$timepoint)

mp_res_df_annot_filt <- mp_res_df_annot %>% 
  filter(Cluster%in%c("Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4", "Cluster_5", "Cluster_6", "Cluster_7", "Cluster_8", "Cluster_9", "Cluster_10", "Cluster_11", "Cluster_12")) %>% 
  dplyr::select(Cluster:Sample, SampleID)

nmf_intersect_KEEP <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/malignant_nrun10/nmf_intersect_KEEP.RDS")
inds_new <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/malignant_nrun10/inds_new.RDS")

tmp <- nmf_intersect_KEEP[inds_new,inds_new]
tmp_filtered <- tmp[rownames(tmp)%in%mp_res_df_annot_filt$Sample, colnames(tmp)%in%mp_res_df_annot_filt$Sample]
nmf_intersect_meltI_NEW_compare <- reshape2::melt(tmp_filtered)

nmf_intersect_meltI_NEW <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/malignant_nrun10/nmf_intersect_meltI_NEW.RDS")

#mp_res_df_annot_filt$Sample <- factor(mp_res_df_annot_filt$Sample, levels=levels(nmf_intersect_meltI_NEW$Var1))

nmf_intersect_meltI_NEW$present_in_cluster <- ifelse(nmf_intersect_meltI_NEW$Var1%in%mp_res_df_annot_filt$Sample, "yes", "missing")
all(nmf_intersect_meltI_NEW$present_in_cluster=="yes")
missing_indices <- which(nmf_intersect_meltI_NEW$present_in_cluster!="yes")

nmf_intersect_meltI_NEW[367:377, ]

nmf_intersect_meltI_NEW_annot <- nmf_intersect_meltI_NEW %>%
  filter(Var1%in%mp_res_df_annot_filt$Sample, Var2%in%mp_res_df_annot_filt$Sample)

all(nmf_intersect_meltI_NEW_annot$Var1==nmf_intersect_meltI_NEW$Var1[1:dim(nmf_intersect_meltI_NEW_annot)[1]])
all(nmf_intersect_meltI_NEW_annot$Var2==nmf_intersect_meltI_NEW$Var2[1:dim(nmf_intersect_meltI_NEW_annot)[1]])


p <- ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  #scale_color_gradient2(limits=c(1,50), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 25, oob=squish, name="Similarity\n(Jaccard index)") +
  #scale_fill_gradient2(limits=c(1,50), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 25, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), 
         axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

p2 <- ggplot(data = nmf_intersect_meltI_NEW_annot, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  #scale_color_gradient2(limits=c(1,50), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 25, oob=squish, name="Similarity\n(Jaccard index)") +
  #scale_fill_gradient2(limits=c(1,50), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 25, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), 
         axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

p3 <- ggplot(data = nmf_intersect_meltI_NEW_compare, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  #scale_color_gradient2(limits=c(1,50), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 25, oob=squish, name="Similarity\n(Jaccard index)") +
  #scale_fill_gradient2(limits=c(1,50), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 25, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), 
         axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

meta_care <- read.table("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/caremut_sample_identifier_linker.txt", sep = "\t", row.names = 1, header = TRUE)
rownames(meta_care) <- meta_care$SampleID

nmf_annot_type_tumor  <- data.frame("sample"=sub("\\_.*", "", colnames(tmp_filtered)))
nmf_annot_type_tumor$idh_codel_subtype <- meta_care$idh_codel_subtype[match(nmf_annot_type_tumor$sample, rownames(meta_care))]
nmf_annot_type_tumor$idh_codel_subtype <- factor(nmf_annot_type_tumor$idh_codel_subtype, levels=sort(unique(nmf_annot_type_tumor$idh_codel_subtype)))
nmf_annot_type_tumor$timepoint <- meta_care$timepoint[match(nmf_annot_type_tumor$sample, rownames(meta_care))]
nmf_annot_type_tumor$timepoint <- factor(nmf_annot_type_tumor$timepoint, levels=sort(unique(nmf_annot_type_tumor$timepoint)))
nmf_annot_type_tumor$lab <- meta_care$lab[match(nmf_annot_type_tumor$sample, rownames(meta_care))]
nmf_annot_type_tumor$lab <- factor(nmf_annot_type_tumor$lab, levels=sort(unique(nmf_annot_type_tumor$lab)))


#nmf_annot_type_tumor_old <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/malignant_nrun10/nmf_annot_type_tumor.RDS")


p_annotation_bottom <- ggplot(nmf_annot_type_tumor, aes(y="", x=1:nrow(nmf_annot_type_tumor), fill=idh_codel_subtype, color=idh_codel_subtype)) +
  geom_tile() +
  scale_color_manual(values =c( "IDHmut-codel"= "#EF8A62", "IDHmut-noncodel"="#67A9CF"), name="") +
  scale_fill_manual(values =c( "IDHmut-codel"= "#EF8A62", "IDHmut-noncodel"="#67A9CF"), name="") +                                                    
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10),  axis.title = element_text(size=8),  legend.text.align = 0, legend.key.size = unit(0.4, "cm"), legend.key=element_blank(), legend.position=c(1.3,-4), plot.margin = unit(c(0.5,3,-0.6,0.5), "cm")) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  guides(fill = "none", color="none")

p_annotation_middle <- ggplot(nmf_annot_type_tumor, aes(y="", x=1:nrow(nmf_annot_type_tumor), fill=lab, color=lab)) +
  geom_tile() +
  scale_color_manual(values =c( "Iavarone lab"= "#7570b3", "Suva lab"="#d95f02", "Verhaak lab"="#1b9e77"), name="") +
  scale_fill_manual(values =c( "Iavarone lab"= "#7570b3", "Suva lab"="#d95f02", "Verhaak lab"="#1b9e77"), name="") +                                                    
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10),  axis.title = element_text(size=8),  legend.text.align = 0, legend.key.size = unit(0.4, "cm"), legend.key=element_blank(), legend.position=c(1.3,-4), plot.margin = unit(c(0.5,3,-0.6,0.5), "cm")) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  guides(fill = "none", color="none")

p_annotation_top <- ggplot(nmf_annot_type_tumor, aes(y="", x=1:nrow(nmf_annot_type_tumor), fill=timepoint, color=timepoint)) +
  geom_tile() +
  scale_color_manual(values =c( "T1"= "#f1eef6", "T2"="#2b8cbe", "T3"="#045a8d"), name="") +
  scale_fill_manual(values =c( "T1"= "#f1eef6", "T2"="#2b8cbe", "T3"="#045a8d"), name="") +                                                    
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10),  axis.title = element_text(size=8),  legend.text.align = 0, legend.key.size = unit(0.4, "cm"), legend.key=element_blank(), legend.position=c(1.3,-4), plot.margin = unit(c(0.5,3,-0.6,0.5), "cm")) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  guides(fill = "none", color="none")

# Combine plots
pdf(paste0(fig_dir, "similarity_heatmap_mp_care_mut_annot_all_malignant_11clusters.pdf"), width=6, height=4, onefile=FALSE)
egg::ggarrange(p_annotation_middle, p_annotation_bottom, p3, nrow = 3, heights = c(1, 1, 8))
dev.off()

# Combine plots
pdf(paste0(fig_dir, "similarity_heatmap_mp_care_mut_annot_all_malignant_11clusters_subtype.pdf"), width=6, height=4, onefile=FALSE)
egg::ggarrange(p_annotation_bottom, p3, nrow = 2, heights = c(1, 8))
dev.off()

### ### ### ### ###