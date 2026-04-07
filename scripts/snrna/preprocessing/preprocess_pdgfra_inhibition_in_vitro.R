##################################
# 10x scRNA processing of NORLUX patient-derived neurospheres treated with Dasatinib, CP-673451, and DMSO
# Author: Kevin Johnson
# Date Updated: 2026.04.06
##################################

library(tidyverse) # 1.3.1 
library(Seurat) # 4.3.0
library(Matrix) # 1.6-5
library(ggpubr) # 0.4.0 
library(DoubletFinder) # 2.0.3
library(openxlsx) # 4.2.5.2 

# Notes: These results were based on cellranger 9.0.1 and used the 10x Genomics 2020 reference transcriptome.

# Define where figures and results will be stored:
fig_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/figures/perturbation/"
out_data_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/perturbation/pdgfrai/"
setwd("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/")
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/caremut_utils.R")

# Specify directories that have the count matrices. 
in_data_dir <- "/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/cellranger/pdgfrai/9.0.1/"

# Metadata that maps individual file names (e.g., A7) to experimental conditions (T394NS DMSO).
norlux_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/perturbation/20250318_scrnaseq_10xgenomics_metadata_Yale.txt", sep="\t", header = TRUE)
norlux_md_trim <- norlux_md %>% 
  mutate(SampleID = paste0("A", Sample.number),
         cell_line = sapply(strsplit(Sample.name.code, "_"), "[[", 1),
         treatment = sapply(strsplit(Sample.name.code, "_"), "[[", 2),
         dose = sapply(strsplit(Sample.name.code, "_"), "[[", 3)) %>% 
  dplyr::select(SampleID, exp_group = Sample.name.code, cell_line:dose)

## List all samples in each directory.
SAMPLES <- dir(in_data_dir, full.names = T)

# Create function to grab the last element of each vector (varying length vectors) from file path.
last <- function(x) { return( x[length(x)] ) }
SAMPLES_NAMES <- sapply(strsplit(SAMPLES, "/"), last)

seurat_list <- lapply(SAMPLES_NAMES, function(sample_name){
  cur_data <- Read10X(paste0(in_data_dir, "/", sample_name,'/outs/filtered_feature_bc_matrix/'))
  
  ## Adding the sample identifier to the cell barcode.
  colnames(cur_data) <- paste0(sample_name, "-", colnames(cur_data))
  
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3,
    min.features=200,
    project='pdgfrai'
  )
  cur_seurat$SampleID <- sample_name
  return(cur_seurat)
})

# Merge seurat objects. 
seurat_obj <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])

# Clean up large objects
rm(seurat_list)
gc()

# Make consistent with naming scheme from other runs:
seurat_obj@meta.data$CellID <- rownames(seurat_obj@meta.data)
seurat_obj@meta.data$SampleID <- sapply(strsplit(seurat_obj@meta.data$CellID, "-"), "[[", 1)

#########################################################################################
##### QC assessment - pre-filtering out low quality cells
#########################################################################################
# Perform quality control assessment. Assign a percent mitochondrial RNA to the metadata.
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")

pdf(paste0(fig_dir, "qc_mito_violin_pdgfrai.pdf"), width=10, height=6)
VlnPlot(seurat_obj, group.by="SampleID", features = "percent.mt", ncol = 1, pt.size=0) + NoLegend() + ylim(0, 10)
dev.off()

pdf(paste0(fig_dir, "qc_nfeature_rna_violin_pdgfrai.pdf"), width=10, height=6)
VlnPlot(seurat_obj, group.by="SampleID", features = "nFeature_RNA", ncol = 1, pt.size=0) + NoLegend()
dev.off()

# In the past, we set a relatively liberal criteria in tumor snRNAseq data of 500 genes detected so as not to miss lymphocytes, which do not express many genes.
# These data were deeply sequenced and come from single cell cell line models. 
# Based on the distribution above, we can set a higher threshold to select the best cells.
sum(seurat_obj$nFeature_RNA>2000) 
sum(seurat_obj$nFeature_RNA<10000) 
sum(seurat_obj$nFeature_RNA>10000) 
hist(seurat_obj$percent.mt)
abline(v = 7.5, col = "red", lwd = 2)
sum(seurat_obj$percent.mt>7.5) # Inspect how many cells would we lose depending on percent.mt threshold. 5-20% may all be valid.

# Assess the changes in cell number following QC. Filter cells based on general quality assessment in scRNAseq data
pre_qc_cell_num <- dim(seurat_obj)[2]

# Limiting to fewer than 10K genes to remove potential doublets plus setting a threshold for mitochondrial genes and a minimum number of genes.
# Note that these values differ from other experiments due to high data quality.
seurat_obj <- subset(seurat_obj, nFeature_RNA > 2000 & nFeature_RNA < 10000 & percent.mt < 7.5)

# 40,653 cells were removed. 86,396 cells leftover.
dropped_cells <- pre_qc_cell_num-dim(seurat_obj)[2]
sprintf("%s cells were removed due to low quality", dropped_cells)


#########################################################################################
##### Doublet Finder
#########################################################################################

# I followed the following example: https://rpubs.com/kenneditodd/doublet_finder_example

# Split aggregated data by sample
pdgfrai_split <- SplitObject(seurat_obj, split.by = "SampleID") 

# loop through samples to find doublets
for (i in 1:length(pdgfrai_split)) {
  # print the sample we are on
  print(paste0("Sample ",i))
  
  # Pre-process seurat object with standard seurat workflow
  pdgfrai_sample <- NormalizeData(pdgfrai_split[[i]])
  pdgfrai_sample <- FindVariableFeatures(pdgfrai_sample)
  pdgfrai_sample <- ScaleData(pdgfrai_sample)
  pdgfrai_sample <- RunPCA(pdgfrai_sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- pdgfrai_sample[["pca"]]@stdev
  sum.stdv <- sum(pdgfrai_sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # finish pre-processing
  pdgfrai_sample <- RunUMAP(pdgfrai_sample, dims = 1:min.pc)
  pdgfrai_sample <- FindNeighbors(object = pdgfrai_sample, dims = 1:min.pc)              
  pdgfrai_sample <- FindClusters(object = pdgfrai_sample, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep_v3(pdgfrai_sample)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bimodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- pdgfrai_sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(pdgfrai_sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  pdgfrai_sample <- doubletFinder_v3(seu = pdgfrai_sample, 
                                     PCs = 1:min.pc, 
                                     pK = optimal.pk,
                                     nExp = nExp.poi.adj)
  metadata <- pdgfrai_sample@meta.data
  # Get the number of columns in the data frame
  num_columns <- ncol(metadata)
  
  # Change the name of the last column
  colnames(metadata)[num_columns] <- "doublet_finder"
  
  pdgfrai_sample@meta.data <- metadata 
  
  # subset and save
  pdgfrai_singlets <- subset(pdgfrai_sample, doublet_finder == "Singlet")
  pdgfrai_split[[i]] <- pdgfrai_singlets
  remove(pdgfrai_singlets)
}

print("Finished detecting singlets")


# Converge the co-culture splits after identifying the singlets.
seurat_obj_singlets <- merge(x=pdgfrai_split[[1]], 
                             y=pdgfrai_split[2:length(pdgfrai_split)],
                             project = "pdgfrai")

columns_to_remove <- grep("^pANN", colnames(seurat_obj_singlets@meta.data), value = TRUE)
seurat_obj_singlets@meta.data <- seurat_obj_singlets@meta.data[, !colnames(seurat_obj_singlets@meta.data) %in% columns_to_remove]

## Save the files
saveRDS(seurat_obj_singlets, paste0(out_data_dir, "pdgfrai_seurat_obj_singlets_20260226.RDS"))

# seurat_obj_singlets <- readRDS(paste0(out_data_dir, "pdgfrai_seurat_obj_singlets_20260226.RDS"))


#########################################################################################
##### Preprocess with Standard processing  
#########################################################################################

# Add the extra metadata to the Seurat object about the experimental conditions.
seurat_obj_singlets@meta.data$CellID <- rownames(seurat_obj_singlets@meta.data)
seurat_md <- seurat_obj_singlets@meta.data
add_data <- seurat_md %>% 
  left_join(norlux_md_trim, by="SampleID")
# Add CellID back
rownames(add_data) <-  add_data$CellID

# sanity check
all(row.names(seurat_obj_singlets[[]])==add_data$CellID)
row.names(seurat_obj_singlets[[]])==add_data$CellID
head(rownames(seurat_obj_singlets[[]]))
head(add_data$CellID)

row.names(add_data) <- row.names(seurat_obj_singlets[[]])
seurat_obj_singlets <- AddMetaData(seurat_obj_singlets, metadata = add_data)

pdf(paste0(fig_dir, "pass_qc_mito_rna_violin_pdgfrai.pdf"), width=3.5, height=3, useDingbats = FALSE)
VlnPlot(seurat_obj_singlets, group.by="exp_group", features = "percent.mt", ncol = 1, pt.size=0) + 
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    color = "black",
    fill = "white"
  )  +
  geom_hline(yintercept = 7.5, linetype = "dashed", color="red") +
  NoLegend() + ylim(0, 15) + xlab("Sample ID") + ylab("Percent mitochrondrial genes (%)") + theme(plot.title = element_blank()) +
  plot_theme
dev.off()

pdf(paste0(fig_dir, "pass_qc_nfeature_rna_violin_pdgfrai.pdf"), width=3.5, height=3, useDingbats = FALSE)
VlnPlot(seurat_obj_singlets, group.by="exp_group", features = "nFeature_RNA", ncol = 1, pt.size=0) + 
  geom_hline(yintercept = c(2000, 10000), linetype = "dashed", color="red") +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    color = "black",
    fill = "white"
  )  +
  NoLegend() + ylim(0, 10000) + xlab("Sample ID") + ylab("nFeature RNA") + theme(plot.title = element_blank()) +
  plot_theme
dev.off()


# Run through a standard Seurat workflow
pdgfrai_sobj_all <- NormalizeData(seurat_obj_singlets, normalization.method = "LogNormalize", scale.factor = 10000)
pdgfrai_sobj_all <- FindVariableFeatures(pdgfrai_sobj_all, selection.method = "vst", nfeatures = 5000)
# By default only variable features are scaled
all.genes <- rownames(pdgfrai_sobj_all)
pdgfrai_sobj_all <- ScaleData(pdgfrai_sobj_all, features = all.genes)
pdgfrai_sobj_all <- RunPCA(pdgfrai_sobj_all, features = VariableFeatures(object = pdgfrai_sobj_all), npcs=50)
pdgfrai_sobj_all <- FindNeighbors(pdgfrai_sobj_all, dims = 1:20) 
pdgfrai_sobj_all <- FindClusters(pdgfrai_sobj_all, resolution = 0.6) 
pdgfrai_sobj_all <- RunUMAP(pdgfrai_sobj_all, dims = 1:20, min.dist = 0.2, n.neighbors = 40)

options(ggrepel.max.overlaps = Inf)
umap_plot_std <- DimPlot(pdgfrai_sobj_all, label = TRUE, group.by = "exp_group", repel = TRUE, label.size = 2.25) + 
  xlab("UMAP 1") + ylab("UMAP 2") + plot_theme + theme(plot.title = element_blank()) 

pdf(paste0(fig_dir, "pass_qc_standard_preprocesing_umap_pdgfrai.pdf"), width=4, height=4, useDingbats = FALSE)
umap_plot_std
dev.off()

ggsave(paste0(fig_dir, "pass_qc_standard_preprocesing_umap_pdgfrai.png"), umap_plot_std, width = 4.25, height = 3, dpi = 300)

# Extract the gene names to see how they are structured.
gene_names <- rownames(seurat_obj_singlets[["RNA"]])

# Assign cell cycle phase status based on Seurat's scoring approach for the Neftel et al.
neftel_sig <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/misc/neftel_metamodule_genelists.csv", sep=",", header = TRUE)
s.genes <- neftel_sig$G1.S[!is.na(neftel_sig$G1.S)]  
g2m.genes <- neftel_sig$G2.M[!is.na(neftel_sig$G2.M)]  
pdgfrai_sobj_all <- CellCycleScoring(pdgfrai_sobj_all, s.features = s.genes, g2m.features = g2m.genes)

# View cell cycle scores and phase assignments
seurat_md_cc <- pdgfrai_sobj_all@meta.data

cell_cycle_phase <- seurat_md_cc %>% 
  group_by(exp_group, cell_line, Phase) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup()

cell_cycle_phase$cell_line <- sapply(strsplit(cell_cycle_phase$exp_group, "_"), "[[", 1)
cell_cycle_phase$drug <- sapply(strsplit(cell_cycle_phase$exp_group, "_"), "[[", 2)
cell_cycle_phase$dose <- sapply(strsplit(cell_cycle_phase$exp_group, "_"), "[[", 3)
cell_cycle_phase$Phase <- factor(cell_cycle_phase$Phase, levels=c("G1", "S", "G2M"))
cell_cycle_phase$drug_dose <- paste0(cell_cycle_phase$drug, "_",cell_cycle_phase$dose)

pdf(paste0(fig_dir, "pdgfrai_malignant_cellcycle_assignment_stacked.pdf"), width = 8, height = 6, bg = "transparent")
cell_cycle_phase %>% 
  ggplot(aes(x=drug_dose, fill = factor(Phase), y=freq*100)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "Relative cell abundance (%)", fill="Cell Cycle\nPhase", x = "") +
  scale_fill_manual(values = c("G1"= "#67a9cf",
                               "S" = "gray80",
                               "G2M" = "#ef8a62")) +
  facet_grid(.~cell_line, scales = "free_x", space = "free_x") +
  plot_theme +
  theme(panel.grid.major = element_line(),
        axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

pdf(paste0(fig_dir, "pdgfrai_malignant_cellcycle_assignment_scores_downsampled.pdf"), width = 8, height = 6, bg = "transparent")
seurat_md_cc %>%
  group_by(cell_line) %>%
  slice_sample(n = 500) %>%
  ungroup() %>%
  ggplot(aes(x = S.Score, y = G2M.Score)) +
  geom_point(aes(color = Phase)) +
  labs(x="G1/S score", y = "G2/M score", fill = "Cell Cycle\nPhase") +
  scale_fill_manual(values = c("G1" = "#67a9cf",
                               "S" = "gray80",
                               "G2M" = "#ef8a62")) +
  facet_grid(. ~ cell_line, scales = "free_x", space = "free_x") +
  plot_theme
dev.off()

# Data exploration for key features.
options(ggrepel.max.overlaps = Inf)
DimPlot(pdgfrai_sobj_all, label = TRUE, group.by = "RNA_snn_res.0.6", repel = TRUE) + NoLegend()

# Some separation by cell line
DimPlot(pdgfrai_sobj_all, label = TRUE, group.by = "cell_line", repel = TRUE)
# Only a minor difference in drugs and not anything clear from a dose perspective.
DimPlot(pdgfrai_sobj_all, label = TRUE, group.by = "treatment", repel = TRUE)
DimPlot(pdgfrai_sobj_all, label = TRUE, group.by = "dose", repel = TRUE)

# Inspect quality control metrics across the different experimental conditions.
VlnPlot(pdgfrai_sobj_all, features = "percent.mt", group.by = "treatment")
VlnPlot(pdgfrai_sobj_all, features = "percent.mt", group.by = "exp_group")

pdgfrai_sobj_all[["percent.ribo"]] <- PercentageFeatureSet(object = pdgfrai_sobj_all, pattern = "^RPS|^RPL")
VlnPlot(pdgfrai_sobj_all, features = "percent.ribo", group.by = "exp_group")
# A select group of markers: SOX2 (general stemness), CD44 (MES-like), TOP2A (cycling)
FeaturePlot(pdgfrai_sobj_all, features = c("SOX2", "CD44", "TOP2A"))

# PDGFRA gene expression level does not seem to be decreased by the various treatments.
VlnPlot(pdgfrai_sobj_all, features = "PDGFRA", group.by = "exp_group")

# Perform differential expression by the different experimental groups. 
Idents(pdgfrai_sobj_all) <- "exp_group"

T394NS_cp_5uM <- FindMarkers(
  object = pdgfrai_sobj_all,
  ident.1 = "T394NS_CP-673451_5uM",
  ident.2 = "T394NS_DMSO_CTR")

T407NS_cp_5uM <- FindMarkers(
  object = pdgfrai_sobj_all,
  ident.1 = "T407NS_CP-673451_5uM",
  ident.2 = "T407NS_DMSO_CTR")

# Compare the top upregulated hits across cell lines for the CP-673451 5uM.
T394NS_cp_5uM_hits <- T394NS_cp_5uM %>% 
  filter(avg_log2FC > 0)
colnames(T394NS_cp_5uM_hits) <- paste0("T394", colnames(T394NS_cp_5uM_hits))
T394NS_cp_5uM_hits$gene_name <- rownames(T394NS_cp_5uM_hits)

T407NS_cp_5uM_hits <- T407NS_cp_5uM %>% 
  filter(avg_log2FC > 0)
colnames(T407NS_cp_5uM_hits) <- paste0("T407", colnames(T407NS_cp_5uM_hits))
T407NS_cp_5uM_hits$gene_name <- rownames(T407NS_cp_5uM_hits)

# Are there common genes? JUN/FOS/PDEs lots of indicators of stress response and potential resistance
cp_hits <- T394NS_cp_5uM_hits %>% 
  inner_join(T407NS_cp_5uM_hits, by="gene_name")

# Repeat with Dasatinib 5uM.
T394NS_dasatinib_5uM <- FindMarkers(
  object = pdgfrai_sobj_all,
  ident.1 = "T394NS_Dasatinib_5uM",
  ident.2 = "T394NS_DMSO_CTR")
T407NS_dasatinib_5uM <- FindMarkers(
  object = pdgfrai_sobj_all,
  ident.1 = "T407NS_Dasatinib_5uM",
  ident.2 = "T407NS_DMSO_CTR")

T394NS_dasatinib_5uM_hits <- T394NS_dasatinib_5uM %>% 
  filter(avg_log2FC > 0)
colnames(T394NS_dasatinib_5uM_hits) <- paste0("T394", colnames(T394NS_dasatinib_5uM_hits))
T394NS_dasatinib_5uM_hits$gene_name <- rownames(T394NS_dasatinib_5uM_hits)

T407NS_dasatinib_5uM_hits <- T407NS_dasatinib_5uM %>% 
  filter(avg_log2FC > 0)
colnames(T407NS_dasatinib_5uM_hits) <- paste0("T407", colnames(T407NS_dasatinib_5uM_hits))
T407NS_dasatinib_5uM_hits$gene_name <- rownames(T407NS_dasatinib_5uM_hits)

# Are there common genes for Dasatinib 5uM? PDEs are high on the list
dasatinib_hits <- T394NS_dasatinib_5uM_hits %>% 
  inner_join(T407NS_dasatinib_5uM_hits, by="gene_name")

# Common genes across both drugs? ~25%
sum(dasatinib_hits$gene_name%in%cp_hits$gene_name)/length(cp_hits$gene_name)
# Several relevant genes including a 5 PDEs.
dasatinib_hits$gene_name[dasatinib_hits$gene_name%in%cp_hits$gene_name]

# Importantly, are there are any major differences across the cell lines? Looks like PDGFRA is higher in T407NS.
cell_line_dmso_deg <- FindMarkers(
  object = pdgfrai_sobj_all,
  ident.1 = "T394NS_DMSO_CTR",
  ident.2 = "T407NS_DMSO_CTR")
cell_line_dmso_deg$gene_symbol <- rownames(cell_line_dmso_deg)

cell_line_dmso_deg %>% 
  filter(gene_symbol=="PDGFRA")

cell_line_das_deg <- FindMarkers(
  object = pdgfrai_sobj_all,
  ident.1 = "T394NS_Dasatinib_1uM",
  ident.2 = "T407NS_Dasatinib_1uM")
cell_line_das_deg$gene_symbol <- rownames(cell_line_das_deg)

cell_line_cp_deg <- FindMarkers(
  object = pdgfrai_sobj_all,
  ident.1 = "T394NS_CP-673451_1uM",
  ident.2 = "T407NS_CP-673451_1uM")
cell_line_cp_deg$gene_symbol <- rownames(cell_line_cp_deg)


### ### ### ### ### ### ### ### ### ### ### ###
### Signature scoring
### ### ### ### ### ### ### ### ### ### ### ###
# Read in the signatures for human tumor IDH-mutant malignant metaprograms
mut_mp <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/metaprograms/care_mut_selected_malignant_metaprograms.csv", sep = ",", header = T, row.names = 1)
mut_mp_list <- lapply(names(mut_mp), function(col_name) mut_mp[[col_name]])
names(mut_mp_list) <-paste0(colnames(mut_mp), "_MUT")


# Include public signatures in order to score the IDH-mutant hierarchy used elsewhere in the manuscript - if needed.

######### Venteicher IDH-A  #########
venteicher <- readWorkbook("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/misc/venteicher_table_s3.xlsx", startRow = 5, colNames = TRUE)
venteicher_signatures <- venteicher %>% 
  dplyr::select("Venteicher_OC_scRNA"=`Oligo-program.(Fig..2C)`, "Venteicher_AC_scRNA"=`Astro-program.(Fig..2C)`,
                "Venteicher_Stemness_scRNA"=`Stemness.program.(Fig..3C)`)
venteicher_signatures$Venteicher_OC_scRNA <- trimws(venteicher_signatures$Venteicher_OC_scRNA)
venteicher_signatures$Venteicher_AC_scRNA <- trimws(venteicher_signatures$Venteicher_AC_scRNA)
venteicher_signatures$Venteicher_Stemness_scRNA <- trimws(venteicher_signatures$Venteicher_Stemness_scRNA)
venteicher_sig_list <- as.list(venteicher_signatures[1:50,1:3])
venteicher_sig_list <- lapply(venteicher_sig_list, function(x) x[!is.na(x)])

######### Tirosh IDH-O    #########
tirosh <- readWorkbook("//vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/misc/tirosh_nature_2016_supplementary_table_1.xlsx", startRow = 9, colNames = TRUE)
tirosh_signatures <- tirosh %>% 
  dplyr::select("Tirosh_OC_scRNA"=`OC.(PCA-only)`,
                "Tirosh_AC_scRNA"=`AC.(PCA-only)`,
                "Tirosh_Stemness_scRNA"=`stemness`)
tirosh_signatures$Tirosh_OC_scRNA <- trimws(tirosh_signatures$Tirosh_OC_scRNA)
tirosh_signatures$Tirosh_AC_scRNA <- trimws(tirosh_signatures$Tirosh_AC_scRNA)
tirosh_signatures$Tirosh_Stemness_scRNA <- trimws(tirosh_signatures$Tirosh_Stemness_scRNA)
tirosh_signatures_list <- as.list(tirosh_signatures[1:50,1:3])

sigs_list <- c(mut_mp_list, venteicher_sig_list, tirosh_signatures_list)

### ### ### ### ### ### ### ### ### ### ###
## Create a umi list for analysis and name the elements
### ### ### ### ### ### ### ### ### ### ###
umi_data_list <- setNames(lapply(SAMPLES_NAMES, function(sample_name){
  print(sample_name)
  cur_data <- Read10X(paste0(in_data_dir, "/", sample_name,'/outs/filtered_feature_bc_matrix/'), gene.column=2)
  colnames(cur_data) <- paste0(sample_name, "-", colnames(cur_data))
  return(cur_data)
}), SAMPLES_NAMES)


#########################################################################################
##### Assign within-sample metaprogram score - 
#########################################################################################
# Define a metadata table for which we input into the score within samples function.
pdgfrai_md <- pdgfrai_sobj_all@meta.data
table(pdgfrai_md$SampleID)

# Define a function that subsets columns based on CellID
subset_columns <- function(mat, cell_ids) {
  col_idx <- which(colnames(mat) %in% cell_ids)
  mat[, col_idx, drop = FALSE]
}

# Apply the function to each element of the list
umi_data_all <- map(umi_data_list, subset_columns, cell_ids = pdgfrai_md$CellID)
names(umi_data_all)
lapply(umi_data_all, dim)

# Perform the scoring within a single sample. Note that AddModuleScore results differ depending on cell/sample dataset.
# We are only scoring those cells that are included post-processing.
set.seed(123)
mp_scores <- score_within_samples_caremut(umi_data_all, md = pdgfrai_md, sigs = sigs_list)
# saveRDS(mp_scores, paste0(out_data_dir, "pdgfrai_mp_sig_scores_20260406.RDS"))
# mp_scores <- readRDS(paste0(out_data_dir, "pdgfrai_mp_sig_scores_20260406.RDS"))

# Previous run
mp_scores_prior <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/perturbation/pdgfrai/pdgfrai_mp_sig_scores_20260406.RDS")
identical(mp_scores, mp_scores_prior[ ,1:14])

#####################################################################
### Cell state score assignment
#####################################################################
# mdata is an object (tibble) that contains the meta-data for the cells that should be classified (Sample, Treatment, etc.).
# The object also must include the CellID variable that identifies the cells.
# MP_scores (tibble) includes the scores for the meta-program and the CellID variable. All meta-programs score columns should
# start with the literal "MP_".
# Due to the shuffling each variable in the tibble is theoretically normally distributed (or is at least close to being ND).
mdata <- pdgfrai_md %>%
  left_join(mp_scores,
            by = c("CellID", "SampleID")) %>% 
  as_tibble() 

## The approach requires variables to be named in a certain way:
sigs <- mut_mp_list

# We call this function to generate a NULL distribution to facilitate classification from this experimental dataset. According to the configured parameters
# it will sample 5000 cells from the pool of cells, shuffle the expression values while maintaining the mean expression of
# each gene and score the artificial cells for the meta-programs. It will repeat the process 20 times to generate a NULL 
# distribution of 100K cells. It returns a tibble of 100K x n (where n is the number of meta-programs)
set.seed(42)
permuted_data <- generate_null_dist(umi_data_list = umi_data_all,
                                    md = mdata,
                                    sigs = sigs,
                                    n_iter = 20, n_cells = 5000, verbose = T)

# Save output of the permuted data
saveRDS(permuted_data, paste0(out_data_dir, "pdgfrai_class_permuted_data_20260406.RDS"))
# permuted_data <- readRDS( paste0(out_data_dir, "pdgfrai_class_permuted_data_20260406.RDS"))

# Previous run
permuted_data_prior <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/perturbation/pdgfrai/pdgfrai_class_permuted_data_20260406.RDS")
identical(permuted_data, permuted_data_prior)

# The signatures are the cell state programs
state_programs <- sigs

vars <- names(state_programs)

# We fit a normal distribution to each of the variables of the permuted data
scores_nd <- lapply(colnames(permuted_data), function(mp) {
  
  x <- permuted_data %>% pull(mp)
  
  fit <- MASS::fitdistr(x, "normal")
  class(fit)
  
  para <- fit$estimate
  
  tibble(MP = mp, Mean = para[1], SD = para[2])  
})
scores_nd <- do.call(rbind, scores_nd)

mean_vec <- setNames(scores_nd$Mean, scores_nd$MP)
sd_vec <- setNames(scores_nd$SD, scores_nd$MP)

####################################################################################################################################
# Plot the actual scores vs. the NULL distribution for each MP
####################################################################################################################################
library(reshape2)
library(ggdist)
norm_fit <- lapply(colnames(permuted_data), function(mp) tibble(MP = mp,
                                                                Sig = rnorm(n = 100000,
                                                                            mean = mean_vec[mp],
                                                                            sd = sd_vec[mp])))
norm_fit <- do.call(rbind, norm_fit)

dm <- rbind(melt(data = permuted_data,
                 measure.vars = vars) %>% mutate(DataType = "Permuted"),
            melt(data = mdata %>% select(ends_with("_MUT")),
                 measure.vars = vars) %>% mutate(DataType = "Actual"),
            norm_fit %>%
              rename(variable = MP, value = Sig) %>%
              mutate(DataType = "Classifier")) %>%
  mutate(DataType = factor(DataType, c("Permuted", "Actual", "Classifier")))

dm_stats <- dm %>%
  group_by(variable) %>%
  filter(DataType == "Classifier") %>%
  summarise(Q95 = quantile(value, .95), Q99 = quantile(value, .99))

png(paste0(fig_dir, "pdgfrai_mp_classification_permuted_actual.png"), width = 12, height = 8, units = 'in', res = 300, bg = "transparent")
ggplot(data = dm, aes(x = value, y = after_stat(ncount), color = DataType, linetype = DataType)) +
  facet_wrap(~variable, scales = "free_x", nrow = 2) +
  geom_freqpoly(bins = 100, size = 1, show.legend = c(color = T, linetype = F)) +
  scale_color_manual(name = "", values = c("Permuted" = "dodgerblue", "Actual" = "red", "Classifier" = "black")) +
  scale_linetype_manual(values = c("Permuted" = "solid", "Actual" = "solid", "Classifier" = "dashed")) +
  scale_fill_discrete(name = "Data distribution") +
  geom_vline(data = dm_stats, mapping = aes(xintercept = Q99), linetype = "dashed", size = 1) +
  xlab("Program score") +
  ylab("Count (scaled to 1)") + 
  plot_theme
dev.off()

####################################################################################################################################
# Cell state classification
####################################################################################################################################
# Names of all MPs to-be-classified should be included in this vector. 
state_vars <- names(state_programs)

# Melt the data to make the computation easier
state_data <- melt(data = mdata %>%
                     select(CellID, all_of(state_vars)),
                   id.vars = "CellID",
                   variable.name = "Program",
                   value.name = "Score",
                   measure.vars = state_vars)
state_data <- as_tibble(state_data)
state_data$Program <- as.character(state_data$Program)
table(state_data$Program)

# Generate a Z-score for each MP (using the mean and SD of the NULL distribution)
state_data <- state_data %>%
  mutate(Score_z = (Score - mean_vec[Program]) / sd_vec[Program])

# Compute a p-value for each (CellID, MP) pair. This is a one-sided test with the hypothesis that
# the actual score is not greater than expected by chance
state_data <- state_data %>%
  mutate(p.val = pnorm(Score_z, lower.tail = F))

# Correct for multiple testing **within each cell** using the Holm method
state_data <- state_data %>%
  group_by(CellID) %>%
  mutate(p.adj = p.adjust(p.val, "holm"))

# Nominal p-value seems to yield more consistent results for single cell data when using single NUCLEUS metaprogram gene sets.
state_data$p.sig <- state_data$p.val < .05

# Compute the classification statistics for each MP/gene set 
state_stats <- state_data %>%
  group_by(Program) %>%
  summarise(n = sum(p.sig == T), N = n(), Freq = n / N)

# This is the actual classification step. We filter out the statistically insignificant scores and classify the cell
# to the MP with the maximal signal. Classify CC separate from the other states.
# The classification will do the following:
# 1. Do not consider cell cycle MP (at first)
# 2. Restrict only those programs with p-value < 0.05
# 3. Rank by descending score such multiple significant scores will be assigned to the one with the highest Score (AddModuleScore).
state_data_classify <- state_data %>% 
  filter(!Program%in%c("MP_CC_MUT")) %>% 
  group_by(CellID) %>%
  filter(p.sig == T) %>%
  arrange(desc(Score)) %>%
  filter(!duplicated(CellID)) %>%
  ungroup()

# Perform this step separately for Cell Cycle:
state_data_classify_cc <- state_data %>% 
  filter(Program=="MP_CC_MUT") %>% 
  group_by(CellID) %>%
  filter(p.sig == T) %>%
  arrange(desc(Score)) %>%
  filter(!duplicated(CellID)) %>%
  ungroup()


####################################################################################################################################
# Assign the state
####################################################################################################################################
# If you have a large proportion of cells that are "Undifferentiated" (i.e. they did not achieve a significant adjusted p-value for any MP)
# then you can adjust the multiple-adjustment method or classification threshold to less stringent values as we did above.
state_vec <- setNames(rep("Undifferentiated", nrow(mdata)), mdata$CellID)
state_vec[state_data_classify$CellID] <- state_data_classify$Program
table(state_vec)
table(state_vec) / length(state_vec)

mdata$State <- state_vec[mdata$CellID]
table(mdata$State)

# There does not appear to be any major differences across states.
png(paste0(fig_dir, "pdgfrai_caremut_mp_classification_complexity.png"), width = 8, height = 5, units = 'in', res = 300, bg = "transparent")
mdata %>% 
  mutate(State = factor(State, c(names(mut_mp_list), "Undifferentiated"))) %>%
  ggplot(aes(x = State, y = nFeature_RNA)) + 
  ggdist::stat_halfeye(adjust = .5, width = .75, justification = -.2, .width = 0, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.color = NA) +
  coord_cartesian(xlim = c(1.2, NA)) +
  xlab("") +
  ylab("nFeature_RNA") +
  geom_hline(yintercept = median(mdata$nFeature_RNA), linetype = "dashed", size = 1, color = "red") +
  scale_y_continuous(breaks = seq(0, 7000, 1000)) +
  theme_bw() +
  theme(panel.grid.major = element_line())
dev.off()

## Cell cycle is not considered a cellular state but rather a feature (since cells can have a clear identity such as OPC or NPC and still be cycling) 
mdata_out <- mdata %>% 
  mutate(isCC = ifelse(mdata$CellID%in%state_data_classify_cc$CellID, TRUE, FALSE))

# Write out scores and assignment so that it is easier to read back in and reproduce figures.
write.table(mdata_out, file = paste0(out_data_dir, "pdgfrai_caremut_mp_select_state_assignment_20260406.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# Create an easy linker file to add back after calculating cell state frequencies.
linker_info <- mdata_out %>% 
  dplyr::select(SampleID, exp_group:dose) %>% 
  distinct()

# Recode and collapse malignant states from AC1 and AC2.
malignant_pval_freq <- mdata_out %>% 
  mutate(State = recode(State, `MP_OPC_MUT` = "OPC-like",
                        `MP_NPC_MUT` = "NPC-like",
                        `MP_AC1_MUT` = "AC-like", 
                        `MP_AC2_MUT` = "AC-like",
                        `MP_MES_MUT` = "MES-like",
                        `Undifferentiated` = "Undifferentiated")) %>% 
  group_by(SampleID, State) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(SampleID, State,
           fill = list(counts = 0, freq = 0)) %>%
  inner_join(linker_info, by="SampleID")

# Confirm that all sample sum 1 and that all cell types are measured.
malignant_pval_freq %>%
  group_by(SampleID) %>%
  summarise(freq_sum = sum(freq)) %>%
  mutate(pass = abs(freq_sum - 1) < 1e-10) %>%
  { if (all(.$pass)) {
    message("PASS: All samples frequencies sum to 1")
  } else {
    failing <- filter(., !pass)
    message("FAIL: ", nrow(failing), " samples do not sum to 1:")
    print(failing)
  }
  }
table(malignant_pval_freq$State)

# Independently calculate the cell cycle.
pdgfrai_state_freq_cc <- mdata_out %>% 
  group_by(SampleID, isCC) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(SampleID, isCC,
           fill = list(counts = 0, freq = 0)) %>%
  filter(isCC=="TRUE") %>% 
  dplyr::select(-isCC) %>% 
  mutate(State = "Cycling") %>% 
  inner_join(linker_info, by="SampleID") 

malignant_pval_freq <- malignant_pval_freq %>% 
  bind_rows(pdgfrai_state_freq_cc)

malignant_pval_freq$State <- factor(malignant_pval_freq$State, levels = c("MES-like", "AC-like", "OPC-like", "NPC-like", "Undifferentiated", "Cycling"))

pdf(paste0(fig_dir, "pdgfrai_malignant_cells_pval_score_assignment.pdf"), width = 4, height = 3.5, bg = "transparent")
malignant_pval_freq %>% 
  # Omitted the cycling cells so that all frequencies sum to 100%. As a reminder the cycling cells are separately enumerated.
  filter(State!="Cycling") %>% 
  ggplot(aes(x=exp_group, fill = factor(State), y=freq*100)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "Relative cell abundance (%)", fill="Cell State", x = "") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  facet_grid(.~cell_line, scales = "free_x", space = "free_x") +
  plot_theme +
  theme(panel.grid.major = element_line(),
        axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

# There is not a clear reduction in proliferation due to PDGFRAi, perhaps due to the FACS enrichment.
pdf(paste0(fig_dir, "pdgfrai_malignant_cells_pval_cc_assignment.pdf"), width = 8, height = 6, bg = "transparent")
pdgfrai_state_freq_cc %>% 
  ggplot(aes(x=exp_group, fill = factor(State), y=freq*100)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "Relative cell abundance (%)", fill="Is Cycling?", x = "") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6")) +
  facet_grid(.~cell_line, scales = "free_x", space = "free_x") +
  plot_theme +
  theme(panel.grid.major = element_line(),
        axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

malignant_pval_freq$treatment <- factor(malignant_pval_freq$treatment, levels = c("DMSO", "Dasatinib", "CP-673451"))

pdf(paste0(fig_dir, "pdgfrai_malignant_cells_pval_state_differences.pdf"), width = 10, height = 6, bg = "transparent")
malignant_pval_freq %>% 
  ggplot(aes(x=treatment, fill = factor(State), y=freq*100)) +
  geom_boxplot() +
  geom_point(aes(shape=dose)) +
  labs(y = "Relative cell abundance (%)", fill="Cell State", x = "") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90",
                             "Cycling" = "#6BAED6")) +
  facet_grid(cell_line~State, scales = "free_x", space = "free_x") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

pdf(paste0(fig_dir, "pdgfrai_malignant_cells_pval_state_differences_select.pdf"), width = 8, height = 6, bg = "transparent")
malignant_pval_freq %>% 
  filter(State%in%c("Undifferentiated", "AC-like", "MES-like")) %>% 
  ggplot(aes(x=treatment, fill = factor(State), y=freq*100)) +
  geom_boxplot() +
  geom_point(aes(shape=dose)) +
  labs(y = "Relative cell abundance (%)", fill="Cell State", x = "") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90",
                             "Cycling" = "#6BAED6")) +
  facet_grid(cell_line~State, scales = "free_x", space = "free_x") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

# Underpowered with this many samples for statistical assessment.
malignant_pval_freq %>% 
  mutate(treatment_binary = ifelse(treatment=="DMSO", "DMSO", "PDGFRAi")) %>% 
  filter(State%in%c("Undifferentiated", "AC-like")) %>% 
  ggplot(aes(x=treatment_binary, y=freq*100)) +
  geom_boxplot(aes(fill = factor(State))) +
  geom_point(aes(shape=dose)) +
  labs(y = "Relative cell abundance (%)", fill="Cell State", x = "") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90",
                             "Cycling" = "#6BAED6")) +
  facet_grid(.~State, scales = "free_x", space = "free_x") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  stat_compare_means(method="t.test")


### ### ### ### ### ### ### ###
### 2D graph for directionality of changes
### ### ### ### ### ### ### ###
malignant_pval_freq_wide <- malignant_pval_freq %>%
  mutate(freq = freq*100) %>% 
  dplyr::select(SampleID, State, freq) %>% 
  pivot_wider(names_from = State,  values_from = freq) %>% 
  mutate(Undiff_Stem = Undifferentiated+`NPC-like`,
         AC_MES = `AC-like`+`MES-like`) %>% 
  inner_join(linker_info, by="SampleID") 


arrow_df <- malignant_pval_freq_wide %>%
  group_by(cell_line) %>%
  # Identify control row per cell line
  mutate(ctrl_x = AC_MES[which(treatment == "DMSO")],
         ctrl_y = Undiff_Stem[which(treatment == "DMSO")]) %>%
  ungroup() %>%
  filter(treatment != "DMSO") %>%
  mutate(x = ctrl_x, y = ctrl_y, xend = AC_MES, yend = Undiff_Stem)


pdf(paste0(fig_dir, "pdgfrai_stem_acmes_arrows.pdf"), width = 8, height = 4, bg = "transparent", useDingbats = FALSE)
ggplot(malignant_pval_freq_wide, aes(x = AC_MES, y = Undiff_Stem)) +
  geom_point(aes(shape = treatment, color = dose), size = 3) +
  geom_segment(data = arrow_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(type = "closed", length = unit(0.4, "cm")),
               linetype = "dotted",
               size = 0.5,
               color = "black") +
  facet_grid(.~cell_line, space = "free") +
  scale_color_manual(values=c("CTR" = "#2166ac",
                              "1uM"="#fddbc7",
                              "5uM"="#b2182b")) +
  labs(x = "AC/MES-like (%)", y = "Undiff/Stem-like (%)", shape = "Treatment", color = "Dose") +
  plot_theme
dev.off()


arrow_df <- malignant_pval_freq_wide %>%
  group_by(cell_line) %>%
  # Identify control row per cell line
  mutate(ctrl_x = AC_MES[which(treatment == "DMSO")],
         ctrl_y = Undifferentiated[which(treatment == "DMSO")]) %>%
  ungroup() %>%
  filter(treatment != "DMSO") %>%
  mutate(x = ctrl_x, y = ctrl_y, xend = AC_MES, yend = Undifferentiated)

pdf(paste0(fig_dir, "pdgfrai_undifferentiated_acmes_arrows.pdf"), width = 8, height = 4, bg = "transparent", useDingbats = FALSE)
ggplot(malignant_pval_freq_wide %>% 
         mutate(text_to_add = paste(treatment, " ", dose)), aes(x = AC_MES, y = Undifferentiated)) +
  geom_point(aes(color = treatment), size = 2) +
  geom_segment(data = arrow_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 25),
               linetype = "dashed",
               size = 0.25,
               color = "black") +
  geom_text(
    aes(label = text_to_add),
    nudge_x = 2,    
    nudge_y = 0.75,    
    size = 3
  ) +
  facet_grid(.~cell_line, space = "free") +
  scale_color_manual(values=c("DMSO" = "#2166ac",
                              "Dasatinib"="#fddbc7",
                              "CP-673451"="#b2182b")) +
  labs(x = "Astrocyte lineage (AC/MES-like %)", y = "Undifferentiated %", shape = "Group", color = "Condition") +
  plot_theme +
  guides(color = FALSE)
dev.off()

### ### ### ### ### ### ### ### ### ###
### Relative changes - easier interpretation, but larger figure footprint
### ### ### ### ### ### ### ### ### ###
malignant_relative <- malignant_pval_freq %>%
  group_by(cell_line, State) %>%
  mutate(freq_dmsomean = mean(freq[treatment == "DMSO"], na.rm = TRUE),
         freq_relative = (freq - freq_dmsomean) / freq_dmsomean) %>%
  ungroup()

malignant_absolute <- malignant_pval_freq %>%
  group_by(cell_line, State) %>%
  mutate(freq_dmso = freq[treatment == "DMSO"],
         freq_abs = (freq - freq_dmso)*100) %>%
  ungroup()

ggplot(malignant_absolute %>% 
         filter(treatment!="DMSO") %>% 
         filter(State%in%c("Undifferentiated", "AC-like", "MES-like")), aes(x = treatment, fill = factor(State), y = freq_abs)) +
  geom_boxplot() +
  geom_point(aes(shape = dose)) +
  labs(y = "Relative change vs DMSO (%)", fill = "Cell State", x = "") +
  scale_fill_manual(values = c("AC-like" = "#AA2756", 
                               "MES-like" = "#F77D58",
                               "NPC-like" = "#7fbf7b",
                               "OPC-like" = "#E8F5A3",
                               "Undifferentiated" = "gray90",
                               "Cycling" = "#6BAED6")) +
  facet_grid(cell_line ~ State, scales = "free_x", space = "free_x") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0)

ggplot(malignant_absolute %>% 
         filter(treatment!="DMSO") %>% 
         filter(State%in%c("Undifferentiated", "AC-like", "MES-like")), aes(x = treatment, y = freq_abs, fill = dose)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(State ~ cell_line, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("1uM" = "#a6cee3",  # light blue
                               "5uM" = "#1f78b4")) +
  labs(y = "Change in cell abundance vs DMSO (%)", fill = "Dose", x = "Treatment") +
  plot_theme +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(malignant_absolute %>% 
         filter(treatment!="DMSO") %>% 
         filter(State%in%c("Undifferentiated", "AC-like", "MES-like")), aes(x = treatment, y = freq_abs, fill = dose)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(cell_line ~ State, scales = "fixed", space = "free_x") +
  scale_fill_manual(values = c("1uM" = "#a6cee3",  # light blue
                               "5uM" = "#1f78b4")) +
  labs(y = "Change in cell abundance vs DMSO (%)", fill = "Dose", x = "Treatment") +
  plot_theme +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(fig_dir, "pdgfrai_relative_state_change_dasatinib.pdf"), width = 8, height = 6, bg = "transparent")
ggplot(malignant_absolute %>% 
         filter(treatment=="Dasatinib") %>% 
         filter(State%in%c("Undifferentiated", "AC-like")), aes(x = treatment, y = freq_abs, fill = dose)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(cell_line ~ State, scales = "fixed", space = "free_x") +
  scale_fill_manual(values = c("1uM" = "#a6cee3",  # light blue
                               "5uM" = "#1f78b4")) +
  labs(y = "Change in cell abundance vs DMSO (%)", fill = "Dose", x = "Treatment") +
  plot_theme +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf(paste0(fig_dir, "pdgfrai_relative_state_change_cp-673451.pdf"), width = 8, height = 6, bg = "transparent")
ggplot(malignant_absolute %>% 
         filter(treatment=="CP-673451") %>% 
         filter(State%in%c("Undifferentiated", "AC-like", "MES-like")), aes(x = treatment, y = freq_abs, fill = dose)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(cell_line ~ State, scales = "fixed", space = "free_x") +
  scale_fill_manual(values = c("1uM" = "#a6cee3",  # light blue
                               "5uM" = "#1f78b4")) +
  labs(y = "Change in cell abundance vs DMSO (%)", fill = "Dose", x = "Treatment") +
  plot_theme +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

malignant_combined <- malignant_absolute %>%
  filter(treatment == "CP-673451") %>%
  filter(State %in% c("Undifferentiated", "AC-like", "MES-like")) %>%
  mutate(New_State = ifelse(State %in% c("AC-like", "MES-like"), "AC/MES-like", "Undifferentiated")) %>%
  group_by(cell_line, dose, treatment, New_State) %>%
  summarise(freq_abs = sum(freq_abs)) %>% 
  ungroup()

ggplot(malignant_combined, aes(x = treatment, y = freq_abs, fill = dose)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(cell_line ~ New_State, scales = "fixed", space = "free_x") +
  scale_fill_manual(values = c("1uM" = "#a6cee3",  # light blue
                               "5uM" = "#1f78b4")) +
  labs(y = "Change in cell abundance vs DMSO (%)", fill = "Dose", x = "Treatment") +
  plot_theme +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### END ###