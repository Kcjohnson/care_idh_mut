##################################
# 10x scRNA processing of MGG152 IDH-mutant gliomaspheres with CRISPR targeting CDKN2A from the Suva/Cahill labs
# Author: Kevin Johnson
# Date Updated: 2026.04.06
##################################

library(tidyverse) # 1.3.1 
library(Seurat) # 4.3.0
library(Matrix) # 1.6-5
library(ggpubr) # 0.4.0 
library(DoubletFinder) # 2.0.3
library(harmony) # 1.2.0 
library(openxlsx) # 4.2.5.2 

# Notes: These results were based on cellranger 9.0.1 to align with other in vitro data that required recent versions of cellranger 
# The Cell Ranger count method used the 10x Genomics 2020 reference transcriptome (refdata-gex-GRCh38-2020-A)
# MGG152 Parental, sgRNA1, sgRNA3, and non-targeting sgRNA.

# Define where figures and results will be stored:
fig_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/figures/perturbation/"
out_data_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/perturbation/cdkn2a/"
setwd("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/")
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")

# Specify directories that have the count matrices. 
in_data_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/cellranger/cdkn2a"

# List all samples in each directory. 
SAMPLES <- dir(in_data_dir, full.names = T)

# Create function to grab the last element of each vector (varying length vectors) from file path.
last <- function(x) { return( x[length(x)] ) }
SAMPLES_NAMES <- sapply(strsplit(SAMPLES, "/"), last)

seurat_list <- lapply(SAMPLES_NAMES, function(sample_name){
  cur_data <- Read10X(paste0(in_data_dir, "/", sample_name,'/outs/filtered_feature_bc_matrix/'))
  
  # Adding the sample identifier to the cell barcode.
  colnames(cur_data) <- paste0(sample_name, "-", colnames(cur_data))
  
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3,
    min.features=200,
    project='cdkn2a'
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

pdf(paste0(fig_dir, "preqc_mito_violin_cdkn2a.pdf"), width=10, height=6)
VlnPlot(seurat_obj, group.by="SampleID", features = "percent.mt", ncol = 1, pt.size=0) + NoLegend() + ylim(0, 10)
dev.off()

pdf(paste0(fig_dir, "preqc_nfeature_rna_violin_cdkn2a.pdf"), width=10, height=6)
VlnPlot(seurat_obj, group.by="SampleID", features = "nFeature_RNA", ncol = 1, pt.size=0) + NoLegend()
dev.off()

# In the past, we set a relatively liberal criteria in tumor snRNAseq data of 500 genes detected so as not to miss lymphocytes, which do not express many genes.
# These data were deeply sequenced and come from single cell cell line models. 
# Based on the distribution above, we can set a higher threshold to select the best cells.
sum(seurat_obj$nFeature_RNA>1000 & seurat_obj$nFeature_RNA<2500) 
sum(seurat_obj$nFeature_RNA<2500) 
sum(seurat_obj$nFeature_RNA>10000) 
hist(seurat_obj$percent.mt)
abline(v = 7.5, col = "red", lwd = 2)
sum(seurat_obj$percent.mt>7.5) # Inspect how many cells would we lose depending on percent.mt threshold. 5-20% may all be valid.

# There is a second peak in these data with lower gene counts. Filtering out that population with more stringent threshold.
hist(seurat_obj$nFeature_RNA)
abline(v = 2500, col = "red", lwd = 2)

# Assess the changes in cell number following QC.
pre_qc_cell_num <- dim(seurat_obj)[2]

# Limiting to fewer than 10K genes to remove doublets plus setting a threshold for mitochondrial genes and a minimum number of genes.
# Note that these values differ from other experiments due to high data quality.
# Here, when initially processing the data, there seemed to be co-clustering across experiments for cells with lower gene counts.
seurat_obj <- subset(seurat_obj, nFeature_RNA > 2500 & nFeature_RNA < 10000 & percent.mt < 7.5)

# 19,740 cells were removed. 18,741 cells left over.
dropped_cells <- pre_qc_cell_num-dim(seurat_obj)[2]
sprintf("%s cells were removed due to low quality", dropped_cells)

pdf(paste0(fig_dir, "postqc_nfeature_rna_violin_cdkn2a.pdf"), width=10, height=6)
VlnPlot(seurat_obj, group.by="SampleID", features = "nFeature_RNA", ncol = 1, pt.size=0) + NoLegend()
dev.off()

# The distribution is such that many observations sit between 4-7.5%.
pdf(paste0(fig_dir, "postqc_mito_violin_cdkn2a.pdf"), width=10, height=6)
VlnPlot(seurat_obj, group.by="SampleID", features = "percent.mt", ncol = 1, pt.size=0) + NoLegend()
dev.off()

#########################################################################################
##### Doublet Finder
#########################################################################################

# I followed the following example: https://rpubs.com/kenneditodd/doublet_finder_example

# Split aggregated data by sample
cdkn2a_split <- SplitObject(seurat_obj, split.by = "SampleID") 

# loop through samples to find doublets
for (i in 1:length(cdkn2a_split)) {
  # print the sample we are on
  print(paste0("Sample ",i))
  
  # Pre-process seurat object with standard seurat workflow
  cdkn2a_sample <- NormalizeData(cdkn2a_split[[i]])
  cdkn2a_sample <- FindVariableFeatures(cdkn2a_sample)
  cdkn2a_sample <- ScaleData(cdkn2a_sample)
  cdkn2a_sample <- RunPCA(cdkn2a_sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- cdkn2a_sample[["pca"]]@stdev
  sum.stdv <- sum(cdkn2a_sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # finish pre-processing
  cdkn2a_sample <- RunUMAP(cdkn2a_sample, dims = 1:min.pc)
  cdkn2a_sample <- FindNeighbors(object = cdkn2a_sample, dims = 1:min.pc)              
  cdkn2a_sample <- FindClusters(object = cdkn2a_sample, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep_v3(cdkn2a_sample)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bimodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- cdkn2a_sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(cdkn2a_sample@meta.data)) ## Assuming 7.5% doublet formation rate
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  cdkn2a_sample <- doubletFinder_v3(seu = cdkn2a_sample, 
                                     PCs = 1:min.pc, 
                                     pK = optimal.pk,
                                     nExp = nExp.poi.adj)
  metadata <- cdkn2a_sample@meta.data
  # Get the number of columns in the data frame
  num_columns <- ncol(metadata)
  
  # Change the name of the last column
  colnames(metadata)[num_columns] <- "doublet_finder"
  
  cdkn2a_sample@meta.data <- metadata 
  
  # subset and save
  cdkn2a_singlets <- subset(cdkn2a_sample, doublet_finder == "Singlet")
  cdkn2a_split[[i]] <- cdkn2a_singlets
  remove(cdkn2a_singlets)
}

print("Finished detecting singlets")


# Converge the CDKN2A splits after identifying the singlets.
seurat_obj_singlets <- merge(x=cdkn2a_split[[1]], 
                             y=cdkn2a_split[2:length(cdkn2a_split)],
                             project = "cdkn2a")

columns_to_remove <- grep("^pANN", colnames(seurat_obj_singlets@meta.data), value = TRUE)
seurat_obj_singlets@meta.data <- seurat_obj_singlets@meta.data[, !colnames(seurat_obj_singlets@meta.data) %in% columns_to_remove]

## Save the files
saveRDS(seurat_obj_singlets, paste0(out_data_dir, "cdkn2a_seurat_obj_singlets_20250617_min2500genes.RDS"))

# seurat_obj_singlets <- readRDS(paste0(out_data_dir, "cdkn2a_seurat_obj_singlets_20250617_min2500genes.RDS"))

# A few additional plots to visualize post-QC steps
seurat_obj_singlets@meta.data$plot_id <- seurat_obj_singlets@meta.data$SampleID
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="gRNA_new"] <- "sgNTC" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="gRNA1"] <- "sgRNA1" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="gRNA3"] <- "sgRNA3" 
seurat_obj_singlets@meta.data$plot_id <- factor(seurat_obj_singlets@meta.data$plot_id, levels= c("Parental", "sgNTC", "sgRNA1", "sgRNA3"))

pdf(paste0(fig_dir, "pass_qc_nfeature_rna_violin_cdkn2a.pdf"), width=3.5, height=3, useDingbats = FALSE)
VlnPlot(seurat_obj_singlets, group.by="plot_id", features = "nFeature_RNA", ncol = 1, pt.size=0) + 
  geom_hline(yintercept = c(2500, 10000), linetype = "dashed", color="red") +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    color = "black",
    fill = "white"
  )  +
  NoLegend() + ylim(0, 10000) + xlab("Sample ID") + ylab("nFeature RNA") + theme(plot.title = element_blank()) +
  plot_theme
dev.off()

pdf(paste0(fig_dir, "pass_qc_mito_rna_violin_cdkn2a.pdf"), width=3.5, height=3, useDingbats = FALSE)
VlnPlot(seurat_obj_singlets, group.by="plot_id", features = "percent.mt", ncol = 1, pt.size=0) + 
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

# Create an image for the total number of cells passing quality control for each model.
cluster_counts <- table(seurat_obj_singlets@meta.data$SampleID)
cluster_df <- as.data.frame(cluster_counts)
colnames(cluster_df) <- c("SampleID", "count")
cluster_df$SampleID <- as.character(cluster_df$SampleID)
cluster_df$SampleID[cluster_df$SampleID=="gRNA_new"] <- "sgNTC" 
cluster_df$SampleID[cluster_df$SampleID=="gRNA1"] <- "sgRNA1" 
cluster_df$SampleID[cluster_df$SampleID=="gRNA3"] <- "sgRNA3" 

pdf(paste0(fig_dir, "pass_qc_cell_counts_cdkn2a.pdf"), width=8, height=6)
ggplot(cluster_df, aes(x = SampleID, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(x = "SampleID", y = "Cells passing QC") +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#########################################################################################
##### Preprocess with standard processing  
#########################################################################################
# Run through a standard Seurat workflow
cdkn2a_sobj_all <- NormalizeData(seurat_obj_singlets, normalization.method = "LogNormalize", scale.factor = 10000)
cdkn2a_sobj_all <- FindVariableFeatures(cdkn2a_sobj_all, selection.method = "vst", nfeatures = 5000)
# By default only variable features are scaled
all.genes <- rownames(cdkn2a_sobj_all)
cdkn2a_sobj_all <- ScaleData(cdkn2a_sobj_all, features = all.genes)
cdkn2a_sobj_all <- RunPCA(cdkn2a_sobj_all, features = VariableFeatures(object = cdkn2a_sobj_all), npcs=50)
cdkn2a_sobj_all <- FindNeighbors(cdkn2a_sobj_all, dims = 1:20) 
cdkn2a_sobj_all <- FindClusters(cdkn2a_sobj_all, resolution = 0.6) 
cdkn2a_sobj_all <- RunUMAP(cdkn2a_sobj_all, dims = 1:20, min.dist = 0.2, n.neighbors = 40)

# Rename the SampleID to be more directly interpretable compared with the raw data.
cdkn2a_sobj_all@meta.data$exp_id <- cdkn2a_sobj_all@meta.data$SampleID
cdkn2a_sobj_all@meta.data$exp_id[cdkn2a_sobj_all@meta.data$exp_id=="gRNA_new"] <- "sgNTC"
cdkn2a_sobj_all@meta.data$exp_id[cdkn2a_sobj_all@meta.data$exp_id=="gRNA1"] <- "sgRNA1"
cdkn2a_sobj_all@meta.data$exp_id[cdkn2a_sobj_all@meta.data$exp_id=="gRNA3"] <- "sgRNA3"
cdkn2a_sobj_all@meta.data$exp_group <- ifelse(cdkn2a_sobj_all@meta.data$exp_id%in%c("sgRNA1", "sgRNA3"), "CDKN2A-/-", "Control")

# Data exploration for key features.
options(ggrepel.max.overlaps = Inf)
DimPlot(cdkn2a_sobj_all, label = TRUE, group.by = "RNA_snn_res.0.6",repel = TRUE) + NoLegend()

png(paste0(fig_dir, "cdkn2a_standard_preprocessing_sample_id.png"), width=10, height=8, res=300, units='in')
DimPlot(cdkn2a_sobj_all, label = TRUE, group.by = "exp_id", repel = TRUE)
dev.off()

umap_plot <- DimPlot(cdkn2a_sobj_all, label = TRUE, group.by = "exp_id", repel = TRUE, label.size = 3) +
  labs(title = "Conditions (n = 17K cells)") + 
  plot_theme

ggsave(paste0(fig_dir, "cdkn2a_standard_preprocessing_sample_id.pdf"), umap_plot, width = 4, height = 3, dpi = 300)

# Inspect the plot to see how the key features are distributed.
FeaturePlot(cdkn2a_sobj_all, features = c("CDKN2A",  "TOP2A",  "percent.mt", "nFeature_RNA"))
FeaturePlot(cdkn2a_sobj_all, features = c("CDK6",  "SOX11",  "SOX4"))

Idents(cdkn2a_sobj_all) <- "exp_id"
png(paste0(fig_dir, "mgg152_cdkn2a_expression_sampleid.png"), width=8, height=6, res=300, units='in')
VlnPlot(cdkn2a_sobj_all, features = c("CDKN2A"))
dev.off()

# Inspect the marker genes for the different experimental conditions.
sgRNA1_markers <- FindMarkers(object = cdkn2a_sobj_all, 
                       ident.1 = "sgRNA1",
                       ident.2 = "Parental",
                       logfc.threshold = 0.5,
                       min.pct = 0.25)
colnames(sgRNA1_markers) <- paste0("sg1_", colnames(sgRNA1_markers))
sgRNA1_markers$gene_name <- rownames(sgRNA1_markers)
rownames(sgRNA1_markers) <- NULL
sgRNA3_markers <- FindMarkers(object = cdkn2a_sobj_all, 
                              ident.1 = "sgRNA3",
                              ident.2 = "Parental",
                              logfc.threshold = 0.5,
                              min.pct = 0.25)
colnames(sgRNA3_markers) <- paste0("sg3_", colnames(sgRNA3_markers))
sgRNA3_markers$gene_name <- rownames(sgRNA3_markers)
rownames(sgRNA3_markers) <- NULL

# SOX4 (NPC-like marker) is upregulated in both sgRNAs targeting CDKN2A.
shared_up_markers <- sgRNA3_markers %>% 
  inner_join(sgRNA1_markers, by="gene_name") %>% 
  filter(sg3_avg_log2FC > 0, sg1_avg_log2FC > 0)

# CDKN2A is among the top downregulated genes.
shared_down_markers <- sgRNA3_markers %>% 
  inner_join(sgRNA1_markers, by="gene_name") %>% 
  filter(sg3_avg_log2FC < 0, sg1_avg_log2FC < 0)

# Comare with the non-targeting guide RNA
sgRNA1_vs_ntg_markers <- FindMarkers(object = cdkn2a_sobj_all, 
                              ident.1 = "sgRNA1",
                              ident.2 = "sgNTC",
                              logfc.threshold = 0.5,
                              min.pct = 0.25)
colnames(sgRNA1_vs_ntg_markers) <- paste0("sg1_vs_ntg_", colnames(sgRNA1_vs_ntg_markers))
sgRNA1_vs_ntg_markers$gene_name <- rownames(sgRNA1_vs_ntg_markers)

sgRNA3_vs_ntg_markers <- FindMarkers(object = cdkn2a_sobj_all, 
                                     ident.1 = "sgRNA3",
                                     ident.2 = "sgNTC",
                                     logfc.threshold = 0.5,
                                     min.pct = 0.25)
colnames(sgRNA3_vs_ntg_markers) <- paste0("sg3_vs_ntg_", colnames(sgRNA3_vs_ntg_markers))
sgRNA3_vs_ntg_markers$gene_name <- rownames(sgRNA3_vs_ntg_markers)

shared_up_cdkn2a_ntg_markers <- sgRNA3_vs_ntg_markers %>% 
  inner_join(sgRNA1_vs_ntg_markers, by="gene_name") %>% 
  filter(sg3_vs_ntg_avg_log2FC > 0, sg1_vs_ntg_avg_log2FC > 0)
shared_down_cdkn2a_ntg_markers <- sgRNA3_vs_ntg_markers %>% 
  inner_join(sgRNA1_vs_ntg_markers, by="gene_name") %>% 
  filter(sg3_vs_ntg_avg_log2FC < 0, sg1_vs_ntg_avg_log2FC < 0)

# Several genes were upregulated across both sgRNA vs parental and sgRNA vs NTC sgRNA.
shared_up_markers$gene_name[shared_up_markers$gene_name%in%shared_up_cdkn2a_ntg_markers$gene_name]
shared_down_markers$gene_name[shared_down_markers$gene_name%in%shared_down_cdkn2a_ntg_markers$gene_name]

# The parental line has upregulated ribosomal protein genes. sgNTC has an enrichment for SOX4 indicating some selection for this phenotype when going through the targeting process.
parental_vs_ntg_markers <- FindMarkers(object = cdkn2a_sobj_all, 
                                     ident.1 = "Parental",
                                     ident.2 = "sgNTC",
                                     logfc.threshold = 0.5,
                                     min.pct = 0.25)
colnames(parental_vs_ntg_markers) <- paste0("parental_vs_ntg_", colnames(parental_vs_ntg_markers))
parental_vs_ntg_markers$gene_name <- rownames(parental_vs_ntg_markers)


### ### ### ### ### ### ### ### ### ### ###
## Create a umi list for analysis and name the elements
### ### ### ### ### ### ### ### ### ### ###
umi_data_list <- setNames(lapply(SAMPLES_NAMES, function(sample_name){
  print(sample_name)
  cur_data <- Read10X(paste0(in_data_dir, "/", sample_name,'/outs/filtered_feature_bc_matrix/'), gene.column=2)
  colnames(cur_data) <- paste0(sample_name, "-", colnames(cur_data))
  return(cur_data)
}), SAMPLES_NAMES)

### ### ### ### ### ### ### ### ### ### ### ###
### Signature scoring
### ### ### ### ### ### ### ### ### ### ### ###
# Read in the signatures for in vivo IDH-mutant malignant metaprograms
mut_mp <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/metaprograms/care_mut_selected_malignant_metaprograms.csv", sep = ",", header = T, row.names = 1)
mut_mp_list <- lapply(names(mut_mp), function(col_name) mut_mp[[col_name]])
names(mut_mp_list) <-paste0(colnames(mut_mp), "_MUT")


# Include public signatures in order to score the IDH-mutant hierarchy used elsewhere in the manuscript

# Venteicher astrocytoma
venteicher <- readWorkbook("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/misc/venteicher_table_s3.xlsx", startRow = 5, colNames = TRUE)
venteicher_signatures <- venteicher %>% 
  dplyr::select("Venteicher_OC_scRNA"=`Oligo-program.(Fig..2C)`, "Venteicher_AC_scRNA"=`Astro-program.(Fig..2C)`,
                "Venteicher_Stemness_scRNA"=`Stemness.program.(Fig..3C)`)
venteicher_signatures$Venteicher_OC_scRNA <- trimws(venteicher_signatures$Venteicher_OC_scRNA)
venteicher_signatures$Venteicher_AC_scRNA <- trimws(venteicher_signatures$Venteicher_AC_scRNA)
venteicher_signatures$Venteicher_Stemness_scRNA <- trimws(venteicher_signatures$Venteicher_Stemness_scRNA)
venteicher_sig_list <- as.list(venteicher_signatures[1:50,1:3])
venteicher_sig_list <- lapply(venteicher_sig_list, function(x) x[!is.na(x)])

# Tirosh oligodendroglioma
tirosh <- readWorkbook("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/misc/tirosh_nature_2016_supplementary_table_1.xlsx", startRow = 9, colNames = TRUE)
tirosh_signatures <- tirosh %>% 
  dplyr::select("Tirosh_OC_scRNA"=`OC.(PCA-only)`,
                "Tirosh_AC_scRNA"=`AC.(PCA-only)`,
                "Tirosh_Stemness_scRNA"=`stemness`)
tirosh_signatures$Tirosh_OC_scRNA <- trimws(tirosh_signatures$Tirosh_OC_scRNA)
tirosh_signatures$Tirosh_AC_scRNA <- trimws(tirosh_signatures$Tirosh_AC_scRNA)
tirosh_signatures$Tirosh_Stemness_scRNA <- trimws(tirosh_signatures$Tirosh_Stemness_scRNA)
tirosh_signatures_list <- as.list(tirosh_signatures[1:50,1:3])

# Load the Neftel et al signatures for the different cell cycle phases
neftel_sig <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/misc/neftel_metamodule_genelists.csv", sep=",", header = TRUE)
neftel_signatures_list <- as.list(neftel_sig[1:50, ])
neftel_signatures_list <- lapply(neftel_signatures_list, function(x) x[!is.na(x)])

sigs_list <- c(mut_mp_list, venteicher_sig_list, tirosh_signatures_list, neftel_signatures_list)

#########################################################################################
##### Assign within-sample metaprogram score - 
#########################################################################################
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/caremut_utils.R")

# Define a metadata table for which we input into the score within samples function.
cdkn2a_md <- cdkn2a_sobj_all@meta.data
cdkn2a_md$CellID <- rownames(cdkn2a_md)

# Define a function that subsets columns based on CellID
subset_columns <- function(mat, cell_ids) {
  col_idx <- which(colnames(mat) %in% cell_ids)
  mat[, col_idx, drop = FALSE]
}

# Apply the function to each element of the list
umi_data_all <- map(umi_data_list, subset_columns, cell_ids = cdkn2a_md$CellID)
names(umi_data_all)
lapply(umi_data_all, dim)

# Perform the scoring within a single sample. Note that AddModuleScore results differ depending on cell/sample dataset.
# We are only scoring those cells that are included post-processing.
set.seed(123)
mp_scores <- score_within_samples_caremut(umi_data_all, md = cdkn2a_md, sigs = sigs_list)
saveRDS(mp_scores, paste0(out_data_dir, "cdkn2a_mp_sig_scores_20260406.RDS"))

# mp_scores <- readRDS(paste0(out_data_dir, "cdkn2a_mp_sig_scores_20260406.RDS"))

#####################################################################
### Cell state score assignment
#####################################################################
# mdata is an object (tibble) that contains the meta-data for the cells that should be classified (Sample, Treatment, etc.).
# The object also must include the CellID variable that identifies the cells.
# MP_scores (tibble) includes the scores for the meta-program and the CellID variable. All meta-programs score columns should
# start with the literal "MP_".
# Due to the shuffling each variable in the tibble is theoretically normally distributed (or is at least close to being ND).
mdata <- cdkn2a_md %>%
  left_join(mp_scores,
            by = c("CellID", "SampleID")) %>% 
  as_tibble() 

## The approach requires variables to be named in a certain way:
sigs <- mut_mp_list


# We call this function to generate a NULL distribution to facilitate classification. According to the configured parameters
# it will sample 5000 cells from the pool of cells, shuffle the expression values while maintaining the mean expression of
# each gene and score the artificial cells for the meta-programs. It will repeat the process 20 times to generate a NULL 
# distribution of 100K cells. It returns a tibble of 100K x n (where n is the number of meta-programs. It's 6 MPs here).
# For permuted data in CDKN2A experiments, use reference controls
mdata_controls <- mdata %>% 
  filter(plot_id%in%c("Parental", "sgNTC"))

set.seed(42)
permuted_data <- generate_null_dist(umi_data_list = umi_data_all,
                                    md = mdata_controls,
                                    sigs = sigs,
                                    n_iter = 20, n_cells = 5000, verbose = T)

# Save output of the permuted data
# saveRDS(permuted_data, paste0(out_data_dir, "cdkn2a_class_permuted_data_20260406.RDS"))
# permuted_data <- readRDS( paste0(out_data_dir, "cdkn2a_class_permuted_data_20260406.RDS"))

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

png(paste0(fig_dir, "cdkn2a_mp_classification_permuted_actual.png"), width = 12, height = 8, units = 'in', res = 300, bg = "transparent")
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
# Names of all MPs to-be-classified should be included in this vector (can be used to exclude MPs that reflect artifact/low quality etc.)
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
state_data$p.sig <- state_data$p.val < 0.05

# Compute the classification statistics for each MP/gene set 
state_stats <- state_data %>%
  group_by(Program) %>%
  summarise(n = sum(p.sig == T), N = n(), Freq = n / N)

# This is the actual classification step. We filter out the statistically insignificant scores and classify the cell
# to the MP with the maximal signal. Classify CC separate from the other states.
# The classification will do the following:
# 1. Not consider cell cycle MP (at first)
# 2. Restrict only those programs with nominal p-value < 0.05
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
state_vec <- setNames(rep("Undifferentiated", nrow(mdata)), mdata$CellID)
state_vec[state_data_classify$CellID] <- state_data_classify$Program
table(state_vec)
table(state_vec) / length(state_vec)

mdata$State <- state_vec[mdata$CellID]
table(mdata$State)

# Is there any major difference in gene counts across the different classifications? These data have very high coverage.
png(paste0(fig_dir, "cdkn2a_caremut_mp_classification_complexity.png"), width = 8, height = 5, units = 'in', res = 300, bg = "transparent")
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
  mutate(isCC = ifelse(mdata$CellID%in%state_data_classify_cc$CellID, TRUE, FALSE)) %>% 
  mutate(SampleID = recode(SampleID, `gRNA1` = "sgRNA1",
                           `gRNA3` = "sgRNA3",
                           `gRNA_new` = "sgNTC",
                           `Parental` = "Parental"))

# Write out scores and assignment so that it is easier to read back in.
write.table(mdata_out, file = paste0(out_data_dir, "cdkn2a_caremut_mp_select_state_assignment.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# There is a clear increase in stemness for the sgRNAs, including NTC with the Venteicher signature but not Tirosh signature.
ggplot(mdata_out, aes(x=SampleID, y=Venteicher_Stemness_scRNA)) +
  geom_boxplot() +
  stat_compare_means(method="kruskal") + 
  plot_theme

ggplot(mdata_out, aes(x=SampleID, y=Tirosh_Stemness_scRNA)) +
  geom_boxplot() +
  stat_compare_means(method="kruskal") + 
  plot_theme

# sgRNA1 had a larger effect in Venteicher while sgRNA3 had a higher signature in Tirosh 
mdata_out %>% 
  group_by(SampleID) %>% 
  summarise(venteicher_stemness = median(Venteicher_Stemness_scRNA),
            tirosh_stemness = median(Tirosh_Stemness_scRNA))

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
           fill = list(counts = 0, freq = 0)) 


# Independently calculate the cell cycle.
cdkn2a_state_freq_cc <- mdata_out %>% 
  group_by(SampleID, isCC) %>% 
  summarise(counts = n()) %>% 
  mutate(freq = counts / sum(counts)) %>% 
  ungroup() %>% 
  complete(SampleID, isCC,
           fill = list(counts = 0, freq = 0)) %>%
  filter(isCC=="TRUE") %>% 
  dplyr::select(-isCC) %>% 
  mutate(State="Cycling") 

malignant_pval_freq <- malignant_pval_freq %>% 
  bind_rows(cdkn2a_state_freq_cc)

malignant_pval_freq$State <- factor(malignant_pval_freq$State, levels = rev(c("MES-like", "AC-like", "OPC-like", "NPC-like", "Undifferentiated", "Cycling")))
malignant_pval_freq$exp_group <- ifelse(malignant_pval_freq$SampleID%in%c("sgNTC","Parental"), "Control", "sgRNA")

malignant_pval_freq_wide <- malignant_pval_freq %>%
  mutate(freq = freq*100) %>% 
  dplyr::select(SampleID, State, freq) %>% 
  pivot_wider(names_from = State,  values_from = freq) %>% 
  mutate(Undiff_Stem = Undifferentiated+`NPC-like`,
         AC_MES = `MES-like`+`AC-like`) %>% 
  mutate(exp_group = ifelse(SampleID%in%c("sgNTC","Parental"), "Control", "CDKN2A-/-"))


malignant_pval_freq_wide <- malignant_pval_freq_wide %>% 
  mutate(exp_type = recode(SampleID, `sgRNA1` = "CDKN2A-/-",
                  `sgRNA3` = "CDKN2A-/-",
                  `sgNTC` = "NTC",
                  `Parental` = "Parental")) 


parent <- malignant_pval_freq_wide[malignant_pval_freq_wide$SampleID == "Parental", ]

arrows_df <- malignant_pval_freq_wide %>%
  filter(SampleID != "Parental") %>%
  mutate(y = parent$Undiff_Stem,
         x = parent$AC_MES,
         yend = Undiff_Stem,
         xend = AC_MES)

pdf(paste0(fig_dir, "cdkn2a_stem_acmes_arrows_manuscript.pdf"), width = 2.25, height = 2.25, bg = "transparent", useDingbats = FALSE)
ggplot(malignant_pval_freq_wide, aes(x = AC_MES, y = Undiff_Stem)) +
  geom_point(aes(color = exp_group), size = 2) +
  scale_color_manual(values=c("Control" = "#2166ac", 
                              "CDKN2A-/-"="#b2182b")) +
  geom_segment(data = arrows_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(type = "closed", length = unit(0.2, "cm"), angle = 25),
               linetype = "dashed",
               size = 0.25,  # thin line
               color = "black") +
  geom_text(
    aes(label = exp_type),
    nudge_x = -2.25,    
    nudge_y = 0,    
    size = 2.25
  ) +
  labs(x = "Astrocyte lineage (AC/MES-like %)", y = "Stem-like (Undiff./NPC-like %)", shape = "Group", color = "Condition") +
  plot_theme +
  ylim(85,93) +
  xlim(0,15) +
  guides(color=FALSE)
dev.off()

# The Undifferentiated values are higher for this experiment compared with co-culture. Likely due to the presence of more cell state diversity in co-culture model, including irradiation. Plus, different laboratories.
# The UMAP for CDKN2A indicates more transcriptional heterogeneity across sgRNA targeting CDKN2A-/-, but it may not manifest as cellular states.
malignant_pval_freq$State <- factor(malignant_pval_freq$State, levels = c("MES-like", "AC-like", "OPC-like", "NPC-like", "Undifferentiated", "Cycling"))

pdf(paste0(fig_dir, "cdkn2a_malignant_cells_pval_score_assignment.pdf"), width = 4, height = 3, bg = "transparent")
malignant_pval_freq %>% 
  # Omitted the cycling cells so that all frequencies sum to 100%.
  filter(State!="Cycling") %>% 
  ggplot(aes(x=SampleID, fill = factor(State), y=freq*100)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "Relative cell abundance (%)", fill="Cell State", x = "") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  plot_theme +
  theme(panel.grid.major = element_line(),
        axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

# All replicates have a high percentage of cells that are cycling.
pdf(paste0(fig_dir, "cdkn2a_malignant_cells_pval_cc_assignment.pdf"), width = 8, height = 6, bg = "transparent")
cdkn2a_state_freq_cc %>% 
  ggplot(aes(x=SampleID, fill = factor(State), y=freq*100)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "Relative cell abundance (%)", fill="Is Cycling?", x = "") +
  scale_fill_manual(values=c("Cycling" = "#6BAED6")) +
  plot_theme +
  theme(panel.grid.major = element_line(),
        axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

### END ###