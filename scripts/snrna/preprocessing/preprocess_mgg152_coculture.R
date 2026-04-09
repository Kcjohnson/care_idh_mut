##################################
# 10x scRNA processing of MGG152 and mouse macrophage co-culture irradiation experiments
# Author: Kevin Johnson
# Date Updated: 2026.04.05
##################################

library(tidyverse)
library(RColorBrewer)
library(viridis)
library(Seurat)
library(harmony)
library(Matrix)
library(ggpubr)
library(openxlsx)
library(DoubletFinder)

# Notes: These results were based on cellranger 9.0.1 and used the hybrid human and mouse transcriptome
# There are 4 experimental conditions, these were repeated across 4 biological replicates (that is, "-1"/RV11, "-2"/RV12, etc.).
# Batches 1 and 2 were processed for 10x and sequenced together. Batches 3 and 4 were processed and sequenced together. Batch3/4 had lower sequencing.

fig_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/figures/perturbation/"
out_data_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/perturbation/coculture/"
setwd("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/")
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/caremut_utils.R")

## Specify directories that have the count matrices. These are spread out across a few directories because we used 10x OCM library preps.
in_data_dir_1 <- "/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/cellranger/coculture/RV11/outs/per_sample_outs"
in_data_dir_2 <- "/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/cellranger/coculture/RV12/outs/per_sample_outs"
in_data_dir_3 <- "/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/cellranger/coculture/LZ1/outs/per_sample_outs"
in_data_dir_4 <- "/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/cellranger/coculture/LZ2/outs/per_sample_outs"

## List all samples in each directory.
SAMPLES_batch1 <- dir(in_data_dir_1, full.names = T)
SAMPLES_batch2 <- dir(in_data_dir_2, full.names = T)
SAMPLES_batch3 <- dir(in_data_dir_3, full.names = T)
SAMPLES_batch4 <- dir(in_data_dir_4, full.names = T)

# Create function to grab the last element of each vector (varying length vectors) from file path.
last <- function(x) { return( x[length(x)] ) }
SAMPLES_batch1_coculture <- sapply(strsplit(SAMPLES_batch1, "/"), last)
SAMPLES_batch2_coculture <- sapply(strsplit(SAMPLES_batch2, "/"), last)
SAMPLES_batch3_coculture <- sapply(strsplit(SAMPLES_batch3, "/"), last)
SAMPLES_batch4_coculture <- sapply(strsplit(SAMPLES_batch4, "/"), last)

SAMPLES <- c(SAMPLES_batch1_coculture, SAMPLES_batch2_coculture, SAMPLES_batch3_coculture, SAMPLES_batch4_coculture)

### ### ### ### ### ### ### ###
# Create seurat object:
### ### ### ### ### ### ### ###
# First 10x batch of MGG152 cells for 4 conditions
seurat_list_batch1 <- lapply(SAMPLES_batch1_coculture, function(sample_name){
  cur_data <- Read10X(paste0(in_data_dir_1, "/", sample_name,  "/count/sample_filtered_feature_bc_matrix"))
  
  ## Adding the sample identifier to the cell barcode.
  colnames(cur_data) <- paste0(sample_name, "-", colnames(cur_data))
  
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3,
    min.features=200,
    project='coculture'
  )
  cur_seurat$SampleID <- sample_name
  return(cur_seurat)
})

# Second batch submitted alongside the first batch.
seurat_list_batch2 <- lapply(SAMPLES_batch2_coculture, function(sample_name){
  cur_data <- Read10X(paste0(in_data_dir_2, "/", sample_name,  "/count/sample_filtered_feature_bc_matrix"))
  
  ## Adding the sample identifier to the cell barcode.
  colnames(cur_data) <- paste0(sample_name, "-", colnames(cur_data))
  
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3,
    min.features=200,
    project='coculture'
  )
  cur_seurat$SampleID <- sample_name
  return(cur_seurat)
})

# Third and fourth batch submitted together but had the same 10x issue
seurat_list_batch3 <- lapply(SAMPLES_batch3_coculture, function(sample_name){
  cur_data <- Read10X(paste0(in_data_dir_3, "/", sample_name,  "/count/sample_filtered_feature_bc_matrix"))
  
  ## Adding the sample identifier to the cell barcode.
  colnames(cur_data) <- paste0(sample_name, "-", colnames(cur_data))
  
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3,
    min.features=200,
    project='coculture'
  )
  cur_seurat$SampleID <- sample_name
  return(cur_seurat)
})

seurat_list_batch4 <- lapply(SAMPLES_batch4_coculture, function(sample_name){
  cur_data <- Read10X(paste0(in_data_dir_4, "/", sample_name,  "/count/sample_filtered_feature_bc_matrix"))
  
  ## Adding the sample identifier to the cell barcode.
  colnames(cur_data) <- paste0(sample_name, "-", colnames(cur_data))
  
  cur_seurat <- CreateSeuratObject(
    counts = cur_data,
    min.cells=3,
    min.features=200,
    project='coculture'
  )
  cur_seurat$SampleID <- sample_name
  return(cur_seurat)
})
# Merge seurat objects. Confirm that the gene names are as expected.
seurat_obj_batch1 <- merge(x=seurat_list_batch1[[1]], y=seurat_list_batch1[2:length(seurat_list_batch1)])
seurat_obj_batch2 <- merge(x=seurat_list_batch2[[1]], y=seurat_list_batch2[2:length(seurat_list_batch2)])
seurat_obj_batch3 <- merge(x=seurat_list_batch3[[1]], y=seurat_list_batch3[2:length(seurat_list_batch2)])
seurat_obj_batch4 <- merge(x=seurat_list_batch4[[1]], y=seurat_list_batch4[2:length(seurat_list_batch2)])

# Clean up large objects
rm(seurat_list_batch1)
rm(seurat_list_batch2)
rm(seurat_list_batch3)
rm(seurat_list_batch4)
gc()

seurat_obj_12 <- merge(x=seurat_obj_batch1, y=seurat_obj_batch2)
seurat_obj_34 <- merge(x=seurat_obj_batch3, y=seurat_obj_batch4)
seurat_obj <- merge(x=seurat_obj_12, y=seurat_obj_34)

## Save the files
saveRDS(seurat_obj, paste0(out_data_dir, "coculture_seurat_obj_unfiltered_20250510.RDS"))
# seurat_obj <- readRDS(paste0(out_data_dir, "coculture_seurat_obj_unfiltered_20250510.RDS"))

# Clean up large objects
rm(seurat_obj_batch1)
rm(seurat_obj_batch2)
rm(seurat_obj_batch3)
rm(seurat_obj_batch4)
rm(seurat_obj_12)
rm(seurat_obj_34)
gc()


#########################################################################################
##### QC assessment
#########################################################################################
# Each sample has an assignment of gem to a species. Let's inspect how it perform.
csv_files_dir1 <- list.files(path = in_data_dir_1, pattern = "gem_classification.csv", recursive = TRUE, full.names = T)
csv_files_dir2 <- list.files(path = in_data_dir_2, pattern = "gem_classification.csv", recursive = TRUE, full.names = T)
csv_files_dir3 <- list.files(path = in_data_dir_3, pattern = "gem_classification.csv", recursive = TRUE, full.names = T)
csv_files_dir4 <- list.files(path = in_data_dir_4, pattern = "gem_classification.csv", recursive = TRUE, full.names = T)

csv_files_path <- rbind(csv_files_dir1, csv_files_dir2, csv_files_dir3, csv_files_dir4)

get_sample_id <- function(path) {
  matches <- regmatches(path, regexpr("per_sample_outs/[^/]+", path))
  gsub("per_sample_outs/", "", matches)
}

# Read and combine all files with SampleID added
species_class_df <- lapply(csv_files_path, function(path) {
  sample_id <- get_sample_id(path)
  df <- read.csv(path)
  df$SampleID <- sample_id
  return(df)
}) %>% bind_rows()

# Create CellID to merge with other analyses.
species_class_df$CellID <- paste0(species_class_df$SampleID, "-", species_class_df$barcode)

# It seems that the Multiplet detection does not work well for samples without both species. For example, very few cells were assign as human in our 100% MGG152 analyses.
table(species_class_df$SampleID, species_class_df$call)

# Add the extra metadata to the Seurat object about the experimental conditions.
seurat_obj@meta.data$CellID <- rownames(seurat_obj@meta.data)
seurat_md <- seurat_obj@meta.data
add_data <- seurat_md %>% 
  left_join(species_class_df, by=c("CellID", "SampleID"))
row.names(add_data) <- row.names(seurat_obj[[]])
seurat_obj <- AddMetaData(seurat_obj, metadata = add_data)


## Perform quality control. Use both the mitochondrial genes for humans and mouse. This might be helpful in excluding multiplets.
gene_names <- rownames(seurat_obj[["RNA"]])
seurat_obj[["percent_mt_human"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "MT-")
seurat_obj[["percent_mt_mouse"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "mt-")
seurat_obj[["percent.ribo.human"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "GRCh38-RPS|GRCh38-RPL")
seurat_obj[["percent.ribo.mm"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "mm10---Rps|mm10---Rpl")

VlnPlot(seurat_obj, group.by="SampleID", features = "nFeature_RNA", ncol = 1, pt.size=0) + NoLegend() + geom_hline(yintercept = 500, linetype = "dashed", color = "red")
VlnPlot(seurat_obj, group.by="SampleID", features = "percent_mt_human", ncol = 1, pt.size=0) + NoLegend() + geom_hline(yintercept = 10, linetype = "dashed", color = "red")
VlnPlot(seurat_obj, group.by="SampleID", features = "percent_mt_mouse", ncol = 1, pt.size=0) + NoLegend() + geom_hline(yintercept = 10, linetype = "dashed", color = "red")
VlnPlot(seurat_obj, group.by="SampleID", features = "percent.ribo.human", ncol = 1, pt.size=0) + NoLegend() + geom_hline(yintercept = 20, linetype = "dashed", color = "red")
VlnPlot(seurat_obj, group.by="SampleID", features = "percent.ribo.mm", ncol = 1, pt.size=0) + NoLegend() + geom_hline(yintercept = 20, linetype = "dashed", color = "red")

pdf(paste0(fig_dir, "nFeature_rna_cutoff.pdf"), width=6, height=4)
hist(seurat_obj$nFeature_RNA) + abline(v=1000, col='red', lwd=3, lty='dashed')
dev.off()

# In the past, we set a relatively liberal criteria in snRNAseq data of 500 genes detected. Increasing it to 1,000 since we have higher quality data with scRNAseq.
# When I inspected the cells that had been removed between 500-1000 genes it did not appear enriched for many macrophages so keeping it at 1,000.
sum(seurat_obj$nFeature_RNA<500 & seurat_obj$call=="mm10") # 36
# We lose an additional 164 mouse cells - this might be okay
sum(seurat_obj$nFeature_RNA<1000 & seurat_obj$call=="mm10") # 473
sum(seurat_obj$nFeature_RNA>10000) # 803
sum(seurat_obj$percent_mt_human>10) # 4437
sum(seurat_obj$percent_mt_human>10 & seurat_obj$call=="mm10") # 6
sum(seurat_obj$percent_mt_mouse>10 & seurat_obj$call=="mm10") # 53

## Assess the changes in cell number following QC. Filter cells based on general quality assessment in scRNAseq data
pre_qc_cell_num <- dim(seurat_obj)[2]
seurat_obj <- subset(seurat_obj, nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent_mt_human < 10 & percent_mt_mouse <10 & percent.ribo.human < 30 & percent.ribo.mm < 30)
dropped_cells <- pre_qc_cell_num-dim(seurat_obj)[2]
# 10616 cells removed (if we use looser threshold - it is dropping cells across the board)
sprintf("%s cells were removed due to low quality", dropped_cells)

# There are major differences between the two batches, but samples within are similar. May want to separate.
pdf(paste0(fig_dir, "qc_violin_plot_postfilter_all_reps.pdf"), width=10, height=10)
VlnPlot(seurat_obj, group.by="SampleID", features = "nFeature_RNA", ncol = 1, pt.size=0) + NoLegend() + geom_hline(yintercept = 1001, linetype = "dashed", color = "red")
dev.off()


table(seurat_obj$call, seurat_obj$SampleID)

#########################################################################################
##### Doublet Finder
#########################################################################################
# I followed the following example: https://rpubs.com/kenneditodd/doublet_finder_example

# Split aggregated data by sample
coculture_split <- SplitObject(seurat_obj, split.by = "SampleID") 

# loop through samples to find doublets
for (i in 1:length(coculture_split)) {
  # print the sample we are on
  print(paste0("Sample ",i))
  
  # Pre-process seurat object with standard seurat workflow
  coculture_sample <- NormalizeData(coculture_split[[i]])
  coculture_sample <- FindVariableFeatures(coculture_sample)
  coculture_sample <- ScaleData(coculture_sample)
  coculture_sample <- RunPCA(coculture_sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- coculture_sample[["pca"]]@stdev
  sum.stdv <- sum(coculture_sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # finish pre-processing
  coculture_sample <- RunUMAP(coculture_sample, dims = 1:min.pc)
  coculture_sample <- FindNeighbors(object = coculture_sample, dims = 1:min.pc)              
  coculture_sample <- FindClusters(object = coculture_sample, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep_v3(coculture_sample)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bimodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- coculture_sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(coculture_sample@meta.data)) ## Assuming 7.5% doublet formation rate
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  coculture_sample <- doubletFinder_v3(seu = coculture_sample, 
                                       PCs = 1:min.pc, 
                                       pK = optimal.pk,
                                       nExp = nExp.poi.adj)
  metadata <- coculture_sample@meta.data
  # Get the number of columns in the data frame
  num_columns <- ncol(metadata)
  
  # Change the name of the last column
  colnames(metadata)[num_columns] <- "doublet_finder"
  
  coculture_sample@meta.data <- metadata 
  
  # subset and save
  caremut_singlets <- subset(coculture_sample, doublet_finder == "Singlet")
  coculture_split[[i]] <- caremut_singlets
  remove(caremut_singlets)
}

print("Finished detecting singlets")


# Converge the co-culture splits after identifying the singlets.
seurat_obj_singlets <- merge(x=coculture_split[[1]], 
                             y=coculture_split[2:length(coculture_split)],
                             project = "coculture")

columns_to_remove <- grep("^pANN", colnames(seurat_obj_singlets@meta.data), value = TRUE)
seurat_obj_singlets@meta.data <- seurat_obj_singlets@meta.data[, !colnames(seurat_obj_singlets@meta.data) %in% columns_to_remove]

## How many cells were labelled as doublets? 17,438 cells were removed as doublets by DoubletFinder.
seurat_obj
seurat_obj_singlets
table(seurat_obj$SampleID)-table(seurat_obj_singlets$SampleID)

# Save output for singlets so it is not necessary to re-run.
saveRDS(seurat_obj_singlets, paste0(out_data_dir, "coculture_seurat_obj_filt_singlets_20250510.RDS"))
print("Saved Seurat singlets")

# Load the data back in - if returning to the script
# seurat_obj_singlets <- readRDS(paste0(out_data_dir, "coculture_seurat_obj_filt_singlets_20250510.RDS"))

# We had lower sequencing depths in the 3rd and 4th batch.
pdf(paste0(fig_dir, "mgg152_coculture_nfeature_rna.pdf"), width = 10, height = 6, useDingbats = FALSE, bg = "transparent")
VlnPlot(seurat_obj_singlets, group.by="SampleID", features = "nFeature_RNA", ncol = 1, pt.size=0) + NoLegend() + geom_hline(yintercept = 1000, linetype = "dashed", color = "red")
dev.off()

# Percent ribosomal RNA was consistent across batches.
pdf(paste0(fig_dir, "mgg152_coculture_percent_ribo.pdf"), width = 10, height = 6, useDingbats = FALSE, bg = "transparent")
VlnPlot(seurat_obj_singlets, group.by="SampleID", features = "percent.ribo.human", ncol = 1, pt.size=0) + NoLegend()
dev.off()

# As expected the co-culture vs. monoculture models had more signal for mouse ribo/mito
VlnPlot(seurat_obj_singlets, group.by="SampleID", features = "percent.ribo.mm", ncol = 1, pt.size=0) + NoLegend()
VlnPlot(seurat_obj_singlets, group.by="SampleID", features = "percent_mt_mouse", ncol = 1, pt.size=0) + NoLegend()

# A little more variation in the mtDNA with a subtle low peak for all the co-culture experiments.
pdf(paste0(fig_dir, "mgg152_coculture_percent_mito_human.pdf"), width = 10, height = 6, useDingbats = FALSE, bg = "transparent")
VlnPlot(seurat_obj_singlets, group.by="SampleID", features = "percent_mt_human", ncol = 1, pt.size=0) + NoLegend()
dev.off()


# Renaming the different conditions and batches.
seurat_obj_singlets@meta.data$plot_id <- seurat_obj_singlets@meta.data$SampleID
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="152-1"] <- "Rep1_mal" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="152-2"] <-  "Rep2_mal" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="152-3"] <-  "Rep3_mal" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="152-4"] <-  "Rep4_mal" 

seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="152-r-1"] <- "Rep1_mal_RT" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="152-r-2"] <-  "Rep2_mal_RT" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="152-r-3"] <-  "Rep3_mal_RT" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="152-r-4"] <-  "Rep4_mal_RT" 

seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="m152-1"] <- "Rep1_mal+mac" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="m152-2"] <-  "Rep2_mal+mac" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="m152-3"] <-  "Rep3_mal+mac" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="m152-4"] <-  "Rep4_mal+mac" 

seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="m152-r-1"] <- "Rep1_mal+mac_RT" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="m152-r-2"] <-  "Rep2_mal+mac_RT" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="m152-r-3"] <-  "Rep3_mal+mac_RT" 
seurat_obj_singlets@meta.data$plot_id[seurat_obj_singlets@meta.data$SampleID=="m152-r-4"] <-  "Rep4_mal+mac_RT" 

# These are grouped by biological replicate.
pdf(paste0(fig_dir, "pass_qc_nfeature_rna_violin_mgg152_coculture.pdf"), width=3.5, height=3, useDingbats = FALSE)
VlnPlot(seurat_obj_singlets, group.by="plot_id", features = "nFeature_RNA", ncol = 1, pt.size=0) + 
  geom_hline(yintercept = c(1000, 10000), linetype = "dashed", color="red") +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    color = "black",
    fill = "white"
  )  +
  NoLegend() + ylim(0, 10000) + xlab("Sample ID") + ylab("nFeature RNA") + theme(plot.title = element_blank()) +
  plot_theme
dev.off()

pdf(paste0(fig_dir, "pass_qc_mito_rna_violin_mgg152_coculture.pdf"), width=3.5, height=2.5, useDingbats = FALSE)
VlnPlot(seurat_obj_singlets, group.by="plot_id", features = "percent_mt_human", ncol = 1, pt.size=0) + 
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    color = "black",
    fill = "white"
  )  +
  geom_hline(yintercept = 10, linetype = "dashed", color="red") +
  NoLegend() + ylim(0, 15) + xlab("Sample ID") + ylab("Percent mitochrondrial genes (%)") + theme(plot.title = element_blank()) +
  plot_theme
dev.off()

#########################################################################################
##### Preprocess with Standard processing  
#########################################################################################
# Create a variable for clear technical batch variable.
seurat_obj_singlets@meta.data$technical_batch <- paste0("batch", sapply(strsplit(seurat_obj_singlets@meta.data$SampleID, "-"), last))

# Run through a standard Seurat workflow
mgg152_mac_rt_sobj <- NormalizeData(seurat_obj_singlets, normalization.method = "LogNormalize", scale.factor = 10000)
mgg152_mac_rt_sobj <- FindVariableFeatures(mgg152_mac_rt_sobj, selection.method = "vst", nfeatures = 5000)
# By default only variable features are scaled
all.genes <- rownames(mgg152_mac_rt_sobj)
mgg152_mac_rt_sobj <- ScaleData(mgg152_mac_rt_sobj, features = all.genes)
mgg152_mac_rt_sobj <- RunPCA(mgg152_mac_rt_sobj, features = VariableFeatures(object = mgg152_mac_rt_sobj), npcs=50)
mgg152_mac_rt_sobj <- FindNeighbors(mgg152_mac_rt_sobj, dims = 1:20) 
mgg152_mac_rt_sobj <- FindClusters(mgg152_mac_rt_sobj, resolution = 0.6) 
mgg152_mac_rt_sobj <- RunUMAP(mgg152_mac_rt_sobj, dims = 1:20, min.dist = 0.2, n.neighbors = 40)

# First two PCs have quite the influence. Mostly human vs mouse (PC1) and then cycling/non-cycling human (PC2)
ElbowPlot(mgg152_mac_rt_sobj, reduction = "pca", ndims = 30)

# Extract the gene names to see how they are structured for this barnyard experiment.
gene_names <- rownames(seurat_obj_singlets[["RNA"]])

# Data exploration for key features.
options(ggrepel.max.overlaps = Inf)
DimPlot(mgg152_mac_rt_sobj, label = TRUE, group.by = "RNA_snn_res.0.6",repel = TRUE) 
DimPlot(mgg152_mac_rt_sobj, label = TRUE, group.by = "SampleID", repel = TRUE)
# Call variable doesn't seem very reliable in monoculture experiments, which make sense.
DimPlot(mgg152_mac_rt_sobj, label = TRUE, group.by = "call", repel = TRUE)

FeaturePlot(mgg152_mac_rt_sobj, features = "GRCh38-SOX2")
FeaturePlot(mgg152_mac_rt_sobj, features = "mm10---Ptprc")
FeaturePlot(mgg152_mac_rt_sobj, features = "GRCh38-TOP2A")
FeaturePlot(mgg152_mac_rt_sobj, features = "GRCh38-PTPRZ1")
FeaturePlot(mgg152_mac_rt_sobj, features = "GRCh38-CD44")

Idents(mgg152_mac_rt_sobj) <- "RNA_snn_res.0.6"
mgg152_mac_rt_sobj <- FindSubCluster(mgg152_mac_rt_sobj,
                                     "11",
                                     graph.name = "RNA_nn",
                                     subcluster.name = "seurat_final_clusters",
                                     resolution = 0.1,
                                     algorithm = 1
)
DimPlot(mgg152_mac_rt_sobj, label = TRUE, group.by = "seurat_final_clusters",repel = TRUE) 

# Clusters 10 and 13 are the mouse macrophages. Going to remove any potentially human cell in this cluster.
possible_doublet_ids <- rownames(mgg152_mac_rt_sobj@meta.data)[mgg152_mac_rt_sobj@meta.data$RNA_snn_res.0.6%in%c(10, 13) & mgg152_mac_rt_sobj@meta.data$call%in%c("GRCh38", "Multiplet")]
DimPlot(mgg152_mac_rt_sobj, label=T, cells.highlight = list(possible_doublet_ids), cols.highlight = c("darkblue"), cols= "grey")

# Remove the cells by subsetting all cells NOT in the doublet list (458 cells)
mgg152_mac_rt_sobj <- subset(mgg152_mac_rt_sobj, 
                             cells = setdiff(colnames(mgg152_mac_rt_sobj), possible_doublet_ids))

# Create additional identifiers
mgg152_mac_rt_sobj@meta.data$condition <- ifelse(grepl("-r-", mgg152_mac_rt_sobj@meta.data$SampleID,), "irradiation", "control") 
mgg152_mac_rt_sobj@meta.data$culture_conditions <- ifelse(grepl("m152", mgg152_mac_rt_sobj@meta.data$SampleID,), "malignant+macrophage", "malignant_only") 

DimPlot(mgg152_mac_rt_sobj, label = TRUE, group.by = "seurat_final_clusters", repel = TRUE) 
mgg152_mac_rt_sobj@meta.data$species <- ifelse(mgg152_mac_rt_sobj@meta.data$seurat_final_clusters%in%c(10, 13) & mgg152_mac_rt_sobj@meta.data$call%in%c("mm10"), "mm10", "GRCh38")

DimPlot(mgg152_mac_rt_sobj, label = TRUE, group.by = "species",repel = TRUE) 

png(paste0(fig_dir, "mgg152_coculture_umap_batch_culture.png"), width = 16, height = 8, units = 'in', res = 300, bg = "transparent")
DimPlot(mgg152_mac_rt_sobj, label = TRUE, group.by = c("technical_batch", "species"), repel = TRUE) 
dev.off()

png(paste0(fig_dir, "mgg152_coculture_umap_ptprz1_ptprc.png"), width = 8, height = 4, units = 'in', res = 300, bg = "transparent")
FeaturePlot(mgg152_mac_rt_sobj, features = c("GRCh38-PTPRZ1", "mm10---Ptprc"))
dev.off()


seurat_out <- mgg152_mac_rt_sobj@meta.data
# Plot the UMI counts per species for these defined cell types.
pdf(paste0(fig_dir, "mgg152_coculture_species_umi_plot.pdf"), width = 5, height = 4, bg = "transparent")
ggplot(seurat_out, aes(x=GRCh38, y=mm10)) +
  geom_point(aes(color=species)) +
  plot_theme +
  labs(x= "GRCh38 - UMI counts per cell", y= "mm10 - UMI counts per cell")
dev.off()


### ### ### ### ### ### ### ###
# DEG 
### ### ### ### ### ### ### ###
# Create a SampleID and species identifier to be able to perform DE tests across each compartment.
mgg152_mac_rt_sobj@meta.data$SampleID_species <- paste0(mgg152_mac_rt_sobj@meta.data$SampleID, "_", "GRCh38")
mgg152_mac_rt_sobj@meta.data$SampleID_species[mgg152_mac_rt_sobj@meta.data$call=="mm10"] <- paste0(mgg152_mac_rt_sobj@meta.data$SampleID[mgg152_mac_rt_sobj@meta.data$call=="mm10"], "_", "mm10")

Idents(mgg152_mac_rt_sobj) <- "SampleID_species"

# Macrophage-specific differential expression highlights clear, consistent irradiation response ("Cdk1na") in mouse cells but not human. 
rep1_mac_irradiation_de_results <- FindMarkers(
  object = mgg152_mac_rt_sobj,
  ident.1 = "m152-r-1_mm10",
  ident.2 = "m152-1_mm10")
rep1_mac_irradiation_de_results$gene_name <- rownames(rep1_mac_irradiation_de_results)
rep1_mac_irradiation_de_results$batch <- "batch1"

rep2_mac_irradiation_de_results <- FindMarkers(
  object = mgg152_mac_rt_sobj,
  ident.1 = "m152-r-2_mm10",
  ident.2 = "m152-2_mm10")
rep2_mac_irradiation_de_results$gene_name <- rownames(rep2_mac_irradiation_de_results)
rep2_mac_irradiation_de_results$batch <- "batch2"

rep3_mac_irradiation_de_results <- FindMarkers(
  object = mgg152_mac_rt_sobj,
  ident.1 = "m152-r-3_mm10",
  ident.2 = "m152-3_mm10")
rep3_mac_irradiation_de_results$gene_name <- rownames(rep3_mac_irradiation_de_results)
rep3_mac_irradiation_de_results$batch <- "batch3"

rep4_mac_irradiation_de_results <- FindMarkers(
  object = mgg152_mac_rt_sobj,
  ident.1 = "m152-r-4_mm10",
  ident.2 = "m152-4_mm10")
rep4_mac_irradiation_de_results$gene_name <- rownames(rep4_mac_irradiation_de_results)
rep4_mac_irradiation_de_results$batch <- "batch4"

# Cdkn1a is a strong marker of radiation response across all replicates.
macrophage_upregulated <- rep1_mac_irradiation_de_results %>% 
  bind_rows(rep2_mac_irradiation_de_results, rep3_mac_irradiation_de_results, rep4_mac_irradiation_de_results) %>% 
  filter(avg_log2FC > 0) %>% 
  group_by(gene_name) %>% 
  summarise(counts = n()) %>% 
  filter(counts > 3)


### ### ### ### ### ### ### ###
### Signature scoring
### ### ### ### ### ### ### ###
# Read in the signatures for malignant metaprograms
malignant_mp <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/metaprograms/care_mut_selected_malignant_metaprograms.csv", sep = ",", header = T, row.names = 1)
mut_mp_list <- lapply(names(malignant_mp), function(col_name) malignant_mp[[col_name]])
names(mut_mp_list) <-paste0(colnames(malignant_mp), "_MUT")

# Subset to HUMAN cells that passed QC and were not considered doublets.
mgg152_mac_rt_sobj_human <- subset(mgg152_mac_rt_sobj, call!="mm10")

# Make gene names compliant for hybrid transcriptome
malignant_mp_list_compliant <- lapply(mut_mp_list, function(gene_vec) {
  paste0("GRCh38-", gene_vec)
})


# Create new variable identifiers based on sample name.
mgg152_mac_rt_sobj_human@meta.data$condition <- ifelse(grepl("-r-", mgg152_mac_rt_sobj_human@meta.data$SampleID,), "irradiation", "control") 
mgg152_mac_rt_sobj_human@meta.data$species <- ifelse(grepl("m152", mgg152_mac_rt_sobj_human@meta.data$SampleID,), "malignant+macrophage", "malignant_only") 
mgg152_mac_rt_sobj_human@meta.data$exp_group <- paste0(mgg152_mac_rt_sobj_human@meta.data$species, "_", mgg152_mac_rt_sobj_human@meta.data$condition)



#########################################################################################
##### Assign within-sample metaprogram score - 
#########################################################################################
# Define a metadata table for which we input into the score within samples function.
mgg152_md <- mgg152_mac_rt_sobj_human@meta.data

## Create a list for analysis and name the elements
umi_data_list_batch1 <- setNames(lapply(SAMPLES_batch1_coculture, function(sample_name){
  print(sample_name)
  cur_data <- Read10X(paste0(in_data_dir_1, "/", sample_name,'/count/sample_filtered_feature_bc_matrix/'), gene.column=2)
  colnames(cur_data) <- paste0(sample_name, "-", colnames(cur_data))
  return(cur_data)
}), SAMPLES_batch1_coculture)

umi_data_list_batch2 <- setNames(lapply(SAMPLES_batch2_coculture, function(sample_name){
  print(sample_name)
  cur_data <- Read10X(paste0(in_data_dir_2, "/", sample_name,'/count/sample_filtered_feature_bc_matrix/'), gene.column=2)
  colnames(cur_data) <- paste0(sample_name, "-", colnames(cur_data))
  return(cur_data)
}), SAMPLES_batch2_coculture)

umi_data_list_batch3 <- setNames(lapply(SAMPLES_batch3_coculture, function(sample_name){
  print(sample_name)
  cur_data <- Read10X(paste0(in_data_dir_3, "/", sample_name,'/count/sample_filtered_feature_bc_matrix/'), gene.column=2)
  colnames(cur_data) <- paste0(sample_name, "-", colnames(cur_data))
  return(cur_data)
}), SAMPLES_batch3_coculture)

umi_data_list_batch4 <- setNames(lapply(SAMPLES_batch4_coculture, function(sample_name){
  print(sample_name)
  cur_data <- Read10X(paste0(in_data_dir_4, "/", sample_name,'/count/sample_filtered_feature_bc_matrix/'), gene.column=2)
  colnames(cur_data) <- paste0(sample_name, "-", colnames(cur_data))
  return(cur_data)
}), SAMPLES_batch4_coculture)


# Combine all data
umi_data_all <- c(umi_data_list_batch1, umi_data_list_batch2, umi_data_list_batch3, umi_data_list_batch4)
names(umi_data_all)

sig_lists <- c(malignant_mp_list_compliant)


## Perform the scoring within a single sample. Note that AddModuleScore results differ depending on cell/sample dataset.
# We are only scoring those cells that are included post-processing
set.seed(123)
mp_scores <- score_within_samples_caremut(umi_data_all, md = mgg152_md, sigs = sig_lists)

# saveRDS(mp_scores, paste0(out_data_dir, "coculture_mp_sig_scores.RDS"))
# mp_scores <- readRDS(paste0(out_data_dir, "coculture_mp_sig_scores.RDS"))


#####################################################################
### Cell state scores
#####################################################################
# mdata is an object (tibble) that contains the meta-data for the cells that should be classified (Sample, Treatment, etc.).
# The object also must include the CellID variable that identifies the cells.
# MP_scores (tibble) includes the scores for the meta-program and the CellID variable. All meta-programs score columns should
# start with the literal "MP_".
# Due to the shuffling each variable in the tibble is theoretically normally distributed (or is at least close to being ND).
mdata <- mgg152_md %>%
  left_join(mp_scores,
            by = c("CellID", "SampleID")) %>% 
  as_tibble() 

## The approach requires variables to be named in a certain way:
sigs <- sig_lists

# We call this function to generate a NULL distribution to facilitate classification. According to the configured parameters
# it will sample 5000 cells from the pool of cells, shuffle the expression values while maintaining the mean expression of
# each gene and score the artificial cells for the meta-programs. It will repeat the process 20 times to generate a NULL 
# distribution of 100K cells. It returns a tibble of 100K x n (where n is the number of meta-programs)
set.seed(42)
permuted_data <- generate_null_dist(umi_data_list = umi_data_all,
                                    md = mdata,
                                    sigs = sigs,
                                    n_iter = 20, n_cells = 5000, verbose = T)

# Save output of the permuted data
saveRDS(permuted_data, paste0(out_data_dir, "mgg152_coculture_permuted_scores.RDS"))
# permuted_data <- readRDS(paste0(out_data_dir, "mgg152_coculture_permuted_scores.RDS"))

mp_scores <- mp_scores
permuted_data <- permuted_data

# The signatures are the cell state programs. Do not include hierarchy scores.
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

png(paste0(fig_dir, "mgg152_coculture_mp_classification_permuted_actual.png"), width = 12, height = 8, units = 'in', res = 300, bg = "transparent")
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

# We consider the p-value to be significant if the nominal p-value is less than 0.05. 
state_data$p.sig <- state_data$p.val < .05

# Compute the classification statistics for each MP/gene set 
state_stats <- state_data %>%
  group_by(Program) %>%
  summarise(n = sum(p.sig == T), N = n(), Freq = n / N)

# This is the actual classification step. We filter out the statistically insignificant scores and classify the cell
# to the MP with the maximal signal. Classify CC separate from the other states.
# The classification will do the following:
# 1. Not consider cell cycle MP (at first)
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
# then you can adjust the multiple-adjustment method or classification threshold to less stringent values
state_vec <- setNames(rep("Undifferentiated", nrow(mdata)), mdata$CellID)
state_vec[state_data_classify$CellID] <- state_data_classify$Program
table(state_vec)
table(state_vec) / length(state_vec)

mdata$State <- state_vec[mdata$CellID]
# There are considerably more MES-like cells than I had anticipated. This could be due to co-culture conditions.
table(mdata$State)

# Bimodal distribution due to the differences in technical batch
png(paste0(fig_dir, "mgg152_caremut_mp_classification_complexity.png"), width = 10, height = 6, units = 'in', res = 300, bg = "transparent")
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
  theme(panel.grid.major = element_line()) +
  facet_grid(.~SampleID)
dev.off()

# Cell cycle is not considered a cellular state but rather a feature (since cells can have a clear identity such as OPC or NPC and still be cycling) 
mdata_out <- mdata %>% 
  mutate(isCC = ifelse(mdata$CellID%in%state_data_classify_cc$CellID, TRUE, FALSE))
table(mdata_out$State, mdata_out$isCC)

# Write out and read back in for ease of figure regeneration
write.table(mdata_out, file = paste0(out_data_dir, "mgg152_coculture_caremut_mp_select_state_assignment.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# mdata_out <- read.delim(file = paste0(out_data_dir, "mgg152_coculture_caremut_mp_select_state_assignment.txt"), sep = "\t", header = TRUE)

linker_info <- mdata_out %>% 
  dplyr::select(SampleID, species, exp_group, condition, technical_batch) %>% 
  distinct()
rownames(linker_info) <- NULL

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

malignant_pval_freq$State <- factor(malignant_pval_freq$State, levels = rev(c("MES-like", "AC-like", "OPC-like", "NPC-like", "Undifferentiated")))

pdf(paste0(fig_dir, "mgg152_malignant_cells_pval_score_assignment.pdf"), width = 10, height = 8, bg = "transparent")
malignant_pval_freq %>% 
  ggplot(aes(x=condition, fill = factor(State), y=freq*100)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(y = "Relative cell abundance (%)", fill="Cell State", x = "") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  theme_bw() +
  theme(panel.grid.major = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(technical_batch~species, scales = "free_x", space = "free_x") +
  plot_theme
dev.off()


pdf(paste0(fig_dir, "mgg152_malignant_cells_mes_pval_score_assignment_boxplot_paired_ttest.pdf"), width = 5, height = 4, bg = "transparent")
malignant_pval_freq %>% 
  filter(State=="MES-like") %>% 
  ggplot(aes(x=condition, y=freq*100)) +
  geom_boxplot(aes(fill = factor(State))) +
  geom_line(aes(group=technical_batch), color="gray70", linetype=2) +
  geom_point() +
  # T-test for small cohort size.
  stat_compare_means(method="t.test", label="p.format", paired = TRUE) +
  facet_grid(.~species, scales="free") +
  labs(y = "Relative MES-like cell abundance (%)", fill="Cell State", x = "") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(fill=FALSE) 
dev.off()

pdf(paste0(fig_dir, "mgg152_malignant_cells_mes_pval_score_assignment_boxplot_paired_wilcox.pdf"), width = 5, height = 4, bg = "transparent")
malignant_pval_freq %>% 
  filter(State=="MES-like") %>% 
  ggplot(aes(x=condition, y=freq*100)) +
  geom_boxplot(aes(fill = factor(State))) +
  geom_line(aes(group=technical_batch), color="gray70", linetype=2) +
  geom_point() +
  stat_compare_means(method="wilcox", label="p.format", paired = TRUE) +
  facet_grid(.~species, scales="free") +
  labs(y = "Relative MES-like cell abundance (%)", fill="Cell State", x = "") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(fill=FALSE) 
dev.off()

malignant_pval_freq <- malignant_pval_freq %>% 
  mutate(plot_id = case_when(exp_group == "malignant_only_control" ~ "Malignant\nmonoculture",
                             exp_group == "malignant+macrophage_control" ~ "Mal.+Mac.\nco-culture",
                             exp_group == "malignant_only_irradiation" ~ "Mal. monoculture\nirradiation",
                             exp_group == "malignant+macrophage_irradiation" ~ "Mal.+Mac. co-culture\nirradiation",
                             TRUE ~ NA_character_
                             ))
malignant_pval_freq$plot_id <- factor(malignant_pval_freq$plot_id, levels=c("Malignant\nmonoculture", "Mal.+Mac.\nco-culture", "Mal. monoculture\nirradiation", "Mal.+Mac. co-culture\nirradiation"))

mono_vs_co_plot <- malignant_pval_freq %>% 
  filter(State=="MES-like", exp_group%in%c("malignant_only_control", "malignant+macrophage_control")) %>% 
  ggplot(aes(x=plot_id, y=freq*100)) +
  geom_boxplot(aes(fill = factor(State))) +
  geom_point() +
  stat_compare_means(method="wilcox", label="p.format") +
  labs(y = "MES-like cell abundance (%) - no irradiation control", x = "Condition", fill="Cell State") +
  scale_fill_manual(values=c("AC-like" = "#AA2756", 
                             "MES-like"="#F77D58",
                             "NPC-like" = "#7fbf7b",
                             "OPC-like"="#E8F5A3",
                             "Undifferentiated" = "gray90")) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(fill=FALSE) 

pdf(paste0(fig_dir, "mgg152_mes_mono_vs_coculture_control_boxplot.pdf"), width = 4, height = 5, bg = "transparent")
mono_vs_co_plot
dev.off()

# Define a batch controlled pairwise difference in MES-like abundance.
mes_abundance_diff <-  malignant_pval_freq %>% 
  filter(State=="MES-like")  %>%
  select(State, freq, species, condition, technical_batch) %>%
  pivot_wider(names_from = condition, values_from = freq) %>%
  mutate(mes_irradiated_change = irradiation - control,
         State = "MES-like") %>% 
  mutate(plot_id = case_when(species == "malignant_only" ~ "Malignant\nmonoculture",
                             species == "malignant+macrophage" ~ "Mal.+Mac.\nco-culture",
                             TRUE ~ NA_character_
  ))
mes_abundance_diff$plot_id <- factor(mes_abundance_diff$plot_id, levels=c("Malignant\nmonoculture", "Mal.+Mac.\nco-culture"))

# We reach statistical significance if we examine the pairwise difference.
irradiated_control_plot <- ggplot(mes_abundance_diff, aes(x=plot_id, y=mes_irradiated_change*100)) +
  geom_boxplot(aes(fill = factor(State))) +
  scale_fill_manual(values=c("MES-like"="#F77D58")) +
  geom_point() +
  stat_compare_means(method = "wilcox", label="p.format") + 
  plot_theme +
  labs(x = "Condition", y = "MES-like change (%) after irradiation") +
  guides(fill = FALSE)  +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

pdf(paste0(fig_dir, "mgg152_mes_mono_vs_coculture_irradiated_vs_control_boxplot.pdf"), width = 4, height = 5, bg = "transparent")
irradiated_control_plot
dev.off()

library(cowplot)

pdf(paste0(fig_dir, "fig5e_monoculture_coculture_irradiation.pdf"), width = 4, height = 5, bg = "transparent")
plot_grid(mono_vs_co_plot, irradiated_control_plot, ncol = 2)
dev.off()


# Create an average value figure, which should help correct for any outliers or misclassified cells.
malignant_score_df_summary <- mdata_out %>% 
  group_by(SampleID, condition, species, technical_batch) %>% 
  summarise(avg_mes_score = mean(MP_MES_MUT),
            avg_ac_score = mean(MP_AC1_MUT),
            avg_ac2_score = mean(MP_AC2_MUT),
            avg_opc_score = mean(MP_OPC_MUT),
            avg_npc_score = mean(MP_NPC_MUT)) %>% 
  ungroup() %>% 
  mutate(avg_ac = (avg_ac_score+avg_ac2_score)/2)

pdf(paste0(fig_dir, "human_malignant_cells_mes_signature_score.pdf"), width = 7, height = 5, useDingbats = FALSE, bg = "transparent")
ggplot(malignant_score_df_summary, aes(x=condition, y=avg_mes_score)) +
  geom_point(aes(color=species)) +
  geom_line(aes(group=technical_batch), color="gray70", linetype=2) +
  plot_theme +
  scale_color_manual(values=c("malignant_only" = "#6BAED6", 
                              "malignant+macrophage" = "#AA2756")) +
  facet_grid(.~species, scales = "free") +
  guides(color=FALSE) +
  stat_compare_means(method="t.test", paired=TRUE) +
  labs(x = "Treatment condition", y="Mean MES-like metaprogram score", color="Culture conditions", title = "**ONLY HUMAN MALIGNANT CELLS SCORED**") 
dev.off()

### END ###