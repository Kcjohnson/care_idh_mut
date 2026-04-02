##############################
### Run ArchR analyses on CAREmut multiome ATAC data
### Author: Kevin Johnson
### Updated: 2026.03.30
##############################

## ArchR creates several directories automatically when creating arrow files and ArchR projects.
workdir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac"
setwd(workdir)
fig_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/figures/archr/"

## Load necessary packages.
library(dplyr)
library(ArchR) # ArchR_1.0.2
library(parallel)
library(pheatmap)
library(chromVARmotifs) # chromVARmotifs_0.2.0

#### Set-up #####
## Check available cores. ArchR recommends setting total number of cores 1/2 to 3/4 of available cores,
n_cores = detectCores() # varies
addArchRThreads(threads = n_cores/2) 
## Each R session requires that the genome is also specified and must match alignment.
addArchRGenome("hg38")

## Load the ArchR object, The human tumor data was set to a minTSS of 4 and a minimum of 1,000 fragments.
projCARE <- loadArchRProject("Save-AllSamples-2026")

# Remove so-called doublets detected by ArchR. May not perform optimally for all samples with low heterogeneity.
# Filtering 8478 cells from ArchRProject.
projCARE_filt <- filterDoublets(projCARE)

# CAREmut data that was processed by Seurat and cell types were assigned based on RNA.
mut_md <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20260320.txt", sep = "\t", header = TRUE)

# For this project, only ATAC profiles were generated for samples processed by the Verhaak lab.
mut_md_verhaak <- mut_md %>% 
  dplyr::filter(lab=="Verhaak lab")

## Create a data.frame that contains the essential information to be merged and inspected.
atac_df <- data.frame(projCARE_filt$cellNames, projCARE_filt@cellColData$Sample, projCARE_filt@cellColData$TSSEnrichment)

# The cell names are a little different between RNA (Seurat) and ATAC (ArchR). Need to create a common linker.
atac_df$CellID <- gsub("#", "-", atac_df$projCARE_filt.cellNames)

# ~118K cells present in ATACseq data that pass doublet removal and are also found in snRNAseq data that passed QC for IDH-mutant.
sum(atac_df$CellID%in%mut_md_verhaak$CellID)

# Combine the two data.frames by adding on the RNA data and filtering out what's leftover. There are about ~60,000 ATAC nuclei that are removed when focusing on the cells that pass QC for both.
atac_df_filt_rna <- atac_df %>% 
  inner_join(mut_md_verhaak, by="CellID") 

## Identify which cells to keep in the analysis. 
tmp <- getCellNames(projCARE_filt)
rna_cells = tmp[which(tmp%in%atac_df_filt_rna$projCARE_filt.cellNames)]  

## Subset the cells to only those that also have RNA.
projCARE_filt_rna <- subsetCells(ArchRProj = projCARE_filt, cellNames = rna_cells)

## Do these cell names retain the same order?
ifelse(all(getCellNames(projCARE_filt_rna)==atac_df_filt_rna$projCARE_filt.cellNames),
       sprintf("All cell names match. Proceed"), sprintf("Warning! Cell names do not match!"))


## Add a few RNA features/variables to the ArchR object.
projCARE_filt_rna$lab <- atac_df_filt_rna$lab
projCARE_filt_rna$idh_codel_subtype <- atac_df_filt_rna$idh_codel_subtype
projCARE_filt_rna$nFeature_RNA <- atac_df_filt_rna$nFeature_RNA
projCARE_filt_rna$sample_barcode <- atac_df_filt_rna$sample_barcode
projCARE_filt_rna$care_id <- atac_df_filt_rna$care_id

atac_df_filt_rna$patient_id <- sapply(strsplit(atac_df_filt_rna$care_id, "T"), "[[", 1)
atac_df_filt_rna$timepoint <- paste0("T", sapply(strsplit(atac_df_filt_rna$care_id, "T"), "[[", 2))

projCARE_filt_rna$patient_id <- atac_df_filt_rna$patient_id
projCARE_filt_rna$timepoint <- atac_df_filt_rna$timepoint
projCARE_filt_rna$CellType_final <- atac_df_filt_rna$CellType_final

paste0("Memory Size = ", round(object.size(projCARE_filt_rna) / 10^6, 3), " MB")

# Which matrices have already been created? Should be the GeneScoreMatrix and the TileMatrix.
getAvailableMatrices(projCARE_filt_rna)

# nFrags (ATAC) and nFeatures (RNA) are positively correlated: 0.49 rho. Higher quality cells across or similar sequencing depth across both modalities.
cor.test(projCARE_filt_rna$nFrags, projCARE_filt_rna$nFeature_RNA, method = "s")

## A few QC plots for the cohort.
df <- getCellColData(projCARE_filt_rna, select = c("log10(nFrags)", "TSSEnrichment"))
density_plot <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed") + theme(
  axis.text  = element_text(size = 6),
  axis.title = element_text(size = 7),
  title      = element_text(size = 7)
)

pdf(paste0(fig_dir, "all_celltypes_archr_all_cells_nfrags_vs_tss.pdf"), width = 3, height = 3, useDingbats = FALSE, bg = "transparent")
density_plot
dev.off()

p <- plotGroups(
  ArchRProj = projCARE_filt_rna, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin"
)


pdf(paste0(fig_dir, "all_celltypes_archr_all_cells_tss_scores.pdf"), width = 8, height = 5, useDingbats = FALSE, bg = "transparent")
p
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ###
##### addModuleScore from DEGs based on snRNA #######
### ### ### ### ### ### ### ### ### ### ### ### ### ###
marker_genes_top_hits <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/misc/caremut_major_celltype_deg.txt", sep = "\t", header = TRUE)

# Need to filter to those genes present in the dataset.
all_genes <- getFeatures(projCARE_filt_rna)

filtered_markers <- marker_genes_top_hits %>% 
  dplyr::select(feature, group) %>% 
  filter(feature%in%all_genes)

features_by_cell_type <- split(filtered_markers$feature, filtered_markers$group)

# The relative gene activity for each set of DEGs. These can be used to identify
set.seed(1)
projCARE_filt_rna <- addModuleScore(projCARE_filt_rna, features = features_by_cell_type, useMatrix = "GeneScoreMatrix")

df <- getCellColData(projCARE_filt_rna, select = c("CellType_final",  "Module.Malignant", "Module.Oligodendrocyte", "Module.Myeloid", "Module.Lymphocyte", "Module.ExcNeuron"))
ggplot(df %>% data.frame(), aes(x=CellType_final, y=Module.Malignant)) +
  geom_boxplot()
ggplot(df %>% data.frame(), aes(x=CellType_final, y=Module.Oligodendrocyte)) +
  geom_boxplot()
ggplot(df %>% data.frame(), aes(x=CellType_final, y=Module.Myeloid)) +
  geom_boxplot()

### ### ### ### ### ### ### ### ### ### ### ### ### ###
##### Dimensionality reduction and clustering #######
### ### ### ### ### ### ### ### ### ### ### ### ### ###
set.seed(1)

## Important to remember that inaccessible chromatin via ATAC can be "non-accessible" or "not sampled". 1s have information and 0s do not.
## These are largely default parameters:
projCARE_filt_rna <- addIterativeLSI(
  ArchRProj = projCARE_filt_rna,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( # See Seurat::FindClusters. Change parameters depending on analysis goal.
    resolution = c(0.2), 
    sampleCells = 20000, 
    n.start = 10
  ), 
  varFeatures = 30000, 
  dimsToUse = 1:30
)

## Clustering is performed with the same methods from scRNAseq relying on Seurat's functionality here.
## Default parameters are shown. Change parameters depending on selection.
projCARE_filt_rna <- addClusters(
  input = projCARE_filt_rna,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.6
)

## Inspect the cluster IDs and look at the number of cells in each cluster.
## For multiple samples, create a confusion matrix for ATAC-defined clusters and RNA-based cell type.
confusion_mat <- confusionMatrix(paste0(projCARE_filt_rna$Clusters), paste0(projCARE_filt_rna$CellType_final))


## Plot the confusion matrix using pheatmap.
confusion_mat <- confusion_mat / Matrix::rowSums(confusion_mat)
p <- pheatmap::pheatmap(
  mat = as.matrix(confusion_mat), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

## Output the confusion matrix
plotPDF(p, name = "CAREmut-Confusion-Matrix-CellType-Clusters.pdf", ArchRProj = projCARE_filt_rna, addDOC = FALSE, width = 5, height = 5)

# What's the breakdown of proposed agreement between the two approaches
projCARE_filt_rna$CellType_ATAC <- projCARE_filt_rna$Clusters
projCARE_filt_rna$CellType_ATAC <- ifelse(projCARE_filt_rna$CellType_ATAC%in%c("C18", "C11", "C12", "C8", "C15", "C5", "C10", "C9", "C6", "C17","C7", "C16", "C13","C24","C25", "C19","C23","C14","C21"), "Malignant", projCARE_filt_rna$CellType_ATAC)
projCARE_filt_rna$CellType_ATAC <- ifelse(projCARE_filt_rna$CellType_ATAC%in%c("C20"), "Astrocyte", projCARE_filt_rna$CellType_ATAC)
projCARE_filt_rna$CellType_ATAC <- ifelse(projCARE_filt_rna$CellType_ATAC%in%c("C4"), "Lymphocyte", projCARE_filt_rna$CellType_ATAC)
projCARE_filt_rna$CellType_ATAC <- ifelse(projCARE_filt_rna$CellType_ATAC%in%c("C2"), "Neuron", projCARE_filt_rna$CellType_ATAC)
projCARE_filt_rna$CellType_ATAC <- ifelse(projCARE_filt_rna$CellType_ATAC%in%c("C3"), "Endo/Mural", projCARE_filt_rna$CellType_ATAC)
projCARE_filt_rna$CellType_ATAC <- ifelse(projCARE_filt_rna$CellType_ATAC%in%c("C1"), "Myeloid", projCARE_filt_rna$CellType_ATAC)
projCARE_filt_rna$CellType_ATAC <- ifelse(projCARE_filt_rna$CellType_ATAC%in%c("C22"), "Oligodendrocyte", projCARE_filt_rna$CellType_ATAC)
table(projCARE_filt_rna$CellType_ATAC, projCARE_filt_rna$CellType_final)

confusion_mat <- confusionMatrix(paste0(projCARE_filt_rna$CellType_ATAC), paste0(projCARE_filt_rna$CellType_final))
## Plot the confusion matrix using pheatmap.
confusion_mat <- confusion_mat / Matrix::rowSums(confusion_mat)
p <- pheatmap::pheatmap(
  mat = as.matrix(confusion_mat), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

plotPDF(p, name = "CAREmut-Confusion-Matrix-CellType-RNA-ATAC.pdf", ArchRProj = projCARE_filt_rna, addDOC = FALSE, width = 5, height = 5)

## Single nuclei embeddings. Parameter selection will be analysis dependent.
set.seed(1)
projCARE_filt_rna <- addUMAP(
  ArchRProj = projCARE_filt_rna, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.3, 
  metric = "cosine"
)

# Create different plots based on various metadata features.
p1 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "idh_codel_subtype", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "timepoint", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "CellType_final", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "CellType_ATAC", embedding = "UMAP")
p7 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "nFeature_RNA", embedding = "UMAP")
p8 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")

# Output
plotPDF(p1,p2, p3, p4, p5, p6, p7, p8, name = "Plot-UMAP-All-Cell-Types.pdf", ArchRProj = projCARE_filt_rna, addDOC = FALSE, width = 5, height = 5)


p0 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "CellType_final", embedding = "UMAP")
p1 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Module.Malignant", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Module.Astrocyte", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Module.Endothelial", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Module.ExcNeuron", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Module.InhNeuron", embedding = "UMAP")
p6 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Module.Lymphocyte", embedding = "UMAP")
p7 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Module.Malignant", embedding = "UMAP")
p8 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Module.Myeloid", embedding = "UMAP")
p9 <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Module.Oligodendrocyte", embedding = "UMAP")

plotPDF(p0, p1,p2, p3, p4, p5, p6, p7, p8,p9, name = "Plot-UMAP-All-Cell-Types-Modules.pdf", ArchRProj = projCARE_filt_rna, addDOC = FALSE, width = 5, height = 5)

# Produce the final UMAPs to be included in the manuscript.
# Use the color scheme from the snRNA data
cols <- c("#BFBADA", "#FCCDE5", "#BC80BD", "#FFED6F", "#8DD3C7", "#FB8072",
          "#FFFFB3", "#80B1D3", "#B3DE69", "gray80")
names(cols) <- names(table(projCARE_filt_rna$CellType_final)) 

# Use the GLASS color codes for astrocytoma and oligodendroglioma.
projCARE_filt_rna$tumor_type <- ifelse(projCARE_filt_rna$idh_codel_subtype=="IDHmut-codel", "Oligo.", "Astro.")
idh_cols <- c("#800074", "#298C8C")
names(idh_cols) <- names(table(projCARE_filt_rna$tumor_type)) 

p_rna <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "CellType_final", embedding = "UMAP", pal = cols)
p_atac <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "tumor_type", embedding = "UMAP", pal = idh_cols)

p_myeloid_mod <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Module.Myeloid", embedding = "UMAP") +
  labs(x="UMAP1", y="UMAP2", fill="Myeloid\nATAC module score", title="Myeloid") +
  guides(fill="none")
p_oligo_mod <- plotEmbedding(ArchRProj = projCARE_filt_rna, colorBy = "cellColData", name = "Module.Oligodendrocyte", embedding = "UMAP") +
  labs(x="UMAP1", y="UMAP2", fill="ATAC module score",title="Oligodendrocyte")


pdf(paste0(fig_dir, "caremut_atac_celltype_module_score_umap.pdf"), width = 5, height = 7, useDingbats = FALSE, bg = "transparent")
ggAlignPlots(p_myeloid_mod, p_oligo_mod, type = "v")
dev.off()

ggsave(paste0(fig_dir, "caremut_atac_celltype_umap.pdf"), p_rna, width = 5, height = 5, dpi = 300)

patient_tss <- plotGroups(
  ArchRProj = projCARE_filt_rna, 
  groupBy = "patient_id", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  maxCells = 10000,
) + labs(x="")

patient_tss <- patient_tss + theme(
  axis.text  = element_text(size = 6),
  axis.title = element_text(size = 7),
  title      = element_text(size = 7)
)

pdf(paste0(fig_dir, "per_patient_tss_scores.pdf"), width = 8, height = 5, useDingbats = FALSE, bg = "transparent")
patient_tss
dev.off()


patient_nfrags <- plotGroups(
  ArchRProj = projCARE_filt_rna, 
  groupBy = "patient_id", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  maxCells = 10000
) + labs(x="")

patient_nfrags <- patient_nfrags + theme(
  axis.text  = element_text(size = 6),
  axis.title = element_text(size = 7),
  title      = element_text(size = 7)
)

pdf(paste0(fig_dir, "per_patient_nfrags.pdf"), width = 8, height = 5, useDingbats = FALSE, bg = "transparent")
patient_nfrags
dev.off()

# This plot represents 118,180 nuclei
pdf(paste0(fig_dir, "per_patient_tss_vs_nfrags.pdf"), width = 3.5, height = 5, useDingbats = FALSE, bg = "transparent")
ggAlignPlots(patient_nfrags, patient_tss, type = "v")
dev.off()


## Extract the information from these cells to be able to feed back into RNA labels.
table(projCARE_filt_rna$CellType_ATAC, projCARE_filt_rna$CellType_final)

# Most common unresolved (RNA) for ATAC were malignant and astrocyte, as expected.
atac_out <- data.frame(projCARE_filt_rna@cellColData)
write.table(atac_out, file = "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/atac/caremut_all_cell_types_cluster_assignment.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

atac_out_filt <- atac_out %>% 
  mutate(CellType_final = recode(CellType_final, `Mural` = "Endo/Mural",
                                 `Endothelial` = "Endo/Mural",
                                 `InhNeuron` = "Neuron",
                                 `ExcNeuron` = "Neuron")) %>% 
  filter(CellType_final!="Unresolved")

# IMPORTANT - 98.7% of cells would have clusters that are defined as the same cell type (i.e., Malignant==Malignant, Myeloid==Myeloid etc).
# Note: That we have not defined unresolved cells in ATAC dataset, rather it is excluded for this comparison.
sum(atac_out_filt$CellType_final==atac_out_filt$CellType_ATAC)/length(atac_out_filt$CellType_final)

# Save output where new images will be deposited and data will reside. 
saveArchRProject(ArchRProj = projCARE_filt_rna, 
                 outputDirectory = "Save-CAREmut-All-RNA", 
                 load = FALSE, 
                 dropCells = TRUE,
                 overwrite = TRUE) 





### END ###