##############################
### Run ArchR differential gene accessibility on multiome ATAC data for Myeloid cells
### Author: Kevin Johnson
##############################

## Myeloid cell states gene differential gene accessibility

# Specify directories
workdir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/"
setwd(workdir)
fig_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/figures/archr/"

## Load necessary packages.
library(dplyr)
library(ArchR)
library(parallel)
library(pheatmap)
library(chromVARmotifs)
library(BSgenome.Hsapiens.UCSC.hg38)

#### Set-up #####
## Check available cores. ArchR recommends setting total number of cores 1/2 to 3/4 of available cores,
num_cores <- detectCores() # e.g., 36
n_threads <- num_cores/2
addArchRThreads(threads = n_threads) 

## Each R session requires that the genome is also specified and must match alignment.
addArchRGenome("hg38")

# Load the ArchR object for IDHmut final analysis set. Doublets have been removed, only cells that intersect with passed qc for RNA, and RNA annotated cell states.
CARE_filt_rna_all <- loadArchRProject("Save-CAREmut-All-RNA-Filtered")

# Get cell names to be able to merge with myeloid cells that were defined in higher resolution
atac_df <- data.frame(CARE_filt_rna_all$cellNames, CARE_filt_rna_all@cellColData$Sample, CARE_filt_rna_all@cellColData$TSSEnrichment)
atac_df$CellID <- gsub("#", "-", atac_df$CARE_filt_rna_all.cellNames)

# Assigned metaprogram/state score based on p-value calculation
myeloid_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/data/snrna/myeloid_cell_classification.tsv", header = TRUE, sep = "\t")

atac_df_filt <- atac_df %>% 
  inner_join(myeloid_md, by=c("CellID")) 

tmp <- getCellNames(CARE_filt_rna_all)
rna_cells = tmp[which(tmp%in%atac_df_filt$CARE_filt_rna_all.cellNames)]  

# Restrict to only myeloid cells - 11,272
CARE_filt_rna_myeloid <- subsetArchRProject(ArchRProj = CARE_filt_rna_all, cells = rna_cells, outputDirectory = "Save-CAREmut-All-Myeloid", force = TRUE)

# Set up to add this information:
CARE_filt_rna_myeloid$CellID <-  gsub("#", "-", CARE_filt_rna_myeloid$cellNames)

atac_df_filt_sub <- atac_df_filt[atac_df_filt$CellID %in% CARE_filt_rna_myeloid$CellID, ]
atac_df_filt_ord <- atac_df_filt_sub[match(CARE_filt_rna_myeloid$CellID, atac_df_filt_sub$CellID), ]
all(CARE_filt_rna_myeloid$CellID==atac_df_filt_ord$CellID)
# Do these cell names retain the same order?
ifelse(all(getCellNames(CARE_filt_rna_myeloid)==atac_df_filt_ord$atacCellNames),
       sprintf("All cell names match. Proceed"), sprintf("Warning! Cell names do not match!"))

# Add myeloid annotation
CARE_filt_rna_myeloid$CellStateCollapsed <- atac_df_filt_ord$myeloid_state_collapsed
CARE_filt_rna_myeloid$CellState <- atac_df_filt_ord$myeloid_state

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Define markers based on myeloid cell states
markersGS <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_myeloid, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "CellState",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 500
)

## Extract the marker list. I relaxed the thresholds because there were too few features
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")

# Examine the distribution of differentially accessible genes. Relatively few compare with malignant. No marker for Unresolved
lapply(markerList, nrow)
data.frame(markerList$Macrophage)
data.frame(markerList$Inflammatory)
data.frame(markerList$`Microglia`)
data.frame(markerList$`MHC-II-high`)

# Selected marker genes of interest.
markerGenes <- c("CCL3", "CD163L1", "CX3CR1")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 0.5", 
  limits = c(-1.5, 1.5),
  labelMarkers = markerGenes,
  transpose = TRUE,
)


plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap-Myeloid",  width = 3.5, height = 3,  ArchRProj = CARE_filt_rna_myeloid, addDOC = FALSE)

pdf(paste0(fig_dir, "myeloid_genescores_marker_heatmap.pdf"), width = 5, height = 4, , useDingbats = FALSE, bg = "transparent")
heatmapGS
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ###
##### Dimensionality reduction and clustering #######
### ### ### ### ### ### ### ### ### ### ### ### ### ###

## Important to remember that inaccessible chromatin via ATAC can be "non-accessible" or "not sampled". 1s have information and 0s do not.
## Needed to reduce the parameters because fewer number of myeloid cells.
CARE_filt_rna_myeloid <- addIterativeLSI(
  ArchRProj = CARE_filt_rna_myeloid,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI",
  force = TRUE,
  iterations = 2, 
  clusterParams = list( # See Seurat::FindClusters. Change parameters depending on analysis goal.
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 10000, 
  dimsToUse = 1:20
)

## Clustering is performed with the same methods from scRNAseq relying on Seurat's functionality here.
## Default parameters are shown. Change parameters depending on selection.
CARE_filt_rna_myeloid <- addClusters(
  input = CARE_filt_rna_myeloid,
  reducedDims = "IterativeLSI",
  force = TRUE,
  method = "Seurat",
  name = "Clusters",
  resolution = 0.6
)

CARE_filt_rna_myeloid <- addUMAP(
  ArchRProj = CARE_filt_rna_myeloid, 
  reducedDims = "IterativeLSI", 
  force = TRUE,
  name = "UMAP", 
  nNeighbors = 20, 
  minDist = 0.3, 
  metric = "cosine"
)

## Create different plots based on various metadata features.
p1 <- plotEmbedding(ArchRProj = CARE_filt_rna_myeloid, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = CARE_filt_rna_myeloid, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = CARE_filt_rna_myeloid, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = CARE_filt_rna_myeloid, colorBy = "cellColData", name = "CellState", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = CARE_filt_rna_myeloid, colorBy = "cellColData", name = "CellStateCollapsed", embedding = "UMAP")

# Not much in terms of separation across these ~11K cells
plotPDF(p1,p2, p3, p4, p5, name = "Plot-UMAP-Myeloid.pdf", ArchRProj = CARE_filt_rna_myeloid, addDOC = FALSE, width = 5, height = 5)

CARE_filt_rna_myeloid <- saveArchRProject(ArchRProj = CARE_filt_rna_myeloid, outputDirectory = "Save-CAREmut-All-Myeloid", load = TRUE, dropCells = TRUE) 

### END ###