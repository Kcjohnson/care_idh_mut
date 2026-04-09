##############################
### Run ArchR peak calling analysis on CAREmut multiome ATAC data for Myeloid cells
### Author: Kevin Johnson
### Updated: 2024.04.28
##############################

## Part 3: Peak calling on RNA-defined MYELOID states across IDH-mutant tumors

## ArchR creates several directories automatically when creating arrow files and ArchR projects.
workdir <- "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/"
setwd(workdir)

## Load necessary packages.
library(dplyr)
library(ArchR)
library(parallel)
library(pheatmap)
library(chromVARmotifs)
library(BSgenome.Hsapiens.UCSC.hg38)

## Specify directories:
fig_dir <- "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/figures/archr/"

#### Set-up #####
## Check available cores. ArchR recommends setting total number of cores 1/2 to 3/4 of available cores,
num_cores <- detectCores() # e.g., 36
n_threads <- num_cores/2
addArchRThreads(threads = n_threads) 
## Each R session requires that the genome is also specified and must match alignment.
addArchRGenome("hg38")

# Load the ArchR object for IDHmut final analysis set. Doublets have been removed, only cells that intersect with passed qc for RNA, and RNA annotated cell states.
CARE_filt_rna_all <- loadArchRProject("Save-CAREmut-All-RNA")

# Remove "Unresolved" (RNA-based) cells from further analysis since it may cause issues. Unresolved are mostly malignant/astrocytes.
atac_df <- data.frame(CARE_filt_rna_all$cellNames, CARE_filt_rna_all@cellColData$Sample, CARE_filt_rna_all@cellColData$TSSEnrichment)
atac_df$CellID <- gsub("#", "-", atac_df$CARE_filt_rna_all.cellNames)

# Assigned metaprogram/state score based on p-value calculation
myeloid_md <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/processed_data/rna/tme/caremut_myeloid_care_pval_classification_20240303.txt", header = TRUE)

atac_df_filt <- atac_df %>% 
  inner_join(myeloid_md, by=c("CellID")) %>% 
  filter(cell_state!="Unresolved")

tmp <- getCellNames(CARE_filt_rna_all)
rna_cells = tmp[which(tmp%in%atac_df_filt$CARE_filt_rna_all.cellNames)]  

# Subset the cells to only those that also have RNA.
CARE_filt_rna_myeloid <- subsetCells(ArchRProj = CARE_filt_rna_all, cellNames = rna_cells)


# Set up to add this information:
CARE_filt_rna_myeloid$CellID <-  gsub("#", "-", CARE_filt_rna_myeloid$cellNames)


atac_df_filt_sub <- atac_df_filt[atac_df_filt$CellID %in% CARE_filt_rna_myeloid$CellID, ]
atac_df_filt_ord <- atac_df_filt_sub[match(CARE_filt_rna_myeloid$CellID, atac_df_filt_sub$CellID), ]
all(CARE_filt_rna_myeloid$CellID==atac_df_filt_ord$CellID)
CARE_filt_rna_myeloid$CellState <- atac_df_filt_ord$cell_state


# Add grade information
patient_md <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240227.txt", sep="\t", header = TRUE)
sample_md <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240227.txt", sep="\t", header = TRUE)

atac_md_mut <- data.frame(getCellColData(CARE_filt_rna_myeloid))
atac_md_mut_annot <- atac_md_mut %>% 
  left_join(sample_md, by=c("care_id", "sample_barcode", "patient_id", "timepoint", "idh_codel_subtype")) %>% 
  mutate(atacCellNames = rownames(atac_md_mut))

## Do these cell names retain the same order?
ifelse(all(getCellNames(CARE_filt_rna_myeloid)==atac_md_mut_annot$atacCellNames),
       sprintf("All cell names match. Proceed"), sprintf("Warning! Cell names do not match!"))

## Add a few RNA features to the ArchR object.
CARE_filt_rna_myeloid$Grade <- paste0("G", atac_md_mut_annot$grade_num)
CARE_filt_rna_myeloid$subtype_grade <- paste0(CARE_filt_rna_myeloid$idh_codel_subtype, "_", CARE_filt_rna_myeloid$Grade)

# What's the breakdown of sample information:7K (IDH-O) vs 3K (IDH-A)
table(CARE_filt_rna_myeloid$idh_codel_subtype)
table(CARE_filt_rna_myeloid$subtype_grade, CARE_filt_rna_myeloid$timepoint)
table(CARE_filt_rna_myeloid$CellType_final)
table(CARE_filt_rna_myeloid$CellState)
############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Assess differential gene activity scores across the different RNA-based cell types.
#  devtools::install_github('immunogenomics/presto', repos = BiocManager::repositories())

## Define markers based on clusters. There's an average of about 4 clusters per tumor.
markersGS <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_myeloid, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "CellState",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 500
)

## Extract the marker list. These are the default thresholds for ArchR.
# Throughout this analysis, I will set a cutoff criteria of FDR < 0.05 and Log2FC >= 1. The ArchR tutorial uses variable cut-offs and I could not find a good explanation for why.
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")

# Examine the distribution of differentially accessible genes
lapply(markerList, nrow)
data.frame(markerList$`TAM_BMDM-like`)
data.frame(markerList$TAM_Inflammatory)
data.frame(markerList$`TAM_MG-like`)
data.frame(markerList$`TAM_Phagocytic`)

# Save the project so that the annotations can be more easily accessed in the future. 
CARE_filt_rna_myeloid <- saveArchRProject(ArchRProj = CARE_filt_rna_myeloid, outputDirectory = "Save-CAREmut-Myeloid-RNA", load = TRUE, dropCells = TRUE) 

CARE_filt_rna_myeloid <- loadArchRProject("Save-CAREmut-Myeloid-RNA")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", 
  transpose = TRUE)

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_myeloid, addDOC = FALSE)


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
p3 <- plotEmbedding(ArchRProj = CARE_filt_rna_myeloid, colorBy = "cellColData", name = "idh_codel_subtype", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = CARE_filt_rna_myeloid, colorBy = "cellColData", name = "CellState", embedding = "UMAP")
p5 <- plotEmbedding(ArchRProj = CARE_filt_rna_myeloid, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")

## Output
plotPDF(p1,p2, p3, p4, p5, name = "Plot-UMAP-Myeloid.pdf", ArchRProj = CARE_filt_rna_myeloid, addDOC = FALSE, width = 5, height = 5)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
## Make pseudo bulk measurements.
## There's a known issue with this on HPC systems: https://github.com/GreenleafLab/ArchR/issues/248
## Might be solved by setting threads = 1
# Pseudo-bulk refers to a grouping of single cells where the data from each single cell is combined into a single pseudo-sample
CARE_filt_rna_myeloid <- addGroupCoverages(ArchRProj = CARE_filt_rna_myeloid, 
                                            groupBy = "CellState", 
                                            minCells = 100,
                                            maxCells = 500,
                                            threads = getArchRThreads(),
                                            # Overwite the data in the ArchRProject object if the pseudo-bulk replicate information already exists
                                            force = TRUE)

## Is macs2 in the path variable?
pathToMacs2 <- findMacs2()

# Iterative overlap peak merging procedure
CARE_filt_rna_myeloid <- addReproduciblePeakSet(
  ArchRProj = CARE_filt_rna_myeloid, 
  groupBy = "CellState", 
  pathToMacs2 = pathToMacs2,
  threads = getArchRThreads(),
)

## Identifying marker peaks with ArchR.
CARE_filt_rna_myeloid <- addPeakMatrix(CARE_filt_rna_myeloid)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

## Identifying marker PEAKS - features that are unique to a specific cell grouping.
## Tell ArchR to account for biases in data quality via TSSEnrichment and nFrags.
markersPeaks <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_myeloid, 
  useMatrix = "PeakMatrix", 
  groupBy = "CellState",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 500
)


## Extract the marker peaks. Get access the GRanges object via `returnGR = TRUE`.
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")
lapply(markerList, nrow)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = TRUE,
  labelMarkers = NULL
)
heatmapPeaks
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-MyeloidStates", width = 8, height = 6, ArchRProj = CARE_filt_rna_myeloid, addDOC = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
### Plot ATAC enrichment
CARE_filt_rna_myeloid <- addArchRAnnotations(ArchRProj = CARE_filt_rna_myeloid, collection = "ATAC", force=TRUE)

enrichATAC <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CARE_filt_rna_myeloid,
  peakAnnotation = "ATAC",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)


heatmapATAC <- plotEnrichHeatmap(enrichATAC, 
                                 n = 10, 
                                 transpose = TRUE,
                                 cutOff = 1.3,
                                 returnMatrix = FALSE)

# The number in parentheses appears to be the max -log10(adj P value)
plotPDF(heatmapATAC, name = "ATAC-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_myeloid, addDOC = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Add motif annotations.
if("Motif" %ni% names(CARE_filt_rna_myeloid@peakAnnotation)){
  CARE_filt_rna_myeloid <- addMotifAnnotations(ArchRProj = CARE_filt_rna_myeloid, motifSet = "cisbp", name = "Motif")
}


## Performing an enrichment of those SCNA regions.
enrichRegions <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CARE_filt_rna_myeloid,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

heatmapRegions <- plotEnrichHeatmap(enrichRegions, 
                                    transpose = TRUE, 
                                    n = 5,
                                    clusterCols= FALSE)

plotPDF(heatmapRegions, name = "Regions-Enriched-Marker-Peak-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_myeloid, addDOC = FALSE)

# Save the project so that the annotations can be more easily accessed in the future. 
saveArchRProject(ArchRProj = CARE_filt_rna_myeloid, outputDirectory = "Save-CAREmut-Myeloid-RNA", load = FALSE, dropCells = FALSE) 


### END ###