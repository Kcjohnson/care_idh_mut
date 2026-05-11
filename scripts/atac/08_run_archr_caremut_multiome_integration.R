##############################
### Run ArchR RNA and ATAC Multiome integration for peak-to-gene linkage
### Author: Kevin Johnson
##############################

# Create a visualization for peak to gene linkage (RNA+ATAC)

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

# This project was created in "04_run_archr_caremut_peakcall_malignant_states.R".
CARE_filt_rna_malignant_peaks <- loadArchRProject("Save-CAREmut-Malignant-RNA-Peaks")

CARE_filt_rna_malignant_multiome <- saveArchRProject(ArchRProj = CARE_filt_rna_malignant_peaks, 
                                                       outputDirectory = "Save-CAREmut-Malignant-RNA-Peaks-Multiome", 
                                                       load = TRUE) 

# Specify the RNA files to be loaded in - these will be the same as the ATAC libraries but selected for RNA.
# These are found in the same cellranger-arc directories as the atac fragment files
#                # NL01
rna_files <- c("/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL01-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL01-2/outs/filtered_feature_bc_matrix.h5",
                 # NL03
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL03-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL03-2/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL03-3/outs/filtered_feature_bc_matrix.h5",
                 # NL04
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL04-0/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL04-1/outs/filtered_feature_bc_matrix.h5",
                 # NL05
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL05-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL05-2/outs/filtered_feature_bc_matrix.h5",
                 # SN05
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/SN05-0/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/SN05-1/outs/filtered_feature_bc_matrix.h5",
                 # SN07
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/SN07-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/SN07-2/outs/filtered_feature_bc_matrix.h5",
                 # SN17
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/SN17-0/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/SN17-1/outs/filtered_feature_bc_matrix.h5",
                 # NL11
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL11-0/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL11-1/outs/filtered_feature_bc_matrix.h5",
                 # NL12
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL12-0/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL12-1/outs/filtered_feature_bc_matrix.h5",
                 # NL23
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL23-0/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL23-1/outs/filtered_feature_bc_matrix.h5",
                 # NL26
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL26-0/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL26-1/outs/filtered_feature_bc_matrix.h5",
                 # SJ03
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ03-0/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ03-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ03-2/outs/filtered_feature_bc_matrix.h5",
                 # SJ04
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ04-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ04-2/outs/filtered_feature_bc_matrix.h5",
                 # SJ06
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ06-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ06-2/outs/filtered_feature_bc_matrix.h5",
                 # SJ07
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ07-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ07-2/outs/filtered_feature_bc_matrix.h5",
                 # SJ08
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ08-0/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ08-2/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ08-3/outs/filtered_feature_bc_matrix.h5",
                 # SJ10
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ10-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ10-2/outs/filtered_feature_bc_matrix.h5",
                 # SJ12
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ12-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ12-2/outs/filtered_feature_bc_matrix.h5",
                 # SJ13
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ13-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ13-2/outs/filtered_feature_bc_matrix.h5",
                 # SJ15
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ15-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ15-2/outs/filtered_feature_bc_matrix.h5",
                 # SJ17
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ17-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ17-2/outs/filtered_feature_bc_matrix.h5",
                 # SJ20
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ20-1/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ20-2/outs/filtered_feature_bc_matrix.h5",
                 "/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ20-3/outs/filtered_feature_bc_matrix.h5")

# Use the shortened multiome IDs for the sample names.
sample_names <- basename(dirname(dirname(rna_files)))
names(rna_files) <- sample_names

# 1. Import each RNA feature matrix
rna_list <- lapply(names(rna_files), function(smp) {
  import10xFeatureMatrix(
    input = rna_files[[smp]],
    names = smp,
    strictMatch = FALSE
  )
})
names(rna_list) <- names(rna_files)

# 2. Combine them
# You can either add one at a time, or merge into one combined object
# If you want a single seRNA object:
common_genes <- Reduce(intersect, lapply(rna_list, rownames))
rna_counts <- lapply(rna_list, function(se) assay(se, "counts")[common_genes, , drop = FALSE])
rna_combined <- do.call(cbind, rna_counts)

# Use shared gene order and metadata
gene_metadata <- rowData(rna_list[[1]])[common_genes, , drop = FALSE]

# Construct a clean SummarizedExperiment
rna_combined_se <- SummarizedExperiment(
  assays = list(counts = rna_combined),
  rowData = gene_metadata
)

# Quick sanity check
dim(rna_combined_se)
length(rownames(rna_combined_se))
length(colnames(rna_combined_se))

# Use gene ranges from one of the imported 10x objects
rowRanges(rna_combined_se) <- rowRanges(rna_list[[1]])[rownames(rna_combined_se)]

# Most of the cells were also recovered from RNA as well in the cellranger-arc pipeline (99.2%, 70,835 out of 71,365). This is mentioned as it's different from the CellRanger 6.1.2 used for RNA-only analyses.
cellsToKeep <- which(getCellNames(CARE_filt_rna_malignant_multiome) %in% colnames(rna_combined_se))
cellsSample <- getCellNames(CARE_filt_rna_malignant_multiome)[cellsToKeep]

# The cells seem to be distributed across many samples and the missing cells are therefore not some other error.
cellsMissing <- which(!getCellNames(CARE_filt_rna_malignant_multiome) %in% colnames(rna_combined_se))
getCellNames(CARE_filt_rna_malignant_multiome)[cellsMissing]

# Create a subset of the ArchR object to save malignant-only multiome analyses.
projMulti <- subsetArchRProject(ArchRProj = CARE_filt_rna_malignant_multiome, cells = cellsSample, outputDirectory = "Save-CAREmut-Malignant-Multiome-Analysis", force = TRUE)

# Load the ArchR project restricted to malignant cells that passed the cellranger-arc pipeline
projMulti <- loadArchRProject("Save-CAREmut-Malignant-Multiome-Analysis")

# Add the gene expression matrix to this ArchR project
projMulti2 <- addGeneExpressionMatrix(
  input = projMulti,               
  seRNA = rna_combined_se,
  force = TRUE
)

# Confirming that the cell state grouping is properly recorded.
table(projMulti2$CellStateGroup)

# Check to make sure that the expected matrices have already been determined
getAvailableMatrices(projMulti2) 

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# You can either use the TileMatrix for ATAC of GeneScoreMatrix
projMulti2 <- addIterativeLSI(
  ArchRProj = projMulti2, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC"
)

# Repeat for RNA 
projMulti2 <- addIterativeLSI(
  ArchRProj = projMulti2, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)


projMulti2 <- addCombinedDims(projMulti2, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.4, force = TRUE)
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 0.4, force = TRUE)
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_Combined", name = "Clusters_Combined", resolution = 0.4, force = TRUE)

cols <- c("#AA2756", "#F77D58", "#7fbf7b", "#E8F5A3", "gray90")
names(cols)  <- names(table(projMulti2$CellStateGroup)) 

p1 <- plotEmbedding(projMulti2, name = "CellStateGroup", embedding = "UMAP_ATAC", size = 1, labelAsFactors=F, labelMeans=F, pal = cols)
p2 <- plotEmbedding(projMulti2, name = "CellStateGroup", embedding = "UMAP_RNA", size = 1, labelAsFactors=F, labelMeans=F, pal = cols)
p3 <- plotEmbedding(projMulti2, name = "CellStateGroup", embedding = "UMAP_Combined", size = 1, labelAsFactors=F, labelMeans=F, pal = cols)

p <- lapply(list(p1,p2,p3), function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

# The combined RNA+ARAC doesn't change the UMAP visualization much. There seems to be patient-specific clusters (likely due to CNVs) and separation across malignant states.
do.call(cowplot::plot_grid, c(list(ncol = 3),p))


# We don't need to re-run these steps because we've already run these steps in a prior script.
# pathToMacs2 <- findMacs2()
# projMulti2 <- addGroupCoverages(ArchRProj = projMulti2, groupBy = "CellStateGroup", verbose = FALSE, force=TRUE)
# projMulti2 <- addReproduciblePeakSet(ArchRProj = projMulti2, groupBy = "CellStateGroup", pathToMacs2 = pathToMacs2, force=TRUE)
# projMulti2 <- addPeakMatrix(ArchRProj = projMulti2, force=TRUE)

projMulti2 <- addPeak2GeneLinks(ArchRProj = projMulti2, reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix")

# Extract some peak2gene links
p2g <- getPeak2GeneLinks(
  ArchRProj = projMulti2,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = TRUE
)

p2g[[1]]

markerGenes  <- c(
  "AQP4", # AC-like
  "PDGFRA", # OPC-like
  "OLIG1", # OPC-like
  "HOXD11", # Undifferentiated
  "CD44", # AC/MES-like
  "VIM", # AC/MES-like
  "ANXA2", # MES-like
  "DLL3") # OPC-like

p <- plotBrowserTrack(
  ArchRProj = projMulti2, 
  groupBy = "CellStateGroup", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(projMulti2)
)
grid::grid.newpage()
grid::grid.draw(p$AQP4)

grid::grid.newpage()
grid::grid.draw(p$OLIG1)

grid::grid.newpage()
grid::grid.draw(p$PDGFRA)

grid::grid.newpage()
grid::grid.draw(p$CD44)

grid::grid.newpage()
grid::grid.draw(p$VIM)

grid::grid.newpage()
grid::grid.draw(p$ANXA2)

grid::grid.newpage()
grid::grid.draw(p$DLL3)

# Print these to the ArchR project
plotPDF(plotList = p, 
        name = "Plot-Tracks-Key-Genes-with-Peak2GeneLinks.pdf", 
        ArchRProj = projMulti2, 
        addDOC = FALSE, width = 5, height = 5)


# Re-run with returnLoops since this causes errors for some reason downstream
p2g <- getPeak2GeneLinks(
  ArchRProj = projMulti2,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

# Getting both the gene and peak name
p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
p2g$peakName <- (metadata(p2g)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g$idxATAC]

# Confirm that it's working as expected
p2g
metadata(p2g)$seRNA

# While the major malignant cell states  are mostly separated in the hierarchical clustering, the exact clustering is pretty dependent on selection of k and corCutOff (as expected)
p2g_heat <- plotPeak2GeneHeatmap(ArchRProj = projMulti2, 
                               groupBy = "CellStateGroup",
                               palGroup=cols,
                               k = 10,
                               corCutOff = 0.45,             
                               varCutOffATAC = 0.25,
                               varCutOffRNA = 0.25,
                               nPlot = 25000)

p2g_heat

p2g_heat_patient <- plotPeak2GeneHeatmap(ArchRProj = projMulti2, 
                               groupBy = "patient_id",
                               k = 10,
                               corCutOff = 0.45,             
                               varCutOffATAC = 0.25,
                               varCutOffRNA = 0.25,
                               nPlot = 25000)

p2g_heat_patient

plotPDF(p2g_heat, name = "Plot-Peak2GeneLinks-Heatmap.pdf", ArchRProj = projMulti2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p2g_heat_patient, name = "Plot-Peak2GeneLinks-Heatmap-Patient.pdf", ArchRProj = projMulti2, addDOC = FALSE, width = 5, height = 5)


p2g_heat_k8 <- plotPeak2GeneHeatmap(ArchRProj = projMulti2, 
                                    groupBy = "CellStateGroup",
                                    palGroup=cols,
                                    k = 8,
                                    corCutOff = 0.5,             
                                    varCutOffATAC = 0.25,
                                    varCutOffRNA = 0.25,
                                    nPlot = 15000)

p2g_heat_k8

p2g_heat_patient_k8 <- plotPeak2GeneHeatmap(ArchRProj = projMulti2, 
                                    groupBy = "patient_id",
                                    k = 8,
                                    corCutOff = 0.5,             
                                    varCutOffATAC = 0.25,
                                    varCutOffRNA = 0.25,
                                    nPlot = 15000)

p2g_heat_patient_k8

# Versions to include as a supplementary figure
plotPDF(p2g_heat_k8, name = "Plot-Peak2GeneLinks-Heatmap-State-K8-nplot15000.pdf", ArchRProj = projMulti2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p2g_heat_patient_k8, name = "Plot-Peak2GeneLinks-Heatmap-Patient-K8-nplot15000.pdf", ArchRProj = projMulti2, addDOC = FALSE, width = 5, height = 5)

# Check point to save.
projMulti2 <- saveArchRProject(ArchRProj = projMulti2, outputDirectory = "Save-ArchR-Multiome-Analysis", overwrite = TRUE, load = TRUE)

### END ###