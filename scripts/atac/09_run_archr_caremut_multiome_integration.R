##############################
### Run ArchR RNA and ATAC Multiome integration
### Author: Kevin Johnson
### Updated: 2026.04.08
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

## Specify the RNA files to be loaded in - these will be the same as the ATAC libraries but selected for RNA
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

# 2. Combine them (optional — see below)
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

# Most of the cells are there from RNA (99.2%, 70,835 out of 71,365).
cellsToKeep <- which(getCellNames(CARE_filt_rna_malignant_multiome) %in% colnames(rna_combined_se))
cellsSample <- getCellNames(CARE_filt_rna_malignant_multiome)[cellsToKeep]

# The cells seem to be distributed acros many samples and the missing cells are therefore not a mapping error.
cellsMissing <- which(!getCellNames(CARE_filt_rna_malignant_multiome) %in% colnames(rna_combined_se))
getCellNames(CARE_filt_rna_malignant_multiome)[cellsMissing]

# Create a subset of the ArchR object to save malignant-only analyses.
projMulti <- subsetArchRProject(ArchRProj = CARE_filt_rna_malignant_multiome, cells = cellsSample, outputDirectory = "Save-CAREmut-Malignant-Multiome-Analysis", force = TRUE)

# Load the ArchR project restricted to malignant cells that passed the ARC pipeline
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
# Not sure what happened with copying the arrow files to this new project but since this returns "Not All Seqnames Identical", I am moving to "GeneScoreMatrix" approach.
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

#projMulti2 <- addIterativeLSI(
#  ArchRProj = projMulti2,
#  useMatrix = "GeneScoreMatrix", 
#  name = "LSI_ATAC", 
#  iterations = 2, 
#  clusterParams = list( 
#    resolution = c(0.2), 
#    sampleCells = 20000, 
#    n.start = 10
#  ), 
#  varFeatures = 2500, 
#  dimsToUse = 1:20, 
#  force = TRUE
#)

# Repeat for RNA - may go to a smaller number of features since it likely picks up sample specific chromatin/expression features at a point

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

# The combined doesn't change the visualization much.
do.call(cowplot::plot_grid, c(list(ncol = 3),p))


# We don't need to re-run these steps because we've already run these steps in a prior script.
# pathToMacs2 <- findMacs2()
# projMulti2 <- addGroupCoverages(ArchRProj = projMulti2, groupBy = "CellStateGroup", verbose = FALSE, force=TRUE)
# projMulti2 <- addReproduciblePeakSet(ArchRProj = projMulti2, groupBy = "CellStateGroup", pathToMacs2 = pathToMacs2, force=TRUE)
# projMulti2 <- addPeakMatrix(ArchRProj = projMulti2, force=TRUE)

# It says "unused force=TRUE", which is odd because I thought this was used in that last one.
# projMulti2 <- addPeak2GeneLinks(ArchRProj = projMulti2, reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix", force=TRUE)
projMulti2 <- addPeak2GeneLinks(ArchRProj = projMulti2, reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix")

# Inspect some of the markers:
se <- getMarkerFeatures(ArchRProj = projMulti2,
                        groupBy = "CellStateGroup",
                        bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"))

heatmap_gex <- plotMarkerHeatmap(
  seMarker = se, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  nLabel = 4,
  transpose = TRUE
)

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
  "OLIG1",
  "HOXD11",
  "CD44",
  "VIM",
  "ANXA2",
  "DLL3")

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

# This seems to be about the best one could hope for in terms of separation
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
                               nPlot = 20000)

p2g_heat_patient

plotPDF(p2g_heat, name = "Plot-Peak2GeneLinks-Heatmap.pdf", ArchRProj = projMulti2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p2g_heat_patient, name = "Plot-Peak2GeneLinks-Heatmap-Patient.pdf", ArchRProj = projMulti2, addDOC = FALSE, width = 5, height = 5)


# Check point to save
projMulti2 <- saveArchRProject(ArchRProj = projMulti2, outputDirectory = "Save-ArchR-Multiome-Analysis", overwrite = TRUE, load = TRUE)

########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~############
#
########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~############













# The analysis below aims to determine peak2gene relationships for each cell state.
# Need to store the results following each cell state and timepoint.
cell_states <- unique(projMulti2$CellStateGroup)
timepoints  <- c("T1", "T2")

# Keep per-subset results
p2g_links_by_subset <- list()

for (state in cell_states) {
  for (tp in timepoints) {
    
    message("P2G for ", state, " @ ", tp)
    
    cells_subset <- getCellNames(projMulti2)[
      projMulti2$CellStateGroup == state & projMulti2$timepoint == tp
    ]
    
    if (length(cells_subset) < 100) {
      message("Skipping ", state, " ", tp, " (", length(cells_subset), " cells)")
      next
    }
    
    # compute links within this subset
    projMulti2 <- addPeak2GeneLinks(
      ArchRProj  = projMulti2,
      reducedDims = "LSI_Combined", 
      useMatrix = "GeneExpressionMatrix",
      cellsToUse  = cells_subset,
      verbose     = TRUE,
      logFile     = paste0("ArchRLogs/P2G_", state, "_", tp, ".log")
    )
    
    
    # immediately retrieve and stash the links produced by the call above
    p2g_now <- getPeak2GeneLinks(
      ArchRProj     = projMulti2,
      corCutOff     = 0.0,        # pull all; filter later
      varCutOffATAC = 0.0,
      varCutOffRNA  = 0.0,
      returnLoops   = FALSE
    )
    
    p2g_now$GeneSymbol <- mcols(metadata(p2g_now)$geneSet)$name[p2g_now$idxRNA]
    p2g_now$Peak <- (metadata(p2g_now)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g_now$idxATAC]
    
    if (!is.null(p2g_now) && nrow(as.data.frame(p2g_now)) > 0) {
      df <- as.data.frame(p2g_now) %>%
        transmute(
          peakName = Peak,
          geneName = GeneSymbol,
          Correlation = Correlation,
          FDR = FDR
        ) %>%
        mutate(cell_state = state, timepoint = tp)
      key <- paste(state, tp, sep = "_")
      p2g_links_by_subset[[key]] <- df
    }
  }
}

### Compare T2 vs T1 within each cell state ### 
p2g_diff_all <- bind_rows(p2g_links_by_subset)

p2g_diff_list <- lapply(cell_states, function(state) {
  t1 <- p2g_links_by_subset[[paste0(state, "_T1")]]
  t2 <- p2g_links_by_subset[[paste0(state, "_T2")]]
  if (is.null(t1) || is.null(t2)) return(NULL)
  
  full_join(t1, t2, by = c("peakName", "geneName"),
            suffix = c("_T1", "_T2")) %>%
    mutate(
      corChange = Correlation_T2 - Correlation_T1,
      direction = case_when(
        # Peak present only at T2 → new / gained accessibility
        is.na(Correlation_T1) & !is.na(Correlation_T2) ~ "New Peak / Gain",
        
        # Peak present only at T1 → lost accessibility
        !is.na(Correlation_T1) & is.na(Correlation_T2) ~ "Lost Peak / Loss",
        
        # Peak present in both, correlation increased
        !is.na(corChange) & corChange >  0.25 ~ "Increased",
        
        # Peak present in both, correlation decreased
        !is.na(corChange) & corChange < -0.25 ~ "Decreased",
        
        # Otherwise no meaningful change
        TRUE ~ "Stable"
      ),
      cell_state = state
    )
})

p2g_diff_all <- bind_rows(p2g_diff_list)

# Store output
saveRDS(p2g_links_by_subset, "/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/longitudinal_malignant_p2g_links_by_subset.RDS")
saveRDS(p2g_diff_all, "/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/longitudinal_malignant_p2g_differenences.RDS")

# Load back in since this took quite some time.
p2g_links_by_subset <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/longitudinal_malignant_p2g_links_by_subset.RDS")
p2g_diff_all <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/longitudinal_malignant_p2g_differenences.RDS")

# Setting a higher bar for correlation strength since we want to garner information about only the most relevant chromatin-gene relationships.
p2g_diff_all_filtered <- p2g_diff_all %>% 
  filter(Correlation_T1>0.5 | Correlation_T2 > 0.5) 

# Produce a summary for plotting
p2g_summary <- p2g_diff_all_filtered %>%
  dplyr::group_by(cell_state, direction) %>%
  dplyr::summarise(n_links = n(), .groups = "drop") %>%
  dplyr::group_by(cell_state) %>%
  dplyr::mutate(pct_links = 100 * n_links / sum(n_links))

p2g_summary

ggplot(p2g_summary, aes(x=cell_state, y=pct_links, fill=direction)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=c(
    "Increased"="steelblue3",
    "Decreased"="firebrick3",
    "New Peak / Gain"="#2171b5",
    "Lost Peak / Loss"="#a50f15",
    "Stable"="grey80"
  )) +
  labs(x="Cell State", y="Fraction of P2G Links", fill="Direction")

### ### ### ### ### ### ### ### ### ###
# Read in the differentially expressed genes from longitudinal pseudobulk analyses
### ### ### ### ### ### ### ### ### ###

library(dplyr)
library(readr)
library(purrr)
library(stringr)

# Define the files where pseudobulk results were stored
files <- c(
  "/vast/palmer/pi/verhaak/kcj28/care_mut/results/pseudobulk/pseudobulk_ac_t1t2_deg.txt",
  "/vast/palmer/pi/verhaak/kcj28/care_mut/results/pseudobulk/pseudobulk_opc_t1t2_deg.txt",
  "/vast/palmer/pi/verhaak/kcj28/care_mut/results/pseudobulk/pseudobulk_npc_t1t2_deg.txt",
  "/vast/palmer/pi/verhaak/kcj28/care_mut/results/pseudobulk/pseudobulk_mes_t1t2_deg.txt",
  "/vast/palmer/pi/verhaak/kcj28/care_mut/results/pseudobulk/pseudobulk_undiff_t1t2_deg.txt"
)

# Helper to extract cell type from filename
get_celltype <- function(path) {
  str_match(basename(path), "pseudobulk_(.*?)_t1t2_deg")[,2]
}

# Define the mapping  between the different naming schema
celltype_map <- c(
  "ac" = "AC",
  "opc" = "OPC",
  "npc" = "NPC",
  "mes" = "MES",
  "undiff" = "Undifferentiated"
)

# Read and combine all pseudobulk DEG tables
pseudobulk_all <- map_dfr(files, function(f) {
  raw_type <- get_celltype(f)
  read.delim(f, check.names = FALSE) %>%
    mutate(
      cell_type_raw = raw_type,
      cell_type = celltype_map[raw_type],
      source_file = basename(f)
    )
})

# Quick sanity check
table(pseudobulk_all$cell_type)

### Examine overlap and test for enrichment ###

up_genes_by_state <- pseudobulk_all %>%
  filter(log2FoldChange > 0, padj < 0.1) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    up_genes = list(unique(feature)),
    .groups = "drop"
  )

up_genes_by_state

p2g_collapsed <- p2g_diff_all_filtered %>%
  dplyr::group_by(cell_state, geneName) %>%
  dplyr::summarise(
    direction = names(sort(table(direction), decreasing = TRUE))[1],
    .groups = "drop"
  )


overlap_summary <- p2g_collapsed %>%
  dplyr::group_by(cell_state, direction) %>%
  dplyr::summarise(
    n_total = n(),
    n_overlap = sum(geneName %in% unlist(
      up_genes_by_state$up_genes[up_genes_by_state$cell_type == unique(cell_state)]
    )),
    prop_overlap = n_overlap / n_total,
    .groups = "drop"
  ) %>%
  arrange(cell_state, prop_overlap)

ggplot(overlap_summary, aes(x = reorder(direction, prop_overlap), y = prop_overlap, fill = direction)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~cell_state, scales = "free_y") +
  coord_flip() +
  labs(
    x = "P2G grouping",
    y = "Proportion overlapping with T2 upregulated genes",
    title = "Overlap between chromatin change groupings and T2 upregulated DEGs"
  ) +
  theme_minimal()


deg_long <- pseudobulk_all %>%
  filter(!is.na(log2FoldChange), !is.na(padj), group=="T2") %>%
  mutate(deg_direction = ifelse(log2FoldChange > 0, "Up_T2", "Down_T2"),
         state_gene = paste0(cell_type, "_", feature),
         deg_status = dplyr::case_when(
           deg_direction == "Up_T2" & padj<0.1 ~ "Upregulated",
           deg_direction == "Down_T2" & padj<0.1 ~ "Downregulated",
             padj > 0.1 ~ "Stable"))

p2g_collapsed_merge <- p2g_collapsed %>% 
  mutate(state_gene = paste0(cell_state, "_", geneName)) 

deg_long_p2g <- deg_long %>% 
  left_join(p2g_collapsed_merge, by="state_gene") %>% 
  mutate(p2g_direction = ifelse(is.na(direction), "no_peak_to_gene", direction))

cor.test(df_tierB$Correlation_T1, df_tierB$log2FoldChange, method = "spearman")


table(deg_long_p2g$deg_direction, deg_long_p2g$p2g_direction)

p2g_summary <- deg_long_p2g %>%
  dplyr::group_by(deg_status, p2g_direction) %>%
  dplyr::summarise(counts = n()) %>% 
  dplyr::mutate(prop = counts / sum(counts)) %>%
  ungroup()


ggplot(p2g_summary, aes(x = deg_status, y = prop, fill = p2g_direction)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(
    "Increased" = "#1b9e77",
    "Decreased" = "#d95f02",
    "Stable" = "#7570b3",
    "Lost Peak / Loss" = "#e7298a",
    "New Peak / Gain" = "#66a61e",
    "no_peak_to_gene" = "grey70"
  )) +
  labs(
    x = "Differential Expression Category",
    y = "Proportion of genes",
    fill = "Peak-to-Gene Change",
    title = "Association between differential expression and chromatin linkage categories"
  ) +
  theme_minimal(base_size = 13) 

tab <- p2g_summary %>%
  select(deg_status, p2g_direction, counts) %>%
  tidyr::pivot_wider(
    names_from = p2g_direction,
    values_from = counts,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("deg_status") %>%
  as.matrix()

chisq.test(tab)


fisher_results <- p2g_summary %>%
  dplyr::group_by(p2g_direction) %>%
  dplyr::group_modify(~{
    mat <- matrix(
      c(
        .x$counts[.x$deg_status == "Upregulated"],
        sum(.x$counts[.x$deg_status == "Upregulated"]),
        .x$counts[.x$deg_status == "Stable"],
        sum(.x$counts[.x$deg_status == "Stable"])
      ),
      nrow = 2
    )
    test <- fisher.test(mat)
    tibble(
      direction = unique(.x$p2g_direction),
      pval = test$p.value,
      odds_ratio = test$estimate
    )
  }) %>%
  mutate(padj = p.adjust(pval, method = "BH"))

### END ###