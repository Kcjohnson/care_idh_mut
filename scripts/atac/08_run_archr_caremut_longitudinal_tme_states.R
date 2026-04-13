##############################
### Run ArchR peak calling analysis on CAREmut multiome ATAC data longitudinal analysis of TME cell states
### Author: Kevin Johnson
### Updated: 2026.04.08
##############################

## Part 8: Longitudinal per-patient analysis of differential chromatin accessibility per cell state

# Specify directories
workdir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/"
setwd(workdir)
fig_dir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/figures/archr/"
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/plot_theme.R")

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
CARE_filt_rna_all <- loadArchRProject("Save-CAREmut-All-RNA")

# Remove "Unresolved" (RNA-based) cells from further analysis since it may cause issues. Unresolved are mostly malignant/astrocytes, which seem to be most prone to misclassification in glioma.
atac_df <- data.frame(CARE_filt_rna_all@cellColData)

atac_df_filt <- atac_df %>% 
  filter(CellType_final!="Unresolved")
current_rna_atac_cells <- getCellNames(CARE_filt_rna_all)
rna_cells_to_keep = current_rna_atac_cells[which(current_rna_atac_cells%in%rownames(atac_df_filt))]  

#  Create a new subset for the longitudinal analyese
CARE_filt_rna_all_longitudinal <- subsetArchRProject(ArchRProj = CARE_filt_rna_all, cells = rna_cells_to_keep, outputDirectory = "Save-CAREmut-All-CellTypes-Longitudinal", force = TRUE)
CARE_filt_rna_all_longitudinal <- loadArchRProject("Save-CAREmut-All-CellTypes-Longitudinal")

CARE_filt_rna_all_longitudinal <- addCellColData(ArchRProj = CARE_filt_rna_all_longitudinal, data = paste0(CARE_filt_rna_all_longitudinal@cellColData$CellType_final,"_x_",CARE_filt_rna_all_longitudinal@cellColData$Sample), name = "CellTypeBySample", cells = getCellNames(CARE_filt_rna_all_longitudinal), force = TRUE)
table(CARE_filt_rna_all_longitudinal$CellTypeBySample)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Make pseudo bulk measurements for cell state groups by sample
# Pseudo-bulk refers to a grouping of single cells where the data from each single cell is combined into a single pseudo-sample.
CARE_filt_rna_all_longitudinal <- addGroupCoverages(ArchRProj = CARE_filt_rna_all_longitudinal, 
                                             groupBy = "CellTypeBySample",  #  Oligodendrocyte_x_NL04-1
                                             minCells = 40,  # default
                                             maxCells = 500,  # default
                                             threads = getArchRThreads(),
                                             # Overwrite the data in the ArchRProject object if the pseudo-bulk replicate information already exists
                                             force = TRUE)

## Is macs2 in the path variable?
pathToMacs2 <- findMacs2()

# Iterative overlap peak merging procedure
set.seed(123)
CARE_filt_rna_all_longitudinal <- addReproduciblePeakSet(
  ArchRProj = CARE_filt_rna_all_longitudinal, 
  groupBy = "CellTypeBySample", 
  pathToMacs2 = pathToMacs2,
  threads = getArchRThreads(),
)

# Needed to derive marker peaks.
CARE_filt_rna_all_longitudinal <- addPeakMatrix(CARE_filt_rna_all_longitudinal)

# Add add motif annotations.
if("Motif" %ni% names(CARE_filt_rna_all_longitudinal@peakAnnotation)){
  CARE_filt_rna_all_longitudinal <- addMotifAnnotations(ArchRProj = CARE_filt_rna_all_longitudinal, motifSet = "cisbp", name = "Motif")
}
# Check to see what is available.
getAvailableMatrices(CARE_filt_rna_all_longitudinal)
names(CARE_filt_rna_all_longitudinal@peakAnnotation)

# Save this following peak matrix creation so that it is quicker to run in the future.
CARE_filt_rna_all_longitudinal <- saveArchRProject(ArchRProj = CARE_filt_rna_all_longitudinal, outputDirectory = "Save-CAREmut-All-CellTypes-Longitudinal", load = TRUE)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Patient specific longitudinal changes

# Objective: For each patient and cell state loop through longitudinal comparisons so
# that differential peaks can be determined within patient while controlling for cell state.
table(CARE_filt_rna_all_longitudinal$CellTypeBySample)
table(CARE_filt_rna_all_longitudinal$Sample)

# Set up for loop run.
archr_project_name <- "CARE_filt_rna_all_longitudinal"
patient_id <- unique(sapply(strsplit(CARE_filt_rna_all_longitudinal@cellColData$Sample, "-"), "[[", 1))
cell_state = unique(CARE_filt_rna_all_longitudinal@cellColData$CellType_final)

# Create empty vectors and a list to store the results.
up_peaks <- down_peaks <- use_cells <- bgd_cells <- top_up_tfs <- top_down_tfs <- c()
outlist <- list()

atac_md_filt <- as.data.frame(CARE_filt_rna_all_longitudinal@cellColData)
table(atac_md_filt$CellTypeBySample, atac_md_filt$Sample)

# Check to see whether there will be any cases where there are 3 tumor samples per patient. Manually remove any 3 time point samples with low numbers 
table(CARE_filt_rna_all_longitudinal@cellColData$Sample)

# Determine the longitudinal difference in cell states. Since we are interested in more subtle changes, let's use a lower threshold for significant difference.
# Log2FC of 0.5 would represent approximately ~50% increase in accessibility.
set.seed(123)
for (i in 1:length(patient_id)){
  for(j in 1:length(cell_state)) {
    
    # Get the sample identifier.
    samples <- unique(CARE_filt_rna_all_longitudinal@cellColData$Sample)[grep(patient_id[i], unique(CARE_filt_rna_all_longitudinal@cellColData$Sample))]
    
    # Remove samples where the patient has 3 samples with a preference for removing a sample that has less than 100 cells.
    samples <- samples[!samples%in%c("NL03-1", "SJ03-2", "SJ08-3", "SJ20-3")]
    initial = min(sapply(strsplit(samples, "-"), "[[", 2))
    recurrence = max(sapply(strsplit(samples, "-"), "[[", 2))
    
    # Setting the "use" group will represent results as positive fold change values when peaks are higher in that group.
    state_use = paste0(cell_state[j], "_x_", patient_id[i], "-", recurrence)
    # A negative fold change values when peaks are higher in the background group.
    state_bgd = paste0(cell_state[j], "_x_", patient_id[i], "-", initial)
    
    out <- paste0(sprintf("Found use group: %s", state_use), " vs ", sprintf("Found background group: %s", state_bgd))
    print(out)
    
    use_cells = sum(CARE_filt_rna_all_longitudinal@cellColData$CellTypeBySample==sprintf("%s", state_use))
    bgd_cells = sum(CARE_filt_rna_all_longitudinal@cellColData$CellTypeBySample==sprintf("%s", state_bgd))
    
    # Skip to next iteration if number of cells is too low
    if(use_cells < 40) next
    if(bgd_cells < 40) next
    
    out <- paste0(sprintf("Analyzing patient %s for the comparison: ", patient_id[i]),  sprintf("Use group: %s", state_use), " vs ", sprintf("Background group: %s", state_bgd))
    print(out)
    
    start.time <- Sys.time()
    markerState <- getMarkerFeatures(
      ArchRProj = CARE_filt_rna_all_longitudinal, 
      useMatrix = "PeakMatrix",
      groupBy = "CellTypeBySample",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = sprintf("%s", state_use), 
      bgdGroups = sprintf("%s", state_bgd) 
    )
    
    motifsUp <- peakAnnoEnrichment(
      seMarker = markerState,
      ArchRProj = CARE_filt_rna_all_longitudinal,
      peakAnnotation = "Motif",
      # Selecting a loosened threshold for these pairwise tests since they are less well powered due to cell number.
      cutOff = "FDR <= 0.05 & Log2FC >= 0.5"
    )
    
    df_up <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
    df_up <- df_up[order(df_up$mlog10Padj, decreasing = TRUE),]
    df_up$rank <- seq_len(nrow(df_up))
    
    motifsDo <- peakAnnoEnrichment(
      seMarker = markerState,
      ArchRProj = CARE_filt_rna_all_longitudinal,
      peakAnnotation = "Motif",
      cutOff = "FDR <= 0.05 & Log2FC <= -0.5"
    )
    
    df_do <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
    df_do <- df_do[order(df_do$mlog10Padj, decreasing = TRUE),]
    df_do$rank <- seq_len(nrow(df_do))
    
    # Reporting the top TFs up and down for a simple summary.
    top_up_tfs <- paste(df_up$TF[1:10], sep="", collapse=";") 
    top_down_tfs <- paste(df_do$TF[1:10], sep="", collapse=";") 
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
    
    up_peaks <- sum(assays(markerState)$Log2FC >= 0.5 & assays(markerState)$FDR<0.05)
    down_peaks <- sum(assays(markerState)$Log2FC <= -0.5 & assays(markerState)$FDR<0.05)
    
    comparison <- paste(patient_id[i], state_use, state_bgd, sep="_")
    
    pv <- plotMarkers(seMarker = markerState, name = sprintf("%s", state_use), cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
    plotPDF(pv, name = sprintf("%s_Marker-Peaks-Volcano", comparison), width = 5, height = 5, ArchRProj = CARE_filt_rna_all_longitudinal, addDOC = FALSE)
    
    results <- data.frame(patient_id[i], cell_state[j], state_use, use_cells, state_bgd, bgd_cells, up_peaks, down_peaks, top_up_tfs, top_down_tfs)
    outlist <- append(outlist, list(results))
  }
}

# Combine all into one data.frame.
patient_list_results <- do.call(rbind, outlist)

# Be sure to remove any samples that were used as both the case/control groups.
patient_list_results_filtered <- patient_list_results %>% 
  mutate(sample_use = sapply(strsplit(state_use, "_x_"), "[[", 2),
         sample_bgd = sapply(strsplit(state_bgd, "_x_"), "[[", 2)) %>% 
# Make sure no SampleID was used twice due to absence of multiple samples from the same patient.
  filter(sample_use!=sample_bgd) 


# Write out the results so that they can be used downstream. 
write.table(patient_list_results_filtered, "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/archr_care_tme_celltype_longitudinal_peaks_min40cells_per_celltype.txt", sep="\t", row.names = FALSE, col.names = TRUE)


### END ###