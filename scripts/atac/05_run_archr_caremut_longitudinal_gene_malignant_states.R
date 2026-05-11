##############################
### Assess chromatin gene activity (ATAC) across time points (init. vs recur.) controlled for RNA state
### Author: Kevin Johnson
##############################

## Longitudinal per-patient analysis of differential gene chromatin accessibility per cell state

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

CARE_filt_rna_malignant_long_genes <- saveArchRProject(ArchRProj = CARE_filt_rna_malignant_peaks, 
                                                  outputDirectory = "Save-CAREmut-Malignant-RNA-Peaks-LongitudinalGene", 
                                                  load = TRUE) 

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Construct a sample-specific identifier per malignant state (e.g., "AC_x_SAMPLEID")
# This will allow us to compare gene accessibility between two timepoints for the same state
# See following github discussion for treatment vs. control style anylases: https://github.com/GreenleafLab/ArchR/discussions/696
CARE_filt_rna_malignant_long_genes <- addCellColData(ArchRProj = CARE_filt_rna_malignant_long_genes, data = paste0(CARE_filt_rna_malignant_long_genes@cellColData$CellStateGroup,"_x_",CARE_filt_rna_malignant_long_genes@cellColData$Sample), name = "StateBySample", cells = getCellNames(CARE_filt_rna_malignant_long_genes), force = TRUE)
table(CARE_filt_rna_malignant_long_genes$StateBySample)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Patient specific longitudinal changes
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Objective: For each patient and cell state loop through longitudinal comparisons so
# that differential peaks can be determined within patient while comparing within cell state.

# Set up for loop run
archr_project_name <- "CARE_filt_rna_malignant_long_genes"
patient_id <- unique(sapply(strsplit(CARE_filt_rna_malignant_long_genes@cellColData$Sample, "-"), "[[", 1))
cell_state = unique(CARE_filt_rna_malignant_long_genes@cellColData$CellStateGroup)

# Create empty vectors and a list to store the results.
up_genes <- down_genes <- use_cells <- bgd_cells <- c()
outlist <- list()

atac_md_filt <- as.data.frame(CARE_filt_rna_malignant_long_genes@cellColData)
# There are many comparisons that can't be made because there are too few within a particular category.
table(atac_md_filt$CellStateGroup, atac_md_filt$Sample)

# Check to see whether there will be any cases where there are 3 tumor samples per patient. Manually remove any 3 time point samples with low numbers 
table(CARE_filt_rna_malignant_long_genes@cellColData$Sample) # NL03-1 has very few malignant cells

# Determine the longitudinal difference in cell states. Since we are interested in more subtle changes, let's use a lower threshold for significant difference.
# Log2FC of 1 would represent approximately ~100% increase in accessibility.
set.seed(123)
for (i in 1:length(patient_id)){
  for(j in 1:length(cell_state)) {
    
    # Get the sample identifier.
    samples <- unique(CARE_filt_rna_malignant_long_genes@cellColData$Sample)[grep(patient_id[i], unique(CARE_filt_rna_malignant_long_genes@cellColData$Sample))]
    
    # Remove samples with low cell numbers where the patient has 3 samples. Select max interval for all others.
    samples <- samples[!samples %in% "NL03-1"]
    initial = min(sapply(strsplit(samples, "-"), "[[", 2))
    recurrence = max(sapply(strsplit(samples, "-"), "[[", 2))
    
    # Setting the "use" group will represent results as positive fold change values when peaks are higher in that group.
    state_use = paste0(cell_state[j], "_x_", patient_id[i], "-", recurrence)
    # A negative fold change values when peaks are higher in the background group.
    state_bgd = paste0(cell_state[j], "_x_", patient_id[i], "-", initial)

    out <- paste0(sprintf("Found use group: %s", state_use), " vs ", sprintf("Found background group: %s", state_bgd))
    print(out)
    
    use_cells = sum(CARE_filt_rna_malignant_long_genes@cellColData$StateBySample==sprintf("%s", state_use))
    bgd_cells = sum(CARE_filt_rna_malignant_long_genes@cellColData$StateBySample==sprintf("%s", state_bgd))
    
    # Skip to next iteration if number of cells is too low. We can always filter to a greater number of cells, but we need a minimum.
    if(use_cells < 25) next
    if(bgd_cells < 25) next
    
    out <- paste0(sprintf("Analyzing patient %s for the comparison: ", patient_id[i]),  sprintf("Use group: %s", state_use), " vs ", sprintf("Background group: %s", state_bgd))
    print(out)
    
    start.time <- Sys.time()
    # Note that the maximum cells allowed to be used for a single group are 500 by default. This is okay for these analyses because there are some samples 
    # with many more cells at a given time point.
    markerState <- getMarkerFeatures(
      ArchRProj = CARE_filt_rna_malignant_long_genes, 
      useMatrix = "GeneScoreMatrix",
      groupBy = "StateBySample",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = sprintf("%s", state_use), 
      bgdGroups = sprintf("%s", state_bgd) 
    )
    
    
    upregulated <- getMarkers(markerState, cutOff = "FDR <= 0.05 & Log2FC >= 1")
    upregulated_genes <- as.data.frame(upregulated[[state_use]])
    
    downregulated <- getMarkers(markerState, cutOff = "FDR <= 0.05 & Log2FC <= -1")
    downregulated_genes <-  as.data.frame(downregulated[[state_use]])
    
    dag_df <- bind_rows(upregulated_genes, downregulated_genes)
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(time.taken)
    
    up_genes <- sum(assays(markerState)$Log2FC >= 1 & assays(markerState)$FDR<0.05)
    down_genes <- sum(assays(markerState)$Log2FC <= -1 & assays(markerState)$FDR<0.05)
    
    comparison <- paste(patient_id[i], state_use, state_bgd, sep="_")
    
    #
    print("Writing out results.")
    write.table(dag_df, paste0("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/state_controlled_dag/", comparison, "_longitudinal_dag.txt"), sep="\t", row.names = FALSE, col.names = TRUE)
    
    results <- data.frame(patient_id[i], cell_state[j], state_use, use_cells, state_bgd, bgd_cells, up_genes, down_genes)
    outlist <- append(outlist, list(results))
  }
}

# Combine all into one data.frame and write out results
patient_list_results <- do.call(rbind, outlist)

# Filter out any comparison where there was only one time point, which is not valid.
patient_list_results_filtered <- patient_list_results %>% 
  filter(state_use!=state_bgd)

# Inspect where different thresholds should be applied. While I started off with > 25 cells, I increased the number for this analysis.
patient_list_results_filtered %>% 
  filter(use_cells > 49, bgd_cells > 49)

write.table(patient_list_results_filtered, "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/archr_malignant_cell_state_longitudinal_gene_score_min25cells_perstate.txt", sep="\t", row.names = FALSE, col.names = TRUE)


### END ###