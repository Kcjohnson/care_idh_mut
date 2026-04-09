##############################
### Longitudinal peak accessibility differences across malignant states within a patient
### Author: Kevin Johnson
### Updated: 2026.04.08
##############################

## Part 5: Longitudinal per-patient analysis of differential chromatin accessibility per cell state

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

CARE_filt_rna_malignant_peaks_long <- saveArchRProject(ArchRProj = CARE_filt_rna_malignant_peaks, 
                                                       outputDirectory = "Save-CAREmut-Malignant-RNA-Peaks-LongitudinalPeak", 
                                                       load = TRUE) 

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

# Construct a sample-specific identifier per malignant state (e.g., "AC_x_SAMPLEID")
# This will allow us to compare gene accessibility between two timepoints for the same state
CARE_filt_rna_malignant_peaks_long <- addCellColData(ArchRProj = CARE_filt_rna_malignant, data = paste0(CARE_filt_rna_malignant@cellColData$CellStateGroup,"_x_",CARE_filt_rna_malignant@cellColData$Sample), name = "StateBySample", cells = getCellNames(CARE_filt_rna_malignant), force = TRUE)
table(CARE_filt_rna_malignant$StateBySample)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Make pseudo bulk measurements for cell state groups by sample
# Pseudo-bulk refers to a grouping of single cells where the data from each single cell is combined into a single pseudo-sample.
CARE_filt_rna_malignant_peaks_long <- addGroupCoverages(ArchRProj = CARE_filt_rna_malignant, 
                                               groupBy = "StateBySample",  #  AC_x_NL04-1
                                               minCells = 40,  # default
                                               maxCells = 500,  # default
                                               threads = getArchRThreads(),
                                               # Overwrite the data in the ArchRProject object if the pseudo-bulk replicate information already exists
                                               force = TRUE)

## Is macs2 in the path variable?
pathToMacs2 <- findMacs2()

# Iterative overlap peak merging procedure
CARE_filt_rna_malignant_peaks_long <- addReproduciblePeakSet(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "StateBySample", 
  pathToMacs2 = pathToMacs2,
  threads = getArchRThreads(),
)

## Needed to derive marker peaks.
CARE_filt_rna_malignant_peaks_long <- addPeakMatrix(CARE_filt_rna_malignant)

# Save so that it is quicker to run in the future.
CARE_filt_rna_malignant_peaks_long <- saveArchRProject(ArchRProj = CARE_filt_rna_malignant_peaks, 
                                                       outputDirectory = "Save-CAREmut-Malignant-RNA-Peaks-LongitudinalPeak", 
                                                       load = TRUE) 

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
## Add add motif annotations.
if("Motif" %ni% names(CARE_filt_rna_malignant@peakAnnotation)){
  CARE_filt_rna_malignant_peaks_long <- addMotifAnnotations(ArchRProj = CARE_filt_rna_malignant, motifSet = "cisbp", name = "Motif")
}

## Add ENCODE annotations 
CARE_filt_rna_malignant_peaks_long <- addArchRAnnotations(ArchRProj = CARE_filt_rna_malignant, collection = "EncodeTFBS")

### ### ### ### ### ### ### ### ### ### ### 
#Patient specific longitudinal changes
### ### ### ### ### ### ### ### ### ### ### 

## Objective: For each patient and cell state loop through longitudinal comparisons so
# that differential peaks can be determined within patient while controlling for cell state.
table(CARE_filt_rna_malignant$StateBySample)

# Set up for loop run
archr_project_name <- "CARE_filt_rna_malignant"
patient_id <- unique(sapply(strsplit(CARE_filt_rna_malignant@cellColData$Sample, "-"), "[[", 1))
cell_state = unique(CARE_filt_rna_malignant@cellColData$CellStateGroup)
# Create empty vectors and a list to store the results.
up_peaks <- down_peaks <- use_cells <- bgd_cells <- top_up_tfs <- top_down_tfs <- c()
outlist <- list()

atac_md_filt <- as.data.frame(CARE_filt_rna_malignant@cellColData)
table(atac_md_filt$CellStateGroup, atac_md_filt$Sample)

# Check to see whether there will be any cases where there are 3 tumor samples per patient. Manually remove any 3 time point samples with low numbers 
table(CARE_filt_rna_malignant@cellColData$Sample)

# Determine the longitudinal difference in cell states. Since we are interested in more subtle changes, let's use a lower threshold for significant difference.
# Log2FC of 0.5 would represent approximately ~50% increase in accessibility.

for (i in 1:length(patient_id)){
  for(j in 1:length(cell_state)) {
    
    # Get the sample identifier.
    samples <- unique(CARE_filt_rna_malignant@cellColData$Sample)[grep(patient_id[i], unique(CARE_filt_rna_malignant@cellColData$Sample))]
    
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
    
    use_cells = sum(CARE_filt_rna_malignant@cellColData$StateBySample==sprintf("%s", state_use))
    bgd_cells = sum(CARE_filt_rna_malignant@cellColData$StateBySample==sprintf("%s", state_bgd))
    
    # Skip to next iteration if number of cells is too low
    if(use_cells < 40) next
    if(bgd_cells < 40) next
    
    out <- paste0(sprintf("Analyzing patient %s for the comparison: ", patient_id[i]),  sprintf("Use group: %s", state_use), " vs ", sprintf("Background group: %s", state_bgd))
    print(out)
    
    start.time <- Sys.time()
    markerState <- getMarkerFeatures(
      ArchRProj = CARE_filt_rna_malignant, 
      useMatrix = "PeakMatrix",
      groupBy = "StateBySample",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = sprintf("%s", state_use), 
      bgdGroups = sprintf("%s", state_bgd) 
    )
    
    motifsUp <- peakAnnoEnrichment(
      seMarker = markerState,
      ArchRProj = CARE_filt_rna_malignant,
      peakAnnotation = "Motif",
      # Selecting a loosened threshold for these pairwise tests since they are less well powered.
      cutOff = "FDR <= 0.05 & Log2FC >= 0.5"
    )
    
    df_up <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
    df_up <- df_up[order(df_up$mlog10Padj, decreasing = TRUE),]
    df_up$rank <- seq_len(nrow(df_up))
    
    motifsDo <- peakAnnoEnrichment(
      seMarker = markerState,
      ArchRProj = CARE_filt_rna_malignant,
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
    plotPDF(pv, name = sprintf("%s_Marker-Peaks-Volcano", comparison), width = 5, height = 5, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)
    
    results <- data.frame(patient_id[i], cell_state[j], state_use, use_cells, state_bgd, bgd_cells, up_peaks, down_peaks, top_up_tfs, top_down_tfs)
    outlist <- append(outlist, list(results))
  }
}

# Combine all into one data.frame.
patient_list_results <- do.call(rbind, outlist)

sample <- read.table("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/data/metadata/caremut_sample_identifier_linker.txt", sep = "\t", row.names = 1, header = TRUE)
sample_filt <- sample %>% 
  dplyr::select(care_id:SampleID)
patient_md <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240227.txt", sep="\t", header = TRUE)

patient_list_results_filt <- patient_list_results %>% 
  filter(bgd_cells>100, use_cells>100) %>% 
  mutate(SampleID = sapply(strsplit(state_use, "_x_"), "[[", 2)) %>% 
  inner_join(sample_filt, by="SampleID") 

# Write out the results so that they can be used downstream. 
write.table(patient_list_results_filt, "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/archr_care_cell_type_longitudinal_peaks_min100cells_perstate_20240506.txt", sep="\t", row.names = FALSE, col.names = TRUE)
patient_list_results_filt <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/archr_care_cell_type_longitudinal_peaks_min100cells_perstate_20240506.txt", sep="\t", header = TRUE)

res <- patient_list_results_filt %>% 
  left_join(patient_md, by="patient_id")
res$down_peaks <- ifelse(res$down_peaks!=0, res$down_peaks*-1, res$down_peaks)

library(tidyverse)
source("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/scripts/misc/plot_theme.R")

malignant_res <- res %>% 
  dplyr::select(patient_id = patient_id.i., cell_state = cell_state.j., grade_change_t1t2, acquired_genetic_alt_t1t2, idh_codel_subtype, up_peaks, down_peaks) %>% 
  pivot_longer(cols= c(up_peaks:down_peaks),
               names_to = "Type",
               values_to = "Values") 
malignant_res$acquired_genetic_alt_t1t2 <- ifelse(is.na(malignant_res$acquired_genetic_alt_t1t2), "No acquired alt.", "Acquired genetic alt.")

png(paste0(fig_dir, "longitudinal_state_controll_atac_malignant_peaks_codel.png"), width = 8, height = 5, units = 'in', res = 300, bg = "transparent")
ggplot(malignant_res %>% filter(cell_state%in%c("AC", "OPC"), idh_codel_subtype=="IDHmut-codel"), aes(x = patient_id, y = Values, fill = Type)) +
  geom_bar(stat = 'identity', position = 'identity', color = 'black') +
  labs(x = 'Patients', y = 'Longitudinal differential accessibility\nstate-conrolled (IDH-O)', fill="Longitudinal change") +
  scale_fill_manual(values = c('blue', 'red')) +
  plot_theme +
  facet_grid(cell_state~acquired_genetic_alt_t1t2, scales="free") +
  ylim(-10000,10000) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

ggplot(malignant_res %>% filter(cell_state%in%c("AC", "OPC", "Undifferentiated"), patient_id%in%c("SJ04","SJ06", "SJ07", "SJ13", "SJ15")), aes(x = patient_id, y = Values, fill = Type)) +
  geom_bar(stat = 'identity', position = 'identity', color = 'black') +
  labs(x = 'Patients', y = 'Up or down longitudinal ATAC peaks', fill="Longitudinal change") +
  scale_fill_manual(values = c('blue', 'red')) +
  plot_theme +
  facet_grid(cell_state~., scales="free") 

tmp <- res %>% 
  filter(patient_id.i.%in%c("SJ04","SJ06", "SJ07", "SJ13", "SJ15"), cell_state.j.%in%c("AC", "OPC", "Undifferentiated"))

tmp <- malignant_res %>% filter(cell_state%in%c("AC", "OPC", "Undifferentiated"), patient_id%in%c("SJ04","SJ06", "SJ07", "SJ13", "SJ15"))


ggplot(malignant_res, aes(x = patient_id, y = Values, fill = Type)) +
  geom_bar(stat = 'identity', position = 'identity', color = 'black') +
  labs(title = 'Longitudinal differentially accessible peaks (malignant)',
       x = 'Patients', y = 'Up or down longitudinal ATAC peaks', fill="Longitudinal change") +
  scale_fill_manual(values = c('blue', 'red')) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(.~grade_change, scales="free", space="free") 
dev.off()

## Testing
# cutOff	
# A valid-syntax logical statement that defines which marker features from seMarker to use. cutoff can contain any of the assayNames from seMarke

markerState <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "PeakMatrix",
  groupBy = "StateBySample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "OPC_x_SJ13-2", 
  bgdGroups = "OPC_x_SJ13-1",
  
)
markerState2 <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "PeakMatrix",
  groupBy = "StateBySample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "AC_x_SJ13-2", 
  bgdGroups = "AC_x_SJ13-1",
  
)


str(markerState)

volcano_opc <- plotMarkers(seMarker = markerState, name = "OPC_x_SJ04-2", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.5", plotAs = "Volcano")
heatmapGS <- plotMarkerHeatmap(
  seMarker = markerState, 
  cutOff = "FDR <= 0.05 & Log2FC >= ", 
  plotLog2FC = TRUE,
  transpose = TRUE
)
markerList <- getMarkers(markerState, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
markerList2 <- getMarkers(markerState2, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)

#markerList <- getMarkers(markerState, cutOff = "FDR <= 0.05 & Log2FC >= 1")

table(markerList2$`AC_x_SJ13-2`@seqnames@values)
table(markerList$`OPC_x_SJ13-2`@seqnames@values)

intersection <- intersect(markerList$`OPC_x_SJ13-2`, markerList2$`AC_x_SJ13-2`)

table(intersection@seqnames)
chr2_coordinates <- intersection[seqnames(intersection) == "chr2"]

CARE_filt_rna_malignant@geneAnnotation$genes

motifsUp <- peakAnnoEnrichment(
  seMarker = markerState2,
  ArchRProj = CARE_filt_rna_malignant,
  peakAnnotation = "Motif",
  # Selecting a loosened threshold for these pairwise tests since they are less well powered.
  cutOff = "FDR <= 0.05 & Log2FC >= 2"
)


df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))


# Summary: There seems to be an enrichment of Pol2 in samples at recurrence with chromatin change.
# This seems consistent when applied to different cell states from the same patient: AC-T2 vs AC-T1; OPC-T2 vs OPC-T1 yields similar results.
# However, I also noticed that genes with downregulated ATAC peaks are also enriched for these same motifs/TFBS.
# Perhaps, these chromatin changes are poised.

# A few things to check: are the regions the same across cell states (i.e., are AC-T2 vs AC-T1 peaks the same).
# It is concerning that many of the analyses pick up the same thing.
# What about comparing OPC-T2 from one tumor and OPC-T2 from another tumor? These pathways were still enriched.
# Are these differences primarily tracking to the promoter, and thus the enrichment of Pol2/DNMT1?
# I evaluated gene activity scores and proportionally it was not the same. Repeat with all samples? Samples without many changes should show differences.


### END ###