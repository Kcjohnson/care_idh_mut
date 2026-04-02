##############################
### Run ArchR peak calling analysis on CAREmut multiome ATAC data longitudinal analysis of TME cell states
### Author: Kevin Johnson
### Updated: 2025.10.14
##############################

## Part 9: Longitudinal per-patient analysis of differential chromatin accessibility per cell state

## ArchR creates several directories automatically when creating arrow files and ArchR projects.
workdir <- "/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/"
setwd(workdir)

source("/vast/palmer/pi/verhaak/kcj28/care_mut/scripts/misc/plot_theme.R")

## Load necessary packages.
library(dplyr)
library(ArchR)
library(parallel)
library(pheatmap)
library(chromVARmotifs)
library(BSgenome.Hsapiens.UCSC.hg38)

## Specify output directory to drop figures:
fig_dir <- "/vast/palmer/pi/verhaak/kcj28/care_mut/results/figures/archr/"

#### Set-up #####
## Check available cores. ArchR recommends setting total number of cores 1/2 to 3/4 of available cores,
num_cores <- detectCores() # e.g., 36
n_threads <- num_cores/2
addArchRThreads(threads = n_threads) 
## Each R session requires that the genome is also specified and must match alignment.
addArchRGenome("hg38")

# This larger dataset contains both IDHmut and IDHwt tumors that we've processed. Need to restrict to IDHmut-only tumors.
projMONITOR <- loadArchRProject("Save-AllSamples-2024")

projMONITOR_filt <- filterDoublets(projMONITOR)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# CAREmut data that was processed by Seurat and cell types were assigned based on RNA.
mut_md <- read.table("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20240212.txt", sep = "\t", header = TRUE)

# For this analysis, only ATAC profiles were generated for samples processed by the Verhaak lab and defined by RNA as malignant cells
mut_md_verhaak <- mut_md %>% 
  filter(lab=="Verhaak lab", !CellType_final%in%c("Unresolved","Malignant"))

## Create a data.frame that contains the essential information to be merged and inspected.
atac_df <- data.frame(projMONITOR_filt$cellNames, projMONITOR_filt@cellColData$Sample, projMONITOR_filt@cellColData$TSSEnrichment)

# The cell names are a little different between RNA (Seurat) and ATAC (ArchR). Need to create a common linker.
atac_df$CellID <- gsub("#", "-", atac_df$projMONITOR_filt.cellNames)
# 46,155 cells present in ATACseq data that pass doublet removal and are also found in snRNAseq data that passed QC for IDH-mutant.
sum(atac_df$CellID%in%mut_md_verhaak$CellID)

# Combine the two data.frames by adding on the RNA data and filtering out what's leftover.
atac_df_filt_rna <- atac_df %>% 
  inner_join(mut_md_verhaak, by="CellID") 

## Identify which cells to keep in the analysis. 
tmp <- getCellNames(projMONITOR_filt)
rna_cells = tmp[which(tmp%in%atac_df_filt_rna$projMONITOR_filt.cellNames)]  

## Subset the cells to only those that also have RNA.
CARE_filt_rna_tme <- subsetCells(ArchRProj = projMONITOR_filt, cellNames = rna_cells)


myeloid_md <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/processed_data/rna/tme/caremut_myeloid_care_pval_classification_20240303.txt", header = TRUE)
myeloid_md_trim <- myeloid_md %>% 
  mutate(cell_state = recode(cell_state, `Unresolved` = "TAM_Unresolved")) %>% 
  dplyr::select(CellID, cell_state)

mut_md_verhaak_annot <- mut_md_verhaak %>% 
  left_join(myeloid_md_trim, by="CellID") %>% 
  mutate(CellTypeArchR = CellType_final)
mut_md_verhaak_annot$CellTypeArchR[mut_md_verhaak_annot$CellTypeArchR=="Myeloid"] <- mut_md_verhaak_annot$cell_state[mut_md_verhaak_annot$CellTypeArchR=="Myeloid"]

table(mut_md_verhaak_annot$CellTypeArchR, mut_md_verhaak_annot$CellType_final)
table(mut_md_verhaak_annot$CellTypeArchR, mut_md_verhaak_annot$care_id)

CARE_filt_rna_tme$CellID <- gsub("#", "-", CARE_filt_rna_tme$cellNames)

# Restricting to malignant cells and the malignant state classification.
cell_class_filt <- mut_md_verhaak_annot[mut_md_verhaak_annot$CellID %in%CARE_filt_rna_tme$CellID, ]
cell_class_filt_ord <- cell_class_filt[match(CARE_filt_rna_tme$CellID, cell_class_filt$CellID), ]
all(CARE_filt_rna_tme$CellID==cell_class_filt_ord$CellID)
CARE_filt_rna_tme$CellStateGroup <- cell_class_filt_ord$CellTypeArchR


CARE_filt_rna_tme <- saveArchRProject(ArchRProj = CARE_filt_rna_tme, outputDirectory = "Save-CAREmut-TME-RNA-Longitudinal", load = TRUE, dropCells = TRUE) 

CARE_filt_rna_tme <- loadArchRProject("Save-CAREmut-TME-RNA-Longitudinal")

all(CARE_filt_rna_tme$CellID==cell_class_filt_ord$CellID)
CARE_filt_rna_tme$care_id <- cell_class_filt_ord$care_id



# CRITICAL STEP!
# See following github discussion for treatment vs. control style anylases: https://github.com/GreenleafLab/ArchR/discussions/696
CARE_filt_rna_tme <- addCellColData(ArchRProj = CARE_filt_rna_tme, data = paste0(CARE_filt_rna_tme@cellColData$CellStateGroup,"_x_",CARE_filt_rna_tme@cellColData$Sample), name = "StateBySample", cells = getCellNames(CARE_filt_rna_tme), force = TRUE)
table(CARE_filt_rna_tme$StateBySample)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Make pseudo bulk measurements for cell state groups by sample
# Pseudo-bulk refers to a grouping of single cells where the data from each single cell is combined into a single pseudo-sample.
CARE_filt_rna_tme <- addGroupCoverages(ArchRProj = CARE_filt_rna_tme, 
                                             groupBy = "StateBySample",  #  Oligodendrocyte_x_NL04-1
                                             minCells = 40,  # default
                                             maxCells = 500,  # default
                                             threads = getArchRThreads(),
                                             # Overwrite the data in the ArchRProject object if the pseudo-bulk replicate information already exists
                                             force = TRUE)

## Is macs2 in the path variable?
pathToMacs2 <- findMacs2()

# Iterative overlap peak merging procedure
CARE_filt_rna_tme <- addReproduciblePeakSet(
  ArchRProj = CARE_filt_rna_tme, 
  groupBy = "StateBySample", 
  pathToMacs2 = pathToMacs2,
  threads = getArchRThreads(),
)

## Needed to derive marker peaks.
CARE_filt_rna_tme <- addPeakMatrix(CARE_filt_rna_tme)

# Save this following peak matrix creation so that it is easier to run in the future.
CARE_filt_rna_tme <- saveArchRProject(ArchRProj = CARE_filt_rna_tme, outputDirectory = "Save-CAREmut-TME-RNA-Longitudinal", load = TRUE)

CARE_filt_rna_tme <- loadArchRProject("Save-CAREmut-TME-RNA-Longitudinal")

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
## Add add motif annotations.
if("Motif" %ni% names(CARE_filt_rna_tme@peakAnnotation)){
  CARE_filt_rna_tme <- addMotifAnnotations(ArchRProj = CARE_filt_rna_tme, motifSet = "cisbp", name = "Motif")
}

## Add ENCODE annotations 
CARE_filt_rna_tme <- addArchRAnnotations(ArchRProj = CARE_filt_rna_tme, collection = "EncodeTFBS")

# Check to see what is available.
getAvailableMatrices(CARE_filt_rna_tme)

### ### ### ### ### ### ### ### ### ### ### 
#Patient specific longitudinal changes
### ### ### ### ### ### ### ### ### ### ### 

## Objective: For each patient and cell state loop through longitudinal comparisons so
# that differential peaks can be determined within patient while controlling for cell state.
table(CARE_filt_rna_tme$StateBySample)

# Set up for loop run
archr_project_name <- "CARE_filt_rna_tme"
patient_id <- unique(sapply(strsplit(CARE_filt_rna_tme@cellColData$Sample, "-"), "[[", 1))
cell_state = unique(CARE_filt_rna_tme@cellColData$CellStateGroup)

# Create empty vectors and a list to store the results.
up_peaks <- down_peaks <- use_cells <- bgd_cells <- top_up_tfs <- top_down_tfs <- c()
outlist <- list()

atac_md_filt <- as.data.frame(CARE_filt_rna_tme@cellColData)
table(atac_md_filt$CellStateGroup, atac_md_filt$Sample)

# Check to see whether there will be any cases where there are 3 tumor samples per patient. Manually remove any 3 time point samples with low numbers 
table(CARE_filt_rna_tme@cellColData$Sample)

# Determine the longitudinal difference in cell states. Since we are interested in more subtle changes, let's use a lower threshold for significant difference.
# Log2FC of 0.5 would represent approximately ~50% increase in accessibility.

for (i in 1:length(patient_id)){
  for(j in 1:length(cell_state)) {
    
    # Get the sample identifier.
    samples <- unique(CARE_filt_rna_tme@cellColData$Sample)[grep(patient_id[i], unique(CARE_filt_rna_tme@cellColData$Sample))]
    
    # Remove samples with low cell numbers where the patient has 3 samples. Select max interval for all others.
    samples <- samples[!samples%in%c("NL03-1", "SJ03-2", "SJ08-3", "SJ20-3")]
    initial = min(sapply(strsplit(samples, "-"), "[[", 2))
    recurrence = max(sapply(strsplit(samples, "-"), "[[", 2))
    
    # Setting the "use" group will represent results as positive fold change values when peaks are higher in that group.
    state_use = paste0(cell_state[j], "_x_", patient_id[i], "-", recurrence)
    # A negative fold change values when peaks are higher in the background group.
    state_bgd = paste0(cell_state[j], "_x_", patient_id[i], "-", initial)
    
    out <- paste0(sprintf("Found use group: %s", state_use), " vs ", sprintf("Found background group: %s", state_bgd))
    print(out)
    
    use_cells = sum(CARE_filt_rna_tme@cellColData$StateBySample==sprintf("%s", state_use))
    bgd_cells = sum(CARE_filt_rna_tme@cellColData$StateBySample==sprintf("%s", state_bgd))
    
    # Skip to next iteration if number of cells is too low
    if(use_cells < 40) next
    if(bgd_cells < 40) next
    
    out <- paste0(sprintf("Analyzing patient %s for the comparison: ", patient_id[i]),  sprintf("Use group: %s", state_use), " vs ", sprintf("Background group: %s", state_bgd))
    print(out)
    
    start.time <- Sys.time()
    markerState <- getMarkerFeatures(
      ArchRProj = CARE_filt_rna_tme, 
      useMatrix = "PeakMatrix",
      groupBy = "StateBySample",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = sprintf("%s", state_use), 
      bgdGroups = sprintf("%s", state_bgd) 
    )
    
    motifsUp <- peakAnnoEnrichment(
      seMarker = markerState,
      ArchRProj = CARE_filt_rna_tme,
      peakAnnotation = "Motif",
      # Selecting a loosened threshold for these pairwise tests since they are less well powered.
      cutOff = "FDR <= 0.05 & Log2FC >= 0.5"
    )
    
    df_up <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
    df_up <- df_up[order(df_up$mlog10Padj, decreasing = TRUE),]
    df_up$rank <- seq_len(nrow(df_up))
    
    motifsDo <- peakAnnoEnrichment(
      seMarker = markerState,
      ArchRProj = CARE_filt_rna_tme,
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
    plotPDF(pv, name = sprintf("%s_Marker-Peaks-Volcano", comparison), width = 5, height = 5, ArchRProj = CARE_filt_rna_tme, addDOC = FALSE)
    
    results <- data.frame(patient_id[i], cell_state[j], state_use, use_cells, state_bgd, bgd_cells, up_peaks, down_peaks, top_up_tfs, top_down_tfs)
    outlist <- append(outlist, list(results))
  }
}

# Combine all into one data.frame.
patient_list_results <- do.call(rbind, outlist)

sample <- read.table("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/caremut_sample_identifier_linker.txt", sep = "\t", row.names = 1, header = TRUE)
sample_filt <- sample %>% 
  dplyr::select(care_id:SampleID, idh_codel_subtype)
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt", sep="\t", header = TRUE)

patient_list_results_filt <- patient_list_results %>% 
  filter(bgd_cells>50, use_cells>50) %>% 
  mutate(SampleID = sapply(strsplit(state_use, "_x_"), "[[", 2)) %>% 
  inner_join(sample_filt, by="SampleID") %>% 
  filter(state_use!=state_bgd)

# Write out the results so that they can be used downstream. 
write.table(patient_list_results_filt, "/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/archr_care_tme_cell_type_longitudinal_peaks_min50cells_perstate_20251015.txt", sep="\t", row.names = FALSE, col.names = TRUE)


tme <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/archr_care_tme_cell_type_longitudinal_peaks_min50cells_perstate_20251015.txt", header = TRUE, sep = "\t")
malignant <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/archr_care_cell_type_longitudinal_peaks_min100cells_perstate_20240506.txt", header = TRUE, sep = "\t")

# Add source labels
malignant <- malignant %>% mutate(compartment = "Malignant")
tme <- tme %>% mutate(compartment = "TME")

# Combine
df_all <- bind_rows(malignant, tme)

df_long <- summary_df %>%
  dplyr::select(cell_state.j., compartment, total_up, total_down) %>%
  tidyr::pivot_longer(cols = starts_with("total_"), names_to = "direction", values_to = "n_peaks") %>%
  dplyr::mutate(direction = recode(direction, total_up = "Up", total_down = "Down"))

ggplot(df_long, aes(x = reorder(cell_state.j., n_peaks), y = n_peaks, fill = direction)) +
  geom_col(position = "stack") +
  facet_wrap(~compartment, scales = "fixed") +
  coord_flip() +
  scale_fill_manual(values = c("Up" = "#e31a1c", "Down" = "#1f78b4")) +
  labs(
    x = "Cell State",
    y = "Number of Differential Peaks",
    title = "Up- and Down-Regulated Peaks in Malignant vs TME Cell States"
  ) +
  theme_bw(base_size = 13)


# Summarize per patient
patient_summary <- df_all %>%
  dplyr::group_by(compartment, patient_id, cell_state.j.) %>%
  dplyr::summarize(
    up_peaks = sum(up_peaks, na.rm = TRUE),
    down_peaks = sum(down_peaks, na.rm = TRUE),
    total_peaks = up_peaks + down_peaks
  ) %>%
  ungroup()

summary_df <- patient_summary %>%
  dplyr::group_by(compartment, cell_state.j.) %>%
  dplyr::summarize(
    mean_total = mean(total_peaks, na.rm = TRUE),
    sd_total = sd(total_peaks, na.rm = TRUE),
    n_patients = n(),
    mean_up = mean(up_peaks, na.rm = TRUE),
    mean_down = mean(down_peaks, na.rm = TRUE)
  ) %>%
  ungroup()

levels()


patient_summary <- patient_summary %>% 
  mutate(cell_state = recode(cell_state.j., `AC`  = "AC-like",
                             `NPC` = "NPC-like",
                             `OPC` = "OPC-like",
                             `MES` = "MES-like",
                             `TAM_MG-like` = "TAM MG-like",
                             `TAM_Inflammatory` = "TAM Inflammatory",
                             `TAM_BMDM-like` = "TAM BMDM-like"))

patient_summary$cell_state <- factor(patient_summary$cell_state, levels=c("OPC-like", "AC-like", "Undifferentiated", "NPC-like", "MES-like",
                                                                          "Oligodendrocyte", "TAM MG-like", "TAM Inflammatory", "TAM BMDM-like", "Astrocyte", "Endothelial"))

pdf(paste0(fig_dir, "longitudinal_peaks_malignant_vs_tme_by_state_20251023.pdf"), width=6.5, height=4, useDingbats = FALSE, bg = "transparent")
ggplot(patient_summary, aes(x = cell_state, y = total_peaks, fill = compartment)) +
  geom_boxplot(
    outlier.shape = NA,      
    alpha = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  geom_jitter(
    aes(color = compartment),
    size = 2,
    alpha = 0.8,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  ) +
  scale_fill_manual(values = c("Malignant" = "#d73027", "TME" = "#4575b4")) +
  scale_color_manual(values = c("Malignant" = "#d73027", "TME" = "#4575b4")) +
  labs(
    x = "Cell State",
    y = "Differentially Accessible Peaks per Patient",
    fill = "Compartment"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  guides(color=FALSE)
dev.off()

# 
# Run Wilcoxon test
wilcox_res <- wilcox.test(total_peaks ~ compartment, data = patient_summary)

# Extract and format p-value
pval <- wilcox_res$p.value
pval
pval_label <- "P = 1.2e-8"

# Examine the differences across compartments for the number of differentially accessible peaks
table(patient_summary$compartment)

pdf(paste0(fig_dir, "longitudinal_peaks_malignant_vs_tme_by_compartment_20251023.pdf"), width=4, height=4, useDingbats = FALSE, bg = "transparent")
ggplot(patient_summary, aes(x = compartment, y = total_peaks, fill = compartment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(
    aes(color = compartment),
    width = 0.1,
    size = 2.5,
    alpha = 0.8
  ) +
  annotate(
    "text",
    x = 1.5,  # roughly center between boxes
    y = max(patient_summary$total_peaks, na.rm = TRUE) * 1.05,  # place slightly above max
    label = pval_label,
    size = 4.5,
    fontface = "italic"
  ) +
  scale_fill_manual(values = c("Malignant" = "#d73027", "TME" = "#4575b4")) +
  scale_color_manual(values = c("Malignant" = "#d73027", "TME" = "#4575b4")) +
  labs(
    x = "Compartment",
    y = "Differentially Accessible Peaks per Patient"
  ) +
  plot_theme +
  guides(fill=FALSE, color=FALSE)
 dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Patient-specific longitudinal gene changes
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

## Objective: For each patient and cell state loop through longitudinal comparisons so
# that differential gene accessibility can be determined within patient while controlling for cell state.

# Set up for loop run
archr_project_name <- "CARE_filt_rna_tme"
patient_id <- unique(sapply(strsplit(CARE_filt_rna_tme@cellColData$Sample, "-"), "[[", 1))
cell_state = unique(CARE_filt_rna_tme@cellColData$CellStateGroup)

# Create empty vectors and a list to store the results.
up_genes <- down_genes <- use_cells <- bgd_cells <- c()
outlist <- list()

atac_md_filt <- as.data.frame(CARE_filt_rna_tme@cellColData)
table(atac_md_filt$CellStateGroup, atac_md_filt$Sample)

# Check to see whether there will be any cases where there are 3 tumor samples per patient. Manually remove any 3 time point samples with low numbers 
table(CARE_filt_rna_tme@cellColData$Sample)

# Determine the longitudinal difference in cell states. Since we are interested in more subtle changes, let's use a lower threshold for significant difference.
# Log2FC of 1 would represent approximately ~100% increase in accessibility.

for (i in 1:length(patient_id)){
  for(j in 1:length(cell_state)) {
    
    # Get the sample identifier.
    samples <- unique(CARE_filt_rna_tme@cellColData$Sample)[grep(patient_id[i], unique(CARE_filt_rna_tme@cellColData$Sample))]
    
    # Remove samples with low cell numbers where the patient has 3 samples. Select max interval for all others.
    samples <- samples[!samples%in%c("NL03-1", "SJ03-2", "SJ08-3", "SJ20-3")]
    initial = min(sapply(strsplit(samples, "-"), "[[", 2))
    recurrence = max(sapply(strsplit(samples, "-"), "[[", 2))
    
    # Setting the "use" group will represent results as positive fold change values when peaks are higher in that group.
    state_use = paste0(cell_state[j], "_x_", patient_id[i], "-", recurrence)
    # A negative fold change values when peaks are higher in the background group.
    state_bgd = paste0(cell_state[j], "_x_", patient_id[i], "-", initial)
    
    out <- paste0(sprintf("Found use group: %s", state_use), " vs ", sprintf("Found background group: %s", state_bgd))
    print(out)
    
    use_cells = sum(CARE_filt_rna_tme@cellColData$StateBySample==sprintf("%s", state_use))
    bgd_cells = sum(CARE_filt_rna_tme@cellColData$StateBySample==sprintf("%s", state_bgd))
    
    # Skip to next iteration if number of cells is too low
    if(use_cells < 50) next
    if(bgd_cells < 50) next
    
    out <- paste0(sprintf("Analyzing patient %s for the comparison: ", patient_id[i]),  sprintf("Use group: %s", state_use), " vs ", sprintf("Background group: %s", state_bgd))
    print(out)
    
    start.time <- Sys.time()
    # Note that the maximum cells allowed to be used for a single group are 500 by default. This is okay for these analyses because there are some samples 
    # with many more cells at a given time point.
    markerState <- getMarkerFeatures(
      ArchRProj = CARE_filt_rna_tme, 
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
    write.table(dag_df, paste0("/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/state_controlled_dag/tme/", comparison, "_longitudinal_dag.txt"), sep="\t", row.names = FALSE, col.names = TRUE)
    
    results <- data.frame(patient_id[i], cell_state[j], state_use, use_cells, state_bgd, bgd_cells, up_genes, down_genes)
    outlist <- append(outlist, list(results))
  }
}

# Combine all into one data.frame.
patient_gene_list_results <- do.call(rbind, outlist)

patient_gene_list_results_filt <- patient_gene_list_results %>% 
  filter(bgd_cells>50, use_cells>50) %>% 
  mutate(SampleID = sapply(strsplit(state_use, "_x_"), "[[", 2)) %>% 
  inner_join(sample_filt, by="SampleID") %>% 
  filter(state_use!=state_bgd)

# Write out the results so that they can be used downstream. 
write.table(patient_gene_list_results_filt, "/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/archr_care_tme_cell_type_longitudinal_genes_min50cells_perstate_20251015.txt", sep="\t", row.names = FALSE, col.names = TRUE)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Tumor type-specific differences in chromatin accessibility gene scores changes
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

final_atac_df <- data.frame(CARE_filt_rna_tme$cellNames, CARE_filt_rna_tme$CellID, CARE_filt_rna_tme@cellColData$care_id)
colnames(final_atac_df) <- c("cellNames", "CellID", "care_id")
final_atac_df <- final_atac_df %>% 
  left_join(sample_filt, by="care_id")
all(final_atac_df$CellID==CARE_filt_rna_tme$CellID)

CARE_filt_rna_tme$idh_codel_subtype <- final_atac_df$idh_codel_subtype

CARE_filt_rna_tme <- addCellColData(ArchRProj = CARE_filt_rna_tme, data = paste0(CARE_filt_rna_tme@cellColData$CellStateGroup,"_x_",CARE_filt_rna_tme@cellColData$idh_codel_subtype), name = "StateBySubtype", cells = getCellNames(CARE_filt_rna_tme), force = TRUE)
table(CARE_filt_rna_tme$StateBySubtype)


# Compare across subtypes. Note that the maximum cells allowed to be used for a single group are 500 by default, which should be acceptable for most analyses in the TME.

# Cell types to compare (minimum 1,000 cells each): Astrocyte, Oligodendrocyte, ExcNeuron, TAM-MG

markerState <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_tme, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "StateBySubtype",
  testMethod = "wilcoxon",
  maxCells = 1000,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Astrocyte_x_IDHmut-noncodel", 
  bgdGroups = "Astrocyte_x_IDHmut-codel"
)

upregulated <- getMarkers(markerState, cutOff = "FDR <= 0.05 & Log2FC >= 1")
downregulated <- getMarkers(markerState, cutOff = "FDR <= 0.05 & Log2FC <= -1")
upregulated_genes <- as.data.frame(upregulated$`Astrocyte_x_IDHmut-noncodel`)
downregulated_genes <- as.data.frame(downregulated$`Astrocyte_x_IDHmut-noncodel`)
dag_genes_astrocyte <- bind_rows(upregulated_genes, downregulated_genes) %>% 
  mutate(state = "Astrocyte")

markerState_oligo <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_tme, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "StateBySubtype",
  testMethod = "wilcoxon",
  maxCells = 1000,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Oligodendrocyte_x_IDHmut-noncodel", 
  bgdGroups = "Oligodendrocyte_x_IDHmut-codel"
)

upregulated_oligo <- getMarkers(markerState_oligo, cutOff = "FDR <= 0.05 & Log2FC >= 1")
downregulated_oligo <- getMarkers(markerState_oligo, cutOff = "FDR <= 0.05 & Log2FC <= -1")
upregulated_genes_oligo <- as.data.frame(upregulated_oligo$`Oligodendrocyte_x_IDHmut-noncodel`)
downregulated_genes_oligo <- as.data.frame(downregulated_oligo$`Oligodendrocyte_x_IDHmut-noncodel`)
dag_genes_oligodendrocyte <- bind_rows(upregulated_genes_oligo, downregulated_genes_oligo) %>% 
  mutate(state = "Oligodendrocyte")

markerState_neuron <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_tme, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "StateBySubtype",
  testMethod = "wilcoxon",
  maxCells = 1000,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "ExcNeuron_x_IDHmut-noncodel", 
  bgdGroups = "ExcNeuron_x_IDHmut-codel"
)

upregulated_neuron <- getMarkers(markerState_neuron, cutOff = "FDR <= 0.05 & Log2FC >= 1")
downregulated_neuron <- getMarkers(markerState_neuron, cutOff = "FDR <= 0.05 & Log2FC <= -1")
upregulated_genes_neuron <- as.data.frame(upregulated_neuron$`ExcNeuron_x_IDHmut-noncodel`)
downregulated_genes_neuron <- as.data.frame(downregulated_neuron$`ExcNeuron_x_IDHmut-noncodel`)
dag_genes_neuron <- bind_rows(upregulated_genes_neuron, downregulated_genes_neuron) %>% 
  mutate(state = "ExcNeuron")



markerState_microglia <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_tme, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "StateBySubtype",
  testMethod = "wilcoxon",
  maxCells = 1000,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "TAM_MG-like_x_IDHmut-noncodel", 
  bgdGroups = "TAM_MG-like_x_IDHmut-codel"
)

upregulated_microglia <- getMarkers(markerState_microglia, cutOff = "FDR <= 0.05 & Log2FC >= 1")
downregulated_microglia <- getMarkers(markerState_microglia, cutOff = "FDR <= 0.05 & Log2FC <= -1")
upregulated_genes_microglia <- as.data.frame(upregulated_microglia$`TAM_MG-like_x_IDHmut-noncodel`)
downregulated_genes_microglia <- as.data.frame(downregulated_microglia$`TAM_MG-like_x_IDHmut-noncodel`)
dag_genes_microglia <- bind_rows(upregulated_genes_microglia, downregulated_genes_microglia) %>% 
  mutate(state = "TAM_MG-like")


state_dag_genes <- bind_rows(dag_genes_astrocyte,
                             dag_genes_oligodendrocyte,
                             dag_genes_neuron,
                             dag_genes_microglia) %>% 
  mutate(direction = ifelse(Log2FC>=1, "IDH-A upregulated", "IDH-A downregulated"))

write.table(state_dag_genes, "/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/state_controlled_dag/tme/tme_state_controlled_subtype_dag.txt", sep="\t", row.names = FALSE, col.names = TRUE)

state_dag_genes <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/state_controlled_dag/tme/tme_state_controlled_subtype_dag.txt", sep = "\t")

top_genes <- state_dag_genes %>% 
  dplyr::group_by(name, seqnames) %>% 
  dplyr::summarise(gene_hits = n()) %>% 
  filter(gene_hits > 1)

malignant_dags <- read.delim( "/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/state_controlled_dag/state_controlled_subtype_dag.txt")

table(malignant_dags$seqnames, malignant_dags$direction)


top_genes <- malignant_dags %>% 
  dplyr::group_by(name, seqnames) %>% 
  dplyr::summarise(gene_hits = n()) %>% 
  filter(gene_hits > 3)

table(top_genes$seqnames)


cyto <- read.table("/vast/palmer/pi/verhaak/kcj28/reference/ucsc/hg38/cytoBand.txt.gz", sep = "\t", header = FALSE,
                   col.names = c("chrom", "start", "end", "band", "gieStain"))
cyto$start <- as.numeric(cyto$start)
cyto$end <- as.numeric(cyto$end)

library(dplyr)
arms <- cyto %>%
  mutate(arm = sub("([pq]).*", "\\1", band),
         chr_arm = paste0(chrom, arm)) %>% 
  filter(arm%in%c("p", "q")) %>% 
  dplyr::group_by(chrom, chr_arm) %>%  
  dplyr::summarise(
    start = min(start),
    end = max(end)) %>% 
  ungroup()
library(GenomicRanges)

df_fixed <- malignant_dags %>%
  mutate(
    start = pmin(start, end),
    end = pmax(start, end)
  )

gr_df <- GRanges(
  seqnames = df_fixed$seqnames,
  ranges = IRanges(start = df_fixed$start, end = df_fixed$end),
  gene = df_fixed$name, 
  state = df_fixed$state
)

gr_arms <- GRanges(seqnames = arms$chrom, 
                   ranges = IRanges(start = arms$start, end = arms$end),
                   arm = arms$chr_arm)

# find overlaps
hits <- findOverlaps(gr_df, gr_arms)
df_fixed$arm <- NA
df_fixed$arm[queryHits(hits)] <- gr_arms$arm[subjectHits(hits)]

table(df_fixed$arm)

dr_prop <- df_fixed %>%
  # Filter to only 1p or 19q
  dplyr::mutate(on_1p19q = arm %in% c("chr1p", "chr19q")) %>%
  dplyr::group_by(state) %>%
  dplyr::summarise(
    total_genes = n(),
    genes_on_1p19q = sum(on_1p19q, na.rm = TRUE),
    proportion_1p19q = genes_on_1p19q / total_genes
  ) %>%
  dplyr::arrange(desc(proportion_1p19q))

mean(dr_prop$proportion_1p19q)

ggplot(dr_prop, aes(x = state, y = proportion_1p19q)) +
  geom_col() +
  labs(
    x = NULL,
    y = "Proportion of genes on 1p or 19q",
    fill = "Chromosome arm"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )


### TME 


state_fixed <- state_dag_genes %>%
  mutate(
    start = pmin(start, end),
    end = pmax(start, end)
  )

state_gr_df <- GRanges(
  seqnames = state_fixed$seqnames,
  ranges = IRanges(start = state_fixed$start, end = state_fixed$end),
  gene = state_fixed$name, 
  state = state_fixed$state
)

# find overlaps
state_hits <- findOverlaps(state_gr_df, gr_arms)
state_fixed$arm <- NA
state_fixed$arm[queryHits(state_hits)] <- gr_arms$arm[subjectHits(state_hits)]

state_df_prop <- state_fixed %>%
  # Filter to only 1p or 19q
  dplyr::mutate(on_1p19q = arm %in% c("chr1p", "chr19q")) %>%
  dplyr::group_by(state) %>%
  dplyr::summarise(
    total_genes = n(),
    genes_on_1p19q = sum(on_1p19q, na.rm = TRUE),
    proportion_1p19q = genes_on_1p19q / total_genes
  ) %>%
  dplyr::arrange(desc(proportion_1p19q))

all_prop <- dr_prop %>% 
  bind_rows(state_df_prop)

ggplot(all_prop, aes(x = state, y = proportion_1p19q)) +
  geom_col() +
  labs(
    x = NULL,
    y = "Proportion of genes on 1p or 19q",
    fill = "Chromosome arm"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

table(state_fixed$state, state_fixed$direction)

df_fixed$compartment <- "Malignant"
state_fixed$compartment <- "TME"
all_combined <- state_fixed %>% 
  bind_rows(df_fixed)

df_summary <- all_combined %>%
  dplyr::mutate(group = ifelse(arm %in% c("chr1p", "chr19q"), "1p/19q", "Other")) %>%
  dplyr::group_by(state, group, compartment) %>%
  dplyr::summarise(
    n_genes = n(),
    .groups = "drop"
  )

df_summary$group <- factor(df_summary$group, levels=rev(c("1p/19q", "Other")))


pdf(paste0(fig_dir, "astrocytoma_vs_oligodendroglioma_dags.pdf"), width=7, height=5, useDingbats = FALSE, bg = "transparent")
ggplot(df_summary, aes(x = state, y = n_genes, fill = group)) +
  geom_col() +
  scale_fill_manual(values = c("1p/19q" = "#1f78b4", "Other" = "#d9d9d9")) +
  labs(
    x = "Cell state ",
    y = "Differentially accessible genes\nastrocytoma vs. oligodendroglioma",
    fill = "Chr. arm\ngroup"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  facet_grid(.~compartment, scales = "free")
dev.off()

# Visualize the 1p/19q yes/no 
df_1p19q_summary <- all_combined %>%
  dplyr::mutate(group = ifelse(arm %in% c("chr1p", "chr19q"), "1p/19q", "Other")) %>%
  dplyr::group_by(compartment, group) %>%
  dplyr::summarise(
    n_genes = n()) %>% 
  dplyr::mutate(prop = n_genes / sum(n_genes)) %>% 
  ungroup()


### END ###
