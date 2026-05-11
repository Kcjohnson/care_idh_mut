##############################
### Run ArchR peak calling analysis on CAREmut multiome ATAC data for all malignant states - NPC, OPC, Undiff, MES, AC
### Author: Kevin Johnson
##############################

## Generate malignant UMAP and perform differential gene/peak accessibility on RNA-defined MALIGNANT states across IDH-mutant tumors

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
library(harmony)

#### Set-up #####
## Check available cores. ArchR recommends setting total number of cores 1/2 to 3/4 of available cores,
num_cores <- detectCores() # e.g., 36
n_threads <- num_cores/2
addArchRThreads(threads = n_threads) 
## Each R session requires that the genome is also specified and must match alignment.
addArchRGenome("hg38")

# Load the ArchR object for IDHmut final analysis set. Doublets have been removed, only cells that intersect with passed qc for RNA, and RNA annotated cell states.
# This project comes from: `02_run_archr_caremut_all_cell_types.R`. 
CARE_filt_rna_all <- loadArchRProject("Save-CAREmut-All-RNA")
getAvailableMatrices(CARE_filt_rna_all) # GeneScoreMatrix and TileMatrix; 118,180 nuclei still including the Unresolved cells

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# CAREmut snRNA data that was processed by Seurat. Here is the associated metadata.
mut_md <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20260320.txt", sep = "\t", header = TRUE)
# Verhaak lab samples are the only samples with ATAC (via Multiome) data
mut_md_verhaak <- mut_md %>% 
  filter(lab=="Verhaak lab") %>% 
  # Remove SJ02-3 from analysis since it doesn't have any malignant cells that were classified (too few malignant cells overall)
  filter(SampleID!="SJ02-3")

# Restrict to IDHmut malignant cells that are classified by an RNA-based state.
cell_class <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/scoring/caremut_malignant_cell_state_assignment.txt",  sep = "\t", header = TRUE)
cell_class <- cell_class %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "Oligo.",
                                    `IDHmut-noncodel` = "Astro."),
         State = recode(State, `MP_AC1_MUT` = "AC-like",
                        `MP_OPC_MUT` = "OPC-like",
                        `MP_NPC_MUT` = "NPC-like",
                        `MP_MES_MUT` = "MES-like",
                        `MP_AC2_MUT` = "AC-like",
                        "Undifferentiated" = "Undifferentiated")) 

cell_class_trim <- cell_class %>% 
  dplyr::select(CellID, State, isCC)
mut_md_verhaak_state <- mut_md_verhaak %>% 
  left_join(cell_class_trim, by="CellID") %>%  
  mutate(tmp = ifelse(CellType_final=="Malignant", State, CellType_final),
         group = gsub("-like", "", tmp)) %>% 
  dplyr::select(-tmp)

# There is a slight difference in ArchR and Seurat naming convention.
CARE_filt_rna_all$CellID <-  gsub("#", "-", CARE_filt_rna_all$cellNames)

# Restricting to malignant cells and the malignant state classification.
cell_class_filt <- mut_md_verhaak_state[mut_md_verhaak_state$CellID %in% CARE_filt_rna_all$CellID, ]
cell_class_filt_ord <- cell_class_filt[match(CARE_filt_rna_all$CellID, cell_class_filt$CellID), ]

# Sanity checks
all(CARE_filt_rna_all$CellID==cell_class_filt_ord$CellID)
CARE_filt_rna_all$CellID==cell_class_filt_ord$CellID

CARE_filt_rna_all$CellStateGroup <- cell_class_filt_ord$group
CARE_filt_rna_all$scna_burden <- cell_class_filt_ord$scna_burden

# Subset entire CARE IDHmut object to only MALIGNANT cells.
idxSample <- BiocGenerics::which(CARE_filt_rna_all$CellType_final%in%"Malignant")
cellsSample <- CARE_filt_rna_all$cellNames[idxSample]

# Create a subset of the ArchR object to save malignant-only analyses.
CARE_filt_rna_malignant <- subsetArchRProject(
  ArchRProj = CARE_filt_rna_all,
  cells = cellsSample,
  outputDirectory = "Save-CAREmut-Malignant-RNA",
  force = TRUE)

# Sanity check that it was correctly loaded - otherwise it may write results to the old directory. What are the median QC values for malignant cells? 9.139 TSS and 13,927 fragments.
CARE_filt_rna_malignant <- loadArchRProject("Save-CAREmut-Malignant-RNA") # 71,365 malignant nuclei

# Plot frags against TSS enrichment for malignant cells only. Similar to analysis of all cells considered for analysis.
df <- getCellColData(CARE_filt_rna_malignant, select = c("log10(nFrags)", "TSSEnrichment"))
df

p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

pdf(paste0(fig_dir, "archr_malignant_cells_nfrags_vs_tss.pdf"), width = 5, height = 5, useDingbats = FALSE)
p
dev.off()

# Add grade information from clinical tables
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt", sep="\t", header = TRUE)
sample_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240608.txt", sep="\t", header = TRUE)

atac_md_mut <- data.frame(getCellColData(CARE_filt_rna_malignant))
atac_md_mut_annot <- atac_md_mut %>% 
  left_join(sample_md, by=c("care_id", "sample_barcode", "patient_id", "timepoint", "idh_codel_subtype")) %>% 
  mutate(atacCellNames = rownames(atac_md_mut))

# Do these cell names retain the same order?
ifelse(all(getCellNames(CARE_filt_rna_malignant)==atac_md_mut_annot$atacCellNames),
       sprintf("All cell names match. Proceed"), sprintf("Warning! Cell names do not match!"))
getCellNames(CARE_filt_rna_malignant)==atac_md_mut_annot$atacCellNames

# Add a few RNA features to the ArchR object.
CARE_filt_rna_malignant$Grade <- paste0("G", atac_md_mut_annot$grade_num)
CARE_filt_rna_malignant$idh_codel_subtype <- atac_md_mut_annot$idh_codel_subtype
CARE_filt_rna_malignant$tumor_type <-  ifelse(CARE_filt_rna_malignant$idh_codel_subtype=="IDHmut-codel", "Oligo", "Astro")
CARE_filt_rna_malignant$subtype_grade <- paste0(CARE_filt_rna_malignant$tumor_type, "_", CARE_filt_rna_malignant$Grade)
CARE_filt_rna_malignant$hypermutation <- ifelse(atac_md_mut_annot$hypermutation==1, "HM", "Non-HM")
CARE_filt_rna_malignant$hypermutation <- ifelse(is.na(CARE_filt_rna_malignant$hypermutation), "Non-HM", CARE_filt_rna_malignant$hypermutation)

# What's the breakdown of sample information: more oligodendroglioma malignant cells due to quality and sample number (50K vs 21K)
table(CARE_filt_rna_malignant$idh_codel_subtype)
table(CARE_filt_rna_malignant$Sample, CARE_filt_rna_malignant$subtype_grade)
table(CARE_filt_rna_malignant$Sample, CARE_filt_rna_malignant$CellStateGroup)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
## Important to remember that inaccessible chromatin via ATAC can be "non-accessible" or "not sampled". 1s have information and 0s do not.
## These are largely default parameters for LSI and Harmony
set.seed(123)
CARE_filt_rna_malignant <- addIterativeLSI(
  ArchRProj = CARE_filt_rna_malignant,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( # See Seurat::FindClusters. Change parameters depending on analysis goal.
    resolution = c(0.2), 
    sampleCells = 20000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:20, 
  force = TRUE
)

## Clustering is performed with the same methods from scRNAseq relying on Seurat's functionality here.
CARE_filt_rna_malignant <- addClusters(
  input = CARE_filt_rna_malignant,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.2,
  force = TRUE
)

# This will be the "uncorrected" UMAP
CARE_filt_rna_malignant <- addUMAP(
  ArchRProj = CARE_filt_rna_malignant, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 20, 
  minDist = 0.3, 
  metric = "cosine", 
  force = TRUE
)

# Repeat with Harmony batch correction
set.seed(123)
CARE_filt_rna_malignant <- addHarmony(
  ArchRProj = CARE_filt_rna_malignant,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = c("Sample"),
  force=TRUE
)

CARE_filt_rna_malignant <- addClusters(
  input = CARE_filt_rna_malignant,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "HarmonyClusters",
  resolution = 0.2,
  force = TRUE
)

CARE_filt_rna_malignant <- addUMAP(
  ArchRProj = CARE_filt_rna_malignant, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 20, 
  minDist = 0.3, 
  metric = "cosine",
  force = TRUE
)


# Malignant cell state colors
cols <- c("#AA2756", "#F77D58", "#7fbf7b", "#E8F5A3", "gray90")
names(cols)  <- names(table(CARE_filt_rna_malignant$CellStateGroup)) 
# Tumor type colors 
idh_cols <- c("#800074", "#298C8C")
names(idh_cols) <- names(table(CARE_filt_rna_malignant$tumor_type)) 

# Use this to extract the legend from a plot
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#### LSI (Uncorrected) UMAP ####
### Tumor type ###
p_subtype <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "tumor_type", embedding = "UMAP", pal = idh_cols, labelMeans=FALSE) +
  labs(x="", y="", color="Tumor type", title="") + theme(panel.border=element_blank()) +  guides(color = guide_legend(override.aes = list(size = 4)))
legend_grob <- g_legend(p_subtype)
p_subtype <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "tumor_type", embedding = "UMAP", pal = idh_cols, labelMeans=FALSE) +
  labs(x="", y="", color="Tumor type", title="") + theme(panel.border=element_blank())  + guides(color=FALSE)

plotPDF(p_subtype, name = "Unadjusted_malignant_UMAP_subtype.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_subtype.png"), p_subtype, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_subtype_legend.pdf"), legend_grob, width = 3, height = 3, dpi = 300)

### Patient ###
p_sample <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "patient_id", embedding = "UMAP", labelMeans=FALSE) +
  labs(x="", y="", color="Sample", title="") + theme(panel.border=element_blank()) + guides(color = guide_legend(override.aes = list(size = 4)))
legend_grob <- g_legend(p_sample)
p_sample <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "patient_id", embedding = "UMAP", labelMeans=FALSE) +
  labs(x="", y="", color="Sample", title="") + theme(panel.border=element_blank()) + guides(color=FALSE)

plotPDF(p_sample, name = "unadjusted_malignant_umap_sample.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_sample.png"), p_sample, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_sample_legend.pdf"), legend_grob, width = 4, height = 3, dpi = 300)

### Cell state ###
p_cell_state <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "CellStateGroup", embedding = "UMAP", pal = cols, labelMeans=FALSE) +
  labs(x="", y="", color="RNA cell state", title="") + theme(panel.border=element_blank()) +  guides(color = guide_legend(override.aes = list(size = 4))) 
legend_grob <- g_legend(p_cell_state)
p_cell_state <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "CellStateGroup", embedding = "UMAP", pal = cols, labelMeans=FALSE) +
  labs(x="", y="", color="RNA cell state", title="") + theme(panel.border=element_blank())  + guides(color=FALSE)

plotPDF(p_cell_state, name = "unadjusted_malignant_umap_cell_state.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_cell_state.png"), p_cell_state, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_cell_state_legend.pdf"), legend_grob, width = 5, height = 3, dpi = 300)


### Hypermutation ####
hyper_cols <- c("red", "gray80")
names(hyper_cols) <- names(table(CARE_filt_rna_malignant$hypermutation)) 
p_hm <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "hypermutation", embedding = "UMAP", pal = hyper_cols, labelMeans=FALSE) +
  labs(x="", y="", color="Hypermutation status", title="") + theme(panel.border=element_blank()) +  guides(color = guide_legend(override.aes = list(size = 4)))
legend_grob <- g_legend(p_hm)
p_hm <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "hypermutation", embedding = "UMAP", pal = hyper_cols, labelMeans=FALSE) +
  labs(x="", y="", color="Hypermutation status", title="") + theme(panel.border=element_blank())  + guides(color=FALSE)


plotPDF(p_hm, name = "unadjusted_malignant_hypermutation.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_hypermutation.png"), p_hm, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_hypermutation_legend.pdf"), legend_grob, width = 3, height = 3, dpi = 300)

#### snRNA SCNA burden ####
p_scna <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "scna_burden", embedding = "UMAP")  +
  labs(x="", y="", color="SCNA burden", title="") + theme(panel.border=element_blank()) + guides(color = guide_legend(override.aes = list(size = 4)))
legend_grob <- g_legend(p_scna)
p_scna <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "scna_burden", embedding = "UMAP") +
  labs(x="", y="", color="SCNA burden", title="") + theme(panel.border=element_blank())  + guides(fill=FALSE)

plotPDF(p_scna, name = "unadjusted_malignant_scna.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_scna.png"), p_scna, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_scna_legend.pdf"), legend_grob, width = 3, height = 3, dpi = 300)


#### LSI post-Harmony batch correction UMAP ####
### Cell State ####
p_harmony_cell_state <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "CellStateGroup", embedding = "UMAPHarmony", pal = cols, labelMeans=FALSE) +
  labs(x="", y="", color="RNA cell state", title="") + theme(panel.border=element_blank()) +  guides(color = guide_legend(override.aes = list(size = 4)))
legend_grob <- g_legend(p_harmony_cell_state)
p_harmony_cell_state <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "CellStateGroup", embedding = "UMAPHarmony", pal = cols, labelMeans=FALSE) +
  labs(x="", y="", color="RNA cell state", title="") + theme(panel.border=element_blank()) +  guides(color=FALSE)

plotPDF(p_harmony_cell_state, name = "harmony_malignant_umap_cell_state.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "harmony_malignant_umap_cell_state.png"), p_harmony_cell_state, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "harmony_malignant_umap_cell_state_legend.pdf"), legend_grob, width = 5, height = 3, dpi = 300)


### Tumor type ####
p_harmony_type <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "tumor_type", embedding = "UMAPHarmony", pal = idh_cols, labelMeans=FALSE) +
  labs(x="", y="", color="Tumor type", title="") + theme(panel.border=element_blank()) +  guides(color = guide_legend(override.aes = list(size = 4)))
legend_grob <- g_legend(p_harmony_type)
p_harmony_type <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "tumor_type", embedding = "UMAPHarmony", pal = idh_cols, labelMeans=FALSE) +
  labs(x="", y="", color="Tumor type", title="") + theme(panel.border=element_blank()) +  guides(color=FALSE)

plotPDF(p_harmony_type, name = "harmony_malignant_umap_tumor_type.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "harmony_malignant_umap_tumor_type.png"), p_harmony_type, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "harmony_malignant_umap_tumor_type_legend.pdf"), legend_grob, width = 3, height = 3, dpi = 300)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Define gene accessibility score markers based on malignant cell state.
set.seed(1)
markersGS <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "CellStateGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 2000
)

# Extract the marker list. These cutoffs are the default thresholds for ArchR. I tried to use these consistently throughout these analyses
markerListGS <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
lapply(markerListGS, nrow) # AC = 478, MES = 151, NPC = 7, OPC = 126, Undifferentiated = 91

# Examine the distribution of differentially accessible genes
undiff_df <- data.frame(markerListGS$Undifferentiated)
# DLL3, OLIG1, OLIG2
opc_df <- data.frame(markerListGS$OPC)
# No real difference for NPC-like cells
markerListGS$NPC
ac_df <- data.frame(markerListGS$AC)
# CD44
mes_df <- data.frame(markerListGS$MES)

# Selected marker genes of interest.
markerGenes <- c("AQP4", "CD44", "VIM", "TNC", "ANXA2", "SPARCL1",
                 "DLL3", "OLIG1",
                 "PDGFRA",
                 "HOXD11", "MET")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", 
  limits = c(-1.5, 1.5),
  labelMarkers = markerGenes,
  transpose = FALSE,
)

plotPDF(heatmapGS, name = "GeneScores-Malignant-Marker-Heatmap", width = 5.5, height = 6, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

pdf(paste0(fig_dir, "caremut_malignant_state_marker_genes.pdf"), width = 5, height = 4, useDingbats = FALSE, bg = "transparent")
heatmapGS
dev.off()


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Make pseudo bulk measurements for malignant cell state groups
# See github discussion on how to structure groups: https://github.com/GreenleafLab/ArchR/discussions/696
# Pseudo-bulk refers to a grouping of single cells where the data from each single sample is combined into a single pseudo-sample.
set.seed(123)
CARE_filt_rna_malignant <- addGroupCoverages(ArchRProj = CARE_filt_rna_malignant, 
                                               groupBy = "CellStateGroup",  # OPC, MES, etc
                                               minCells = 40,  # default
                                               maxCells = 500,  # default
                                               threads = getArchRThreads(),
                                               # Overwrite the data in the ArchRProject object if the pseudo-bulk replicate information already exists
                                               force = TRUE)

## Is macs2 in the path variable?
pathToMacs2 <- findMacs2()

# Iterative overlap peak merging procedure
CARE_filt_rna_malignant <- addReproduciblePeakSet(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "CellStateGroup", 
  pathToMacs2 = pathToMacs2,
  threads = getArchRThreads(),
)

## Needed to derive marker peaks.
CARE_filt_rna_malignant <- addPeakMatrix(CARE_filt_rna_malignant)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
getAvailableMatrices(CARE_filt_rna_malignant) # GeneScoreMatrix, PeakMatrix, TileMatrix

## Identifying marker peaks - features that are unique to a specific cell grouping.
# Account for biases in data quality via TSSEnrichment and nFrags - these are the default settings on the tutorial and seemed advisable.
set.seed(123)
markersPeaks <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "PeakMatrix", 
  groupBy = "CellStateGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 2000
)

saveRDS(markersPeaks, "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/malignant_cell_state_markersPeaks.RDS")
# markersPeaks <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/malignant_cell_state_markersPeaks.RDS")

# Extract the marker peaks. Get access the GRanges object via `returnGR = TRUE`.
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
lengths(markerList) # AC = 27447, MES = 15011, NPC = 1892, OPC = 13268, Undifferentiated = 6087

# Inspect some of the results
table(markerList$Undifferentiated@seqnames)
table(markerList$OPC@seqnames)
table(markerList$NPC@seqnames)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = FALSE,
  limits = c(-1.5, 1.5),
  nLabel = 1
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-mut-malignant-state", width = 5.5, height = 6, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

pdf(paste0(fig_dir, "caremut_malignant_state_marker_peaks.pdf"), width = 6, height = 4, useDingbats = FALSE, bg = "transparent")
heatmapPeaks
dev.off()

# Plot a few key marker peak's browser tracks.
p_opc <- plotBrowserTrack(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "CellStateGroup", 
  geneSymbol = c("OLIG1"),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)["OPC"],
  upstream = 20000,
  downstream = 20000
)

grid::grid.draw(p_opc$OLIG1)
plotPDF(p_opc, name = "Plot-Tracks-With-OLIG1", width = 5, height = 5, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

p_mes <- plotBrowserTrack(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "CellStateGroup", 
  geneSymbol = c("CD44"),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)["MES"],
  upstream = 20000,
  downstream = 20000
)

grid::grid.draw(p_mes$CD44)
plotPDF(p_mes, name = "Plot-Tracks-With-MES", width = 5, height = 5, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

p_undiff <- plotBrowserTrack(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "CellStateGroup", 
  geneSymbol = c("PDGFRA"),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)["Undifferentiated"],
  upstream = 50000,
  downstream = 50000
)

grid::grid.draw(p_undiff$PDGFRA)
plotPDF(p_undiff, name = "Plot-Tracks-With-PDGFRA", width = 5, height = 5, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

# Save the current object so that it can be re-loaded.
CARE_filt_rna_malignant_peaks <- saveArchRProject(ArchRProj = CARE_filt_rna_malignant, 
                                            outputDirectory = "Save-CAREmut-Malignant-RNA-Peaks", 
                                            load = TRUE) 


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
## Add add motif annotations. These enrichments give us some insights into the potential functionality of peak accessibility differences.
if("Motif" %ni% names(CARE_filt_rna_malignant_peaks@peakAnnotation)){
  CARE_filt_rna_malignant_peaks <- addMotifAnnotations(ArchRProj = CARE_filt_rna_malignant_peaks, motifSet = "cisbp", name = "Motif")
}

# For motifs amongst the open chromatin peak regions.
enrichRegions <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CARE_filt_rna_malignant_peaks,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

df_motif <- data.frame(TF = rownames(enrichRegions), mlog10Padj = assay(enrichRegions)[,2])
df_motif <- df_motif[order(df_motif$mlog10Padj, decreasing = TRUE),]
df_motif$rank <- seq_len(nrow(df_motif))

heatmapRegions <- plotEnrichHeatmap(enrichRegions, 
                                    transpose = TRUE, 
                                    cutOff = 5,
                                    clusterCols= FALSE)

plotPDF(heatmapRegions, name = "Regions-Enriched-Marker-Peak-Motif-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_malignant_peaks, addDOC = FALSE)

heatmapATAC_df <- plotEnrichHeatmap(enrichRegions, 
                                    n = 1500, 
                                    transpose = TRUE,
                                    returnMatrix = TRUE)

df <- data.frame(TF = rownames(enrichRegions), mlog10Padj = assay(enrichRegions))


library(viridisLite)
colnames(heatmapATAC_df) <- sapply(strsplit(colnames(heatmapATAC_df), " "), "[[", 1)
colnames(heatmapATAC_df) <- sapply(strsplit(colnames(heatmapATAC_df), "_"), "[[", 1)
enrichment_of_interest <- c("TCF12", "ASCL1", "CREB5", "TAL1", "JUNB", "FOS", "NFIC","SOX9", "POU2F3")
feature_order <- c("OPC", "NPC", "Undifferentiated", "MES", "AC")
enrichment_order <- c("TCF12", "ASCL1", "TAL1", "POU2F3", "CREB5", "JUNB", "FOS", "NFIC","SOX9")

heatmapATAC_df_filtered <- t(heatmapATAC_df[, colnames(heatmapATAC_df)%in%enrichment_of_interest])
heatmapATAC_df_ordered <- heatmapATAC_df_filtered[enrichment_order,feature_order]


manual_hmap <- ComplexHeatmap::Heatmap(heatmapATAC_df_ordered,                    
                                       show_row_dend = FALSE,
                                       cluster_columns = FALSE,
                                       cluster_rows = FALSE,
                                       show_column_dend = FALSE,
                                       col=viridis(100),
                                       name = "Norm. Enrichment -log10(P-adj) [0-Max]",
                                       heatmap_legend_param = list(
                                         legend_direction = "horizontal",
                                         legend_width = unit(5, "cm")
                                       ))


pdf(paste0(fig_dir, "malignant_states_archr_differential_peak_enrichment_motif_heatmap.pdf"), width = 5, height = 5, useDingbats = FALSE, bg = "transparent")
draw(manual_hmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
## Add a set of background peaks, which are used in computing deviations.
set.seed(123)
CARE_filt_rna_malignant_peaks <- addBgdPeaks(CARE_filt_rna_malignant_peaks)

# Available peak annotation
names(CARE_filt_rna_malignant_peaks@peakAnnotation)

## Compute the per-cell deviations across all of our motif annotations.
## This function has an optional parameter called matrixName that allows us to define the name of deviations.
## The option below creates a deviation matrix in each of the Arrow files called "DevMatrix" or "MotifMatrix". Force indicates whether the matrix listed should be overwritten.
# "Identifying Background Peaks!" will run if not background peaks haven't been added already.
CARE_filt_rna_malignant_peaks <- addDeviationsMatrix(
  ArchRProj = CARE_filt_rna_malignant_peaks, 
  peakAnnotation = "Motif",
  matrixName = "MotifMatrix",
  threads = getArchRThreads(),
  force = TRUE
)

motif_dev_df <- getVarDeviations(CARE_filt_rna_malignant_peaks, name = "MotifMatrix", plot = TRUE)

## How can one extract a subset of motifs for downstream analyses? getFeatures()
motifs <- c("ASCL1", "TCF12", "JUNB", "FOS", "RFX2", "TAL1", "NFIC")

markerMotifs <- getFeatures(CARE_filt_rna_malignant_peaks, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- markerMotifs[grep("z:", markerMotifs)]


p <- plotGroups(ArchRProj = CARE_filt_rna_malignant_peaks, 
                groupBy = "CellStateGroup", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(CARE_filt_rna_malignant_peaks)
)

plotPDF(p, name = "Plot-State-Motifs-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = CARE_filt_rna_malignant_peaks, addDOC = FALSE)

## Trying to extract the relevant TFs from the MotifMatrix. Use Z-scores for downstreams analyses.
motif_df <- getMatrixFromProject(
  ArchRProj = CARE_filt_rna_malignant_peaks,
  useMatrix = "MotifMatrix",
  useSeqnames = "z",
  verbose = TRUE,
  binarize = FALSE
)

motif_df_zscore <- assay(motif_df)
motif_df_zscore_out <- as.data.frame(as.matrix(t(motif_df_zscore)))
rownames(motif_df_zscore_out) <-  gsub("#", "-", rownames(motif_df_zscore_out))
motif_df_zscore_out$CellID <- rownames(motif_df_zscore_out)

# Write out the results so they can be used in additional analyses.
write.table(motif_df_zscore_out, "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/atac/archr_care_malignant_state_tf_motif_activity_zscore.txt", sep="\t", row.names = TRUE, col.names = TRUE)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
### Plot general enrichment of our malignant state peaks among other collections of ATACseq data, including normal and glioma samples.
CARE_filt_rna_malignant_peaks <- addArchRAnnotations(ArchRProj = CARE_filt_rna_malignant_peaks, collection = "ATAC")

enrichATAC <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CARE_filt_rna_malignant_peaks,
  peakAnnotation = "ATAC",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 5, cutOff = 3.5, transpose = TRUE)

# Most of the results are confirmatory (astrocytes with AC-like, opcs/oligodendrocytes with OPC-like, MES-like with GBM and other tumor types)
plotPDF(heatmapATAC, name = "ATAC-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_malignant_peaks, addDOC = FALSE)

# NPC-like cells share a similar chromatin profile with OPCs, but had only one significant enrichment for "Brain_Excitatory_neurons"
df <- data.frame(TF = rownames(enrichATAC), mlog10Padj = assay(enrichATAC)[,3])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)

# Save these results for future inspection.
enrichATAC_out <- assay(enrichATAC)
enrichATAC_out$Bulk <- rownames(enrichATAC_out)
write.table(enrichATAC_out, "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/atac/archr_care_mut_malignant_bulk_atac_enrichment_amongst_markerpeaks.txt", sep="\t", row.names = FALSE, col.names = TRUE)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Pairwise differential expression. In the snRNA data, there was a clear anti-correlation for OPC-like and AC-like cells. Assessing that difference here.
markerOPC <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant_peaks, 
  useMatrix = "PeakMatrix",
  groupBy = "CellStateGroup",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "OPC",
  bgdGroups = "AC"
)

# Approximately equal up- and down-regulated. Total peaks: 281,999 (9349 up and 11,004 down)
volcano_opc_v_ac <- plotMarkers(seMarker = markerOPC, name = "OPC", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "Volcano")

png(paste0(fig_dir, "opc_vs_ac_differential_peaks_volcano.png"), width = 5, height = 5, res = 300, units = "in")
volcano_opc_v_ac
dev.off()

markerMES_v_AC <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant_peaks, 
  useMatrix = "PeakMatrix",
  groupBy = "CellStateGroup",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "MES",
  bgdGroups = "AC"
)

# Total peaks: 281,999 (1808 up and 1474 down). Differences between AC and MES are considerably smaller (as expected) than AC-like vs OPC-like.
# This is also evident in the various UMAP plots.
volcano_mes_v_ac <- plotMarkers(seMarker = markerMES_v_AC, name = "MES", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "Volcano")
volcano_mes_v_ac

png(paste0(fig_dir, "mes_vs_ac_differential_peaks_volcano.png"), width = 5, height = 5, res = 300, units = "in")
volcano_mes_v_ac
dev.off()

MesMotifsUp <- peakAnnoEnrichment(
  seMarker = markerMES_v_AC,
  ArchRProj = CARE_filt_rna_malignant_peaks,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

# Domninated by JUN/FOS and STAT3 further down the list.
df <- data.frame(TF = rownames(MesMotifsUp), mlog10Padj = assay(MesMotifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Save various outputs:
getAvailableMatrices(CARE_filt_rna_malignant_peaks) #  "GeneScoreMatrix" "MotifMatrix" "PeakMatrix" "TileMatrix"  
names(CARE_filt_rna_malignant_peaks@peakAnnotation) # "Motif" "ATAC" 

care_malignant_archr_df <- data.frame(CARE_filt_rna_malignant_peaks@cellColData)
write.table(care_malignant_archr_df, "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/archr_care_mut_malignant_state_atac_metadata.txt", sep="\t", row.names = TRUE, col.names = TRUE)

paste0("Memory Size = ", round(object.size(CARE_filt_rna_malignant) / 10^6, 3), " MB")
CARE_filt_rna_malignant_peaks <- saveArchRProject(ArchRProj = CARE_filt_rna_malignant_peaks, 
                                                  outputDirectory = "Save-CAREmut-Malignant-RNA-Peaks", 
                                                  load = TRUE) 

### END ###