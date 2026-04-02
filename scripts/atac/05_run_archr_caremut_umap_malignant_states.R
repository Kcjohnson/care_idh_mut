##############################
### Run ArchR peak calling analysis on CAREmut multiome ATAC data for all malignant states
### Author: Kevin Johnson
### Updated: 2024.06.010
##############################

## Part 5: Peak calling on RNA-defined MALIGNANT states across IDH-mutant tumors

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
# devtools::install_github("immunogenomics/harmony")
library(harmony)

## Specify output directory to drop figures:
fig_dir <- "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/figures/archr/"

#### Set-up #####
## Check available cores. ArchR recommends setting total number of cores 1/2 to 3/4 of available cores,
num_cores <- detectCores() # e.g., 36
n_threads <- num_cores/2
addArchRThreads(threads = n_threads) 
## Each R session requires that the genome is also specified and must match alignment.
addArchRGenome("hg38")

# Load the ArchR object for IDHmut final analysis set. Doublets have been removed, only cells that intersect with passed qc for RNA, and RNA annotated cell states.
projMONITOR <- loadArchRProject("Save-AllSamples-2024")

## This larger dataset contains both IDHmut and IDHwt tumors that we've processed. Need to restrict to IDHmut-only tumors.
projMONITOR_filt <- filterDoublets(projMONITOR)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# CAREmut data that was processed by Seurat and cell types were assigned based on RNA.
mut_md <- read.table("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20240212.txt", sep = "\t", header = TRUE)
# For this project, only ATAC profiles were generated for samples processed by the Verhaak lab.
mut_md_verhaak <- mut_md %>% 
  filter(lab=="Verhaak lab", CellType_final=="Malignant")

## Create a data.frame that contains the essential information to be merged and inspected.
atac_df <- data.frame(projMONITOR_filt$cellNames, projMONITOR_filt@cellColData$Sample, projMONITOR_filt@cellColData$TSSEnrichment)

# The cell names are a little different between RNA (Seurat) and ATAC (ArchR). Need to create a common linker.
atac_df$CellID <- gsub("#", "-", atac_df$projMONITOR_filt.cellNames)
# 71088 cells present in ATACseq data that pass doublet removal and are also found in snRNAseq data that passed QC for IDH-mutant.
sum(atac_df$CellID%in%mut_md_verhaak$CellID)

# Combine the two data.frames by adding on the RNA data and filtering out what's leftover.
atac_df_filt_rna <- atac_df %>% 
  inner_join(mut_md_verhaak, by="CellID") 

## Identify which cells to keep in the analysis. 
tmp <- getCellNames(projMONITOR_filt)
rna_cells = tmp[which(tmp%in%atac_df_filt_rna$projMONITOR_filt.cellNames)]  

## Subset the cells to only those that also have RNA.
CARE_filt_rna_malignant <- subsetCells(ArchRProj = projMONITOR_filt, cellNames = rna_cells)

# Restrict to IDHmut malignant cells and classify by RNA-based state.
cell_class <- read.table("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/processed_data/rna/classification/caremut_all_select_state_assignment_20240416.txt",  sep = "\t", header = TRUE)
cell_class <- cell_class %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "IDH-O",
                                    `IDHmut-noncodel` = "IDH-A"),
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
  mutate(group = gsub("-like", "", ifelse(CellType_final=="Malignant", State, CellType_final))) 


CARE_filt_rna_malignant$CellID <- gsub("#", "-", CARE_filt_rna_malignant$cellNames)

# Restricting to malignant cells and the malignant state classification.
cell_class_filt <- mut_md_verhaak_state[mut_md_verhaak_state$CellID %in%CARE_filt_rna_malignant$CellID, ]
cell_class_filt_ord <- cell_class_filt[match(CARE_filt_rna_malignant$CellID, cell_class_filt$CellID), ]
all(CARE_filt_rna_malignant$CellID==cell_class_filt_ord$CellID)
CARE_filt_rna_malignant$CellStateGroup <- cell_class_filt_ord$group


#CARE_filt_rna_malignant <- saveArchRProject(ArchRProj = CARE_filt_rna_malignant, outputDirectory = "Save-CAREmut-Malignant-RNA-UMAP", load = TRUE, dropCells = TRUE) 

CARE_filt_rna_malignant <- loadArchRProject("Save-CAREmut-Malignant-RNA-UMAP")

all(CARE_filt_rna_malignant$CellID==cell_class_filt_ord$CellID)
CARE_filt_rna_malignant$care_id <- cell_class_filt_ord$care_id

# Add grade information from clinical tables
patient_md <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt", sep="\t", header = TRUE)
sample_md <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240608.txt", sep="\t", header = TRUE)

atac_md_mut <- data.frame(getCellColData(CARE_filt_rna_malignant))
atac_md_mut_annot <- atac_md_mut %>% 
  left_join(sample_md, by=c("care_id")) %>% 
  mutate(atacCellNames = rownames(atac_md_mut))

## Do these cell names retain the same order?
ifelse(all(getCellNames(CARE_filt_rna_malignant)==atac_md_mut_annot$atacCellNames),
       sprintf("All cell names match. Proceed"), sprintf("Warning! Cell names do not match!"))

## Add a few RNA features to the ArchR object.
CARE_filt_rna_malignant$Grade <- paste0("G", atac_md_mut_annot$grade_num)
CARE_filt_rna_malignant$idh_codel_subtype <- atac_md_mut_annot$idh_codel_subtype.x
CARE_filt_rna_malignant$subtype_grade <- paste0(CARE_filt_rna_malignant$idh_codel_subtype, "_", CARE_filt_rna_malignant$Grade)
CARE_filt_rna_malignant$patient_id <- atac_md_mut_annot$patient_id
CARE_filt_rna_malignant$hypermutation <- ifelse(atac_md_mut_annot$hypermutation==1, "HM", "Non-HM")
CARE_filt_rna_malignant$hypermutation <- ifelse(is.na(CARE_filt_rna_malignant$hypermutation), "Non-HM", CARE_filt_rna_malignant$hypermutation)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
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


CARE_filt_rna_malignant <- saveArchRProject(ArchRProj = CARE_filt_rna_malignant, outputDirectory = "Save-CAREmut-Malignant-RNA-UMAP", load = TRUE, dropCells = TRUE) 

CARE_filt_rna_malignant <- loadArchRProject("Save-CAREmut-Malignant-RNA-UMAP")


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
getAvailableMatrices(CARE_filt_rna_malignant) # GeneScoreMatrix, PeakMatrix, TileMatrix

## Identifying marker peaks - features that are unique to a specific cell grouping.
# Account for biases in data quality via TSSEnrichment and nFrags - these are the default settings on the tutorial and seem advisable.
markersPeaks <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "PeakMatrix", 
  groupBy = "CellStateGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 2000
)

saveRDS(markersPeaks, "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/malignant_cell_state_markersPeaks_20240507.RDS")
markersPeaks <- readRDS("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/malignant_cell_state_markersPeaks_20240507.RDS")

# Extract the marker peaks. Get access the GRanges object via `returnGR = TRUE`.
# Throughout these analyses, I am sticking with the FDR <= 0.05 and Log2FC >= 1 so that it is simple.
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")
lapply(markerList, nrow) # AC = 33454, MES = 14250, NPC = 2309, OPC = 12355, Undifferentiated = 5829

table(markerList$Undifferentiated$seqnames)
table(markerList$OPC$seqnames)
table(markerList$NPC$seqnames)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  limits = c(-1.5, 1.5),
  transpose = FALSE,
  nLabel = 1
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-mut-malignant-state_v2", width = 6, height = 4, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

pdf(paste0(fig_dir, "caremut_malignant_state_marker_peaks_20241004.pdf"), width = 6, height = 4, useDingbats = FALSE, bg = "transparent")
heatmapPeaks
dev.off()



## Define markers based on cell state.
markersGS <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "CellStateGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 2000
)

## Extract the marker list. These are the default thresholds for ArchR.
markerListGS <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
lapply(markerListGS, nrow) # AC = 449, MES = 181, NPC = 5, OPC = 120, Undifferentiated. = 108

# Examine the distribution of differentially accessible genes
undiff_df <- data.frame(markerListGS$Undifferentiated)
opc_df <- data.frame(markerListGS$OPC)
markerListGS$NPC
ac_df <- data.frame(markerListGS$AC)
mes_df <- data.frame(markerListGS$MES)


markerGenes <- c("AQP4", "CD44", "VIM", "TNC", "SPARCL1", "VMP1", "ANXA2",
                 "DLL3", "HOXD11",
                 "SLC4A4", "OLIG1")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", 
  limits = c(-1.5, 1.5),
  labelMarkers = markerGenes,
  transpose = FALSE,
)

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

pdf(paste0(fig_dir, "caremut_malignant_state_marker_genes.pdf"), width = 5, height = 4, useDingbats = FALSE, bg = "transparent")
heatmapGS
dev.off()


p <- plotBrowserTrack(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "CellStateGroup", 
  geneSymbol = c("OLIG1"),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)["OPC"],
  upstream = 20000,
  downstream = 20000
)

grid::grid.draw(p$OLIG1)
plotPDF(p, name = "Plot-Tracks-With-OLIG1", width = 5, height = 5, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# ATAC AddModuleScores for RNA-based metaprograms
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
mut_mp <- read.table("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/malignant_selected/care_mut_selected_malignant_metaprograms.csv", sep = ",", header = T, row.names = 1)
mut_mp_list <- lapply(names(mut_mp), function(col_name) mut_mp[[col_name]])
names(mut_mp_list) <-paste0(colnames(mut_mp), "_MUT")

# Need to filter to those genes present in the dataset.
all_genes <- getFeatures(CARE_filt_rna_malignant)

filtered_mp_list <- lapply(mut_mp_list, function(genes) {
  genes[genes %in% all_genes]
})
filtered_mp_list[["MP_AC1_MUT"]][filtered_mp_list[["MP_AC1_MUT"]]%in%ac_df$name]

# Get the differentially accessible genes
undiff_atac <-markerListGS$Undifferentiated
undiff_atac_sig_hits <- undiff_atac %>% 
  as.data.frame() %>% 
  arrange(Log2FC)
undiff_genes <- list(rev(undiff_atac_sig_hits$name))
names(undiff_genes) <- "undifferentiated_atac"

opc_atac <-markerListGS$OPC
opc_atac_sig_hits <- opc_atac %>% 
  as.data.frame() %>% 
  arrange(Log2FC)
opc_genes <- list(rev(opc_atac_sig_hits$name))
names(opc_genes) <- "opc_atac"

ac_atac <-markerListGS$AC
ac_atac_sig_hits <- ac_atac %>% 
  as.data.frame() %>% 
  arrange(Log2FC)
ac_genes <- list(rev(ac_atac_sig_hits$name))
names(ac_genes) <- "ac_atac"

mes_atac <-markerListGS$MES
mes_atac_sig_hits <- mes_atac %>% 
  as.data.frame() %>% 
  arrange(Log2FC)
mes_genes <- list(rev(mes_atac_sig_hits$name))
names(mes_genes) <- "mes_atac"

atac_genes <- c(mes_genes, ac_genes, opc_genes, undiff_genes)


CARE_filt_rna_malignant <- addModuleScore(CARE_filt_rna_malignant, features = atac_genes, useMatrix = "GeneScoreMatrix")
CARE_filt_rna_malignant <- addModuleScore(CARE_filt_rna_malignant, features = filtered_mp_list, useMatrix = "GeneScoreMatrix")

df <- CARE_filt_rna_malignant@cellColData
df_out <- df %>% as.data.frame()
write.table(df_out, "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/atac_malignant_module_scores.txt", sep="\t", row.names = FALSE, col.names = TRUE)




############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
## Add add motif annotations.
if("Motif" %ni% names(CARE_filt_rna_malignant@peakAnnotation)){
  CARE_filt_rna_malignant <- addMotifAnnotations(ArchRProj = CARE_filt_rna_malignant, motifSet = "cisbp", name = "Motif")
}
#CARE_filt_rna_malignant <- addArchRAnnotations(ArchRProj = CARE_filt_rna_malignant, collection = "EncodeTFBS")

# For motifs amongst the open chromatin peak regions.
enrichRegions <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CARE_filt_rna_malignant,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

df_motif <- assay(enrichRegions)
df_motif$TF <- rownames(enrichRegions)
df_motif <- data.frame(TF = rownames(enrichRegions), mlog10Padj = assay(enrichRegions)[,2])
df_motif <- df_motif[order(df_motif$mlog10Padj, decreasing = TRUE),]
df_motif$rank <- seq_len(nrow(df_motif))


heatmapRegions <- plotEnrichHeatmap(enrichRegions, 
                                    transpose = TRUE, 
                                    cutOff = 20,
                                    n = 10,
                                    pal = paletteContinuous(set = "comet", n = 100),
                                    clusterCols= TRUE)

plotPDF(heatmapRegions, name = "Regions-Enriched-Marker-Peak-Motif-Heatmap", width = 6, height = 4, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)


#
df <- data.frame(TF = rownames(enrichRegions), mlog10Padj = assay(enrichRegions)[,5])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))


## Add a set of background peaks, which are used in computing deviations.
CARE_filt_rna_malignant <- addBgdPeaks(CARE_filt_rna_malignant)

# Available peak annotation
names(CARE_filt_rna_malignant@peakAnnotation)

## Compute the per-cell deviations across all of our motif annotations.
## This function has an optional parameter called matrixName that allows us to define the name of deviations.
## The option below creates a deviation matrix in each of the Arrow files called "DevMatrix" or MotifMatrix. Force indicates whether the matrix listed should be overwritten.
# "Identifying Background Peaks!" will run if not background peaks haven't been added already.
CARE_filt_rna_malignant <- addDeviationsMatrix(
  ArchRProj = CARE_filt_rna_malignant, 
  peakAnnotation = "Motif",
  matrixName = "MotifMatrix",
  threads = getArchRThreads(),
  force = TRUE
)
# Not sure why this occurs.
# In mclapply(..., mc.cores = threads, mc.preschedule = preschedule) :
# 10 parallel function calls did not deliver results

# There may have been issue when saving the project during the last analysis.
motif_dev_df <- getVarDeviations(CARE_filt_rna_malignant, name = "MotifMatrix")
#Error in h5read(ArrowFile, paste0(subGroup, "/Info/FeatureDF")) : 
#  Object 'MotifMatrix/Info/FeatureDF' does not exist in this HDF5 file.

CARE_filt_rna_malignant <- saveArchRProject(ArchRProj = CARE_filt_rna_malignant, outputDirectory = "Save-CAREmut-Malignant-RNA-UMAP", load = TRUE) 


CARE_filt_rna_malignant <- loadArchRProject("Save-CAREmut-Malignant-RNA-UMAP")
############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
## Important to remember that inaccessible chromatin via ATAC can be "non-accessible" or "not sampled". 1s have information and 0s do not.
## These are largely default parameters:
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


CARE_filt_rna_malignant <- addHarmony(
  ArchRProj = CARE_filt_rna_malignant,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = c("Sample"),
  force=TRUE
)

## Clustering is performed with the same methods from scRNAseq relying on Seurat's functionality here.
## Default parameters are shown. Change parameters depending on selection.
CARE_filt_rna_malignant <- addClusters(
  input = CARE_filt_rna_malignant,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.2
)

CARE_filt_rna_malignant <- addClusters(
  input = CARE_filt_rna_malignant,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.2,
  force = TRUE
)

## Single nuclei embeddings. Parameter selection will be analysis dependent.
CARE_filt_rna_malignant <- addUMAP(
  ArchRProj = CARE_filt_rna_malignant, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 20, 
  minDist = 0.3, 
  metric = "cosine"
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


CARE_filt_rna_malignant <- saveArchRProject(ArchRProj = CARE_filt_rna_malignant, outputDirectory = "Save-CAREmut-Malignant-RNA-UMAP", load = TRUE) 

CARE_filt_rna_malignant <- loadArchRProject("Save-CAREmut-Malignant-RNA-UMAP")

# Clusters defined post-Harmony batch correction
p0 <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")

# Examine whether there are major differences across clusters
ClustermarkersGS <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 2000
)
# Inspect results
ClusterMarkerList <- getMarkers(ClustermarkersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
c4 <- as.data.frame(ClusterMarkerList$C4)

# Add some additional metadata to be plotted
all(CARE_filt_rna_malignant$CellID==cell_class_filt_ord$CellID)
CARE_filt_rna_malignant$isCC <- cell_class_filt_ord$isCC
CARE_filt_rna_malignant$scna_burden <- cell_class_filt_ord$scna_burden
CARE_filt_rna_malignant$subtype <- ifelse(CARE_filt_rna_malignant$idh_codel_subtype=="IDHmut-codel", "IDH-O", "IDH-A")

# Use this to extract the legend from a plot
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# Plot the non-batch corrected UMAP across several important features to highlight the effect of genetics:
# Sample
p_sample <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "care_id", embedding = "UMAP") +
  labs(x="", y="", color="Sample", title="") + theme(panel.border=element_blank()) + guides(color=FALSE)

plotPDF(p_sample, name = "unadjusted_malignant_umap_sample.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_sample.png"), p_sample, width = 4, height = 4, dpi = 300)

# Cell state colors
cols <- c("#AA2756", "#F77D58", "#7fbf7b", "#E8F5A3", "gray90")
names(cols)  <- names(table(CARE_filt_rna_malignant$CellStateGroup)) 
# Subtype
idh_cols <- c("#67A9CF", "#EF8A62")
names(idh_cols) <- names(table(CARE_filt_rna_malignant$subtype)) 
p_subtype <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "subtype", embedding = "UMAP", pal = idh_cols) +
  labs(x="", y="", color="Subtype", title="") + theme(panel.border=element_blank()) +  guides(color = guide_legend(override.aes = list(size = 4)))
legend_grob <- g_legend(p_subtype)
p_subtype <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "subtype", embedding = "UMAP", pal = idh_cols) +
  labs(x="", y="", color="Subtype", title="") + theme(panel.border=element_blank())  + guides(color=FALSE)

plotPDF(p_subtype, name = "Unadjusted_malignant_UMAP_subtype.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_subtype.png"), p_subtype, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_subtype_legend.pdf"), legend_grob, width = 3, height = 3, dpi = 300)

# Cell state
p_cell_state <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "CellStateGroup", embedding = "UMAP", pal = cols) +
  labs(x="", y="", color="RNA cell state", title="") + theme(panel.border=element_blank()) +  guides(color = guide_legend(override.aes = list(size = 4))) 
legend_grob <- g_legend(p_cell_state)
p_cell_state <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "CellStateGroup", embedding = "UMAP", pal = cols) +
  labs(x="", y="", color="RNA cell state", title="") + theme(panel.border=element_blank())  + guides(color=FALSE)

plotPDF(p_cell_state, name = "unadjusted_malignant_umap_cell_state.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_cell_state.png"), p_cell_state, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_cell_state_legend.pdf"), legend_grob, width = 3, height = 3, dpi = 300)

# Removing the numbers from the UMAP
p_cell_state <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "CellStateGroup", embedding = "UMAP", pal = cols,
                              labelMeans=FALSE) +
  labs(x="", y="", color="RNA cell state", title="") + theme(panel.border=element_blank()) +  guides(color = FALSE) 
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_cell_state_20241004.png"), p_cell_state, width = 4, height = 4, dpi = 300)

# Hypermutation
hyper_cols <- c("red", "gray80")
names(hyper_cols) <- names(table(CARE_filt_rna_malignant$hypermutation)) 
p_hm <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "hypermutation", embedding = "UMAP", pal = hyper_cols) +
  labs(x="", y="", color="Hypermutation status", title="") + theme(panel.border=element_blank()) +  guides(color = guide_legend(override.aes = list(size = 4)))
legend_grob <- g_legend(p_hm)
p_hm <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "hypermutation", embedding = "UMAP", pal = hyper_cols) +
  labs(x="", y="", color="Hypermutation status", title="") + theme(panel.border=element_blank())  + guides(color=FALSE)


plotPDF(p_hm, name = "unadjusted_malignant_hypermutation.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_hypermutation.png"), p_hm, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_hypermutation_legend.pdf"), legend_grob, width = 3, height = 3, dpi = 300)

# SCNA burden
p_scna <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "scna_burden", embedding = "UMAP")  +
  labs(x="", y="", color="SCNA burden", title="") + theme(panel.border=element_blank()) + guides(color = guide_legend(override.aes = list(size = 4)))
legend_grob <- g_legend(p_scna)
p_scna <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "scna_burden", embedding = "UMAP") +
  labs(x="", y="", color="SCNA burden", title="") + theme(panel.border=element_blank())  + guides(fill=FALSE)

plotPDF(p_scna, name = "unadjusted_malignant_scna.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_scna.png"), p_scna, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "unadjusted_malignant_umap_scna_legend.pdf"), legend_grob, width = 3, height = 3, dpi = 300)


### Harmony corrected UMAP ###
# Cell state
p_cell_state <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "CellStateGroup", embedding = "UMAPHarmony", pal = cols) +
  labs(x="", y="", color="RNA cell state", title="") + theme(panel.border=element_blank()) +  guides(color = guide_legend(override.aes = list(size = 4)))
legend_grob <- g_legend(p_cell_state)
p_cell_state <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "CellStateGroup", embedding = "UMAPHarmony", pal = cols) +
  labs(x="", y="", color="RNA cell state", title="") + theme(panel.border=element_blank())  + guides(color=FALSE)

plotPDF(p_cell_state, name = "harmony_malignant_umap_cell_state.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "harmony_malignant_umap_cell_state.png"), p_cell_state, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "harmony_malignant_umap_cell_state_legend.pdf"), legend_grob, width = 3, height = 3, dpi = 300)

p_cell_state <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "CellStateGroup", embedding = "UMAPHarmony", pal = cols,
                              labelMeans=FALSE) +
  labs(x="", y="", color="RNA cell state", title="") + theme(panel.border=element_blank()) +  guides(color = FALSE) 
ggsave(paste0(fig_dir, "harmony_malignant_umap_cell_state_20241004.png"), p_cell_state, width = 4, height = 4, dpi = 300)


### ###
# RNA metaprogram scores
### ###
scores <- read.table("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/processed_data/rna/classification/caremut_malignant_all_signature_scores_20240416.txt",  sep = "\t", header = TRUE)
scores_filt <- scores[scores$CellID %in%CARE_filt_rna_malignant$CellID, ]
scores_filt_ord <- scores_filt[match(CARE_filt_rna_malignant$CellID, scores_filt$CellID), ]
all(CARE_filt_rna_malignant$CellID==scores_filt_ord$CellID)

CARE_filt_rna_malignant$Stemness_RNA_Tirosh <- scores_filt_ord$Stemness_Tirosh2016
CARE_filt_rna_malignant$Stemness_RNA_Venteicher <- scores_filt_ord$Stemness_Venteicher2017
CARE_filt_rna_malignant$MP_OPC_MUT <- scores_filt_ord$MP_OPC_MUT
CARE_filt_rna_malignant$MP_CC_MUT <- scores_filt_ord$MP_CC_MUT
CARE_filt_rna_malignant$MP_MES_MUT <- scores_filt_ord$MP_MES_MUT
CARE_filt_rna_malignant$MP_AC1_MUT <- scores_filt_ord$MP_AC1_MUT

p_opc_mp <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "MP_OPC_MUT", embedding = "UMAPHarmony") +
  labs(x="", y="", color="OPC-like MP", title="") + theme(panel.border=element_blank()) +  guides(color = guide_legend(override.aes = list(size = 4)))
legend_grob <- g_legend(p_opc_mp)
p_opc_mp <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "MP_OPC_MUT", embedding = "UMAPHarmony") +
  labs(x="", y="", color="OPC-like MP", title="") + theme(panel.border=element_blank())  + guides(fill=FALSE)

plotPDF(p_opc_mp, name = "harmony_malignant_umap_opc.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "harmony_malignant_umap_opc.png"), p_opc_mp, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "harmony_malignant_umap_opc_legend.pdf"), legend_grob, width = 3, height = 3, dpi = 300)

p_ac_mp <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "MP_AC1_MUT", embedding = "UMAPHarmony") +
  labs(x="", y="", color="AC-like MP", title="") + theme(panel.border=element_blank()) +  guides(color = guide_legend(override.aes = list(size = 4)))
legend_grob <- g_legend(p_ac_mp)
p_ac_mp <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "MP_AC1_MUT", embedding = "UMAPHarmony") +
  labs(x="", y="", color="AC-like MP", title="") + theme(panel.border=element_blank())  + guides(fill=FALSE)

plotPDF(p_ac_mp, name = "harmony_malignant_umap_ac.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)
ggsave(paste0(fig_dir, "harmony_malignant_umap_ac.png"), p_ac_mp, width = 4, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "harmony_malignant_umap_ac_legend.pdf"), legend_grob, width = 3, height = 3, dpi = 300)


p_nfrags <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "nFrags", embedding = "UMAPHarmony")
p_tss <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAPHarmony")
ggAlignPlots(p_nfrags, p_tss, type = "h")


p_subtype <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "idh_codel_subtype", embedding = "UMAPHarmony",  pal = idh_cols) +
  labs(x="UMAP1", y="UMAP2", color="Subtype",title="Malignant cells (n = 70K)")

p_subtype_grade <- plotEmbedding(ArchRProj = CARE_filt_rna_malignant, colorBy = "cellColData", name = "subtype_grade", embedding = "UMAPHarmony") +
  labs(x="UMAP1", y="UMAP2", color="Subtype, Grade",title="Malignant cells (n = 70K)")


pdf(paste0(fig_dir, "caremut_atac_malignant_subtype_umap.pdf"), width = 4, height = 4, useDingbats = FALSE, bg = "transparent")
p_subtype
dev.off()

pdf(paste0(fig_dir, "caremut_atac_malignant_subtype_grade_umap.pdf"), width = 4, height = 4, useDingbats = FALSE, bg = "transparent")
p_subtype_grade
dev.off()

### ### ### ### ### ### ### ### ### ###
# Trajectory analyses
### ### ### ### ### ### ### ### ### ###

# AC >> Undifferentiated
ac_undiff_trajectory <- c("Undifferentiated", "MES", "AC")

CARE_filt_rna_malignant <- addTrajectory(
  ArchRProj = CARE_filt_rna_malignant, 
  name = "AC_Undiff", 
  groupBy = "CellStateGroup",
  trajectory = ac_undiff_trajectory, 
  embedding = "UMAPHarmony", 
  force = TRUE
)

ac_undiff_traj_plot <- plotTrajectory(CARE_filt_rna_malignant, trajectory = "AC_Undiff", colorBy = "cellColData", name = "AC_Undiff", embedding = "UMAPHarmony") 

pdf(paste0(fig_dir, "caremut_atac_malignant_ac_undifferentiated_trajectory_umap.pdf"), width = 4, height = 4, useDingbats = FALSE, bg = "transparent")
ac_undiff_traj_plot[[1]] + 
  labs(x="UMAP1", y="UMAP2", fill="Pseudotime",title="Undifferentiated to AC lineage trajectory")
dev.off()

CARE_filt_rna_malignant <- addImputeWeights(CARE_filt_rna_malignant)

traj_ac_undiff <- getTrajectory(ArchRProj = CARE_filt_rna_malignant, name = "AC_Undiff", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
# traj_ac_undiff <- getTrajectory(ArchRProj = CARE_filt_rna_malignant, name = "AC_Undiff", useMatrix = "PeakMatrix", log2Norm = TRUE)


# Possible markers to plot
markers_to_plot <- c("chr21:OLIG2", "chr12:ASCL1", "chr4:PDGFRA", "chr17:GFAP", "chr3:SOX2", "chr16:SOX8", "chr6:ID4","chr7:VGF")

ac_undiff_trajectory <- plotTrajectoryHeatmap(traj_ac_undiff,  
                                            #pal = paletteContinuous(set = "horizonExtra"),
                                            limits = c(-1.5, 1.5),
                                            labelTop = 25)

ac_undiff_trajectory <- plotTrajectoryHeatmap(traj_ac_undiff,  
                                               #pal = paletteContinuous(set = "horizonExtra"),
                                               limits = c(-1.5, 1.5),
                                               labelMarkers = markers_to_plot,
                                               labelTop = 1)

# https://github.com/GreenleafLab/ArchR/discussions/580
test <- getCellColData(CARE_filt_rna_malignant)[!is.na(getCellColData(CARE_filt_rna_malignant)$AC_Undiff),c("CellStateGroup", "AC_Undiff")]
as.matrix(test[order(test$AC_Undiff),"CellStateGroup"])
new_df <- data.frame(cbind(test[order(test$AC_Undiff),"CellStateGroup"], c(1:36337), rep(1,36337)))
new_df$X2 <- as.numeric(as.character(new_df$X2))
p <- ggplot(new_df, aes(X2, X3)) + geom_tile(aes(fill = X1)) + 
  scale_fill_manual(values= c("#AA2756", "#F77D58", "gray90")) + theme_ArchR() + theme(axis.text.x = element_text(angle = 90))



plotTrajectory(CARE_filt_rna_malignant, trajectory = "AC_Undiff", colorBy = "GeneScoreMatrix", name = "GFAP", continuousSet = "blueYellow")

pdf(paste0(fig_dir, "caremut_atac_malignant_undiff_ac_trajectory_heatmap.pdf"), width = 5, height = 5, useDingbats = FALSE, bg = "transparent")
ac_undiff_trajectory
dev.off()

opc_undiff_trajectory <- c("Undifferentiated", "OPC")

CARE_filt_rna_malignant <- addTrajectory(
  ArchRProj = CARE_filt_rna_malignant, 
  name = "OPC_Undiff", 
  groupBy = "CellStateGroup",
  trajectory = opc_undiff_trajectory, 
  embedding = "UMAPHarmony", 
  force = TRUE
)

plotTrajectory(CARE_filt_rna_malignant, trajectory = "OPC_Undiff", colorBy = "GeneScoreMatrix", name = "OLIG2", continuousSet = "blueYellow")

opc_undiff_traj_plot <- plotTrajectory(CARE_filt_rna_malignant, trajectory = "OPC_Undiff", colorBy = "cellColData", name = "OPC_Undiff", embedding = "UMAPHarmony") 

pdf(paste0(fig_dir, "caremut_atac_malignant_opc_undifferentiated_trajectory_umap.pdf"), width = 4, height = 4, useDingbats = FALSE, bg = "transparent")
opc_undiff_traj_plot[[1]] + 
  labs(x="UMAP1", y="UMAP2", fill="Pseudotime",title="Undifferentiated to OPC lineage trajectory")
dev.off()

 CARE_filt_rna_malignant <- addImputeWeights(CARE_filt_rna_malignant)

traj_opc_undiff <- getTrajectory(ArchRProj = CARE_filt_rna_malignant, name = "OPC_Undiff", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)

df_opc_undiff_trajectory <- assay(traj_opc_undiff)

# Possible markers to plot
markers_to_plot <- c("chr16:SOX8", "chr21:OLIG2", "chr4:PDGFRA", "chr16:RBFOX1", "chr9:NOTCH1", "chr18:TCF4", "chr3:SOX2", "chr21:CRYAA")

opc_undiff_trajectory <- plotTrajectoryHeatmap(traj_opc_undiff,  
                                              #pal = paletteContinuous(set = "horizonExtra"),
                                              limits = c(-1.5, 1.5),
                                              labelTop = 25)

opc_undiff_trajectory <- plotTrajectoryHeatmap(traj_opc_undiff,  
                                               #pal = paletteContinuous(set = "horizonExtra"),
                                               limits = c(-1.5, 1.5),
                                               labelMarkers = markers_to_plot,
                                               labelTop = 1)

pdf(paste0(fig_dir, "caremut_atac_malignant_undiff_opc_trajectory_heatmap.pdf"), width = 5, height = 5, useDingbats = FALSE, bg = "transparent")
opc_undiff_trajectory
dev.off()



# NPC >> OPC trajectory
opc_trajectory <- c("NPC", "OPC")

CARE_filt_rna_malignant <- addTrajectory(
  ArchRProj = CARE_filt_rna_malignant, 
  name = "StemU", 
  groupBy = "CellStateGroup",
  trajectory = opc_trajectory, 
  embedding = "UMAPHarmony", 
  force = TRUE
)

opc_traj_plot <- plotTrajectory(CARE_filt_rna_malignant, trajectory = "StemU", colorBy = "cellColData", name = "StemU", embedding = "UMAPHarmony") 

pdf(paste0(fig_dir, "caremut_atac_malignant_npc_opc_trajectory_umap.pdf"), width = 4, height = 4, useDingbats = FALSE, bg = "transparent")
opc_traj_plot[[1]] + 
  labs(x="UMAP1", y="UMAP2", fill="Pseudotime",title="NPC to OPC lineage trajectory")
dev.off()

CARE_filt_rna_malignant <- addImputeWeights(CARE_filt_rna_malignant)


trajStem <- getTrajectory(ArchRProj = CARE_filt_rna_malignant, name = "StemU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)

df_stem_trajectory <- assay(trajStem)

# Possible markers to plot
markers_to_plot <- c("chr16:SOX8", "chr2:GAD1", "chr21:OLIG2", "chr7:CNTNAP2", "chr4:PDGFRA", "chr16:RBFOX1", "chr9:NOTCH1", "chr18:TCF4")

npc_opc_trajectory <- plotTrajectoryHeatmap(trajStem,  
                            #pal = paletteContinuous(set = "horizonExtra"),
                            limits = c(-1.5, 1.5),
                            labelTop = 25)

pdf(paste0(fig_dir, "caremut_atac_malignant_npc_opc_trajectory_heatmap.pdf"), width = 5, height = 5, useDingbats = FALSE, bg = "transparent")
npc_opc_trajectory
dev.off()

plotPDF(p2, name = "Plot-UMAP-Malignant-Trajectory.pdf", ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE, width = 5, height = 5)



# NPC >> AC trajectory
trajectory <- c("NPC", "AC")

CARE_filt_rna_malignant <- addTrajectory(
  ArchRProj = CARE_filt_rna_malignant, 
  name = "AC_U", 
  groupBy = "CellStateGroup",
  trajectory = trajectory, 
  embedding = "UMAPHarmony", 
  force = TRUE
)

p_traj_npc_ac <- plotTrajectory(CARE_filt_rna_malignant, trajectory = "AC_U", colorBy = "cellColData", name = "AC_U", continuousSet = "horizonExtra",  embedding = "UMAPHarmony")
p_traj_npc_ac[[1]]

pdf(paste0(fig_dir, "caremut_atac_malignant_npc_ac_trajectory_umap.pdf"), width = 4, height = 4, useDingbats = FALSE, bg = "transparent")
p_traj_npc_ac[[1]] + 
  labs(x="UMAP1", y="UMAP2", fill="Pseudotime",title="NPC to AC lineage trajectory")
dev.off()

traj_ac <- getTrajectory(ArchRProj = CARE_filt_rna_malignant, name = "AC_U", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)

npc_ac_trajectory <- plotTrajectoryHeatmap(traj_ac,  
                                           #pal = paletteContinuous(set = "horizonExtra"),
                                           limits = c(-1.5, 1.5),
                                           labelTop = 10)


pdf(paste0(fig_dir, "caremut_atac_malignant_npc_ac_trajectory_heatmap.pdf"), width = 5, height = 5, useDingbats = FALSE, bg = "transparent")
npc_ac_trajectory
dev.off()


## A few QC plots for the cohort.
df <- getCellColData(CARE_filt_rna_malignant, select = c("log10(nFrags)", "TSSEnrichment"))
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

pdf(paste0(fig_dir, "malignant_cells_density_nfrags_vs_tss.pdf"), width = 5, height = 5, useDingbats = FALSE, bg = "transparent")
p
dev.off()

p1 <- plotGroups(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "patient_id", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  maxCells = 10000,
) + labs(x="")

table(CARE_filt_rna_malignant$patient_id, CARE_filt_rna_malignant$Sample)
pdf(paste0(fig_dir, "per_patient_tss_scores.pdf"), width = 8, height = 5, useDingbats = FALSE, bg = "transparent")
p1
dev.off()


p2 <- plotGroups(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "patient_id", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  maxCells = 10000
) 

table(CARE_filt_rna_malignant$patient_id, CARE_filt_rna_malignant$Sample)
pdf(paste0(fig_dir, "per_patient_nfrags.pdf"), width = 8, height = 5, useDingbats = FALSE, bg = "transparent")
p2
dev.off()

pdf(paste0(fig_dir, "per_patient_tss_vs_nfrags.pdf"), width = 8, height = 5, useDingbats = FALSE, bg = "transparent")
ggAlignPlots(p2, p1, type = "v")
dev.off()


p1 <- plotGroups(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "patient_id", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  maxCells = 10000,
) + labs(x="")

table(CARE_filt_rna_malignant$patient_id, CARE_filt_rna_malignant$Sample)
pdf(paste0(fig_dir, "per_patient_tss_scores.pdf"), width = 8, height = 5, useDingbats = FALSE, bg = "transparent")
p1
dev.off()


p2 <- plotGroups(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "patient_id", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  maxCells = 10000
) 

table(CARE_filt_rna_malignant$patient_id, CARE_filt_rna_malignant$Sample)
pdf(paste0(fig_dir, "per_patient_nfrags.pdf"), width = 8, height = 5, useDingbats = FALSE, bg = "transparent")
p2
dev.off()

pdf(paste0(fig_dir, "per_patient_tss_vs_nfrags.pdf"), width = 8, height = 5, useDingbats = FALSE, bg = "transparent")
ggAlignPlots(p2, p1, type = "v")
dev.off()

# Per sample
p1 <- plotGroups(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "care_id", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  maxCells = 10000,
) + labs(x="")

pdf(paste0(fig_dir, "per_sample_tss_scores.pdf"), width = 8, height = 5, useDingbats = FALSE, bg = "transparent")
p1
dev.off()


p2 <- plotGroups(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "care_id", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  maxCells = 10000
) 

table(CARE_filt_rna_malignant$patient_id, CARE_filt_rna_malignant$Sample)
pdf(paste0(fig_dir, "per_sample_nfrags.pdf"), width = 8, height = 5, useDingbats = FALSE, bg = "transparent")
p2
dev.off()

pdf(paste0(fig_dir, "per_sample_tss_vs_nfrags.pdf"), width = 8, height = 5, useDingbats = FALSE, bg = "transparent")
ggAlignPlots(p2, p1, type = "v")
dev.off()


### END ###