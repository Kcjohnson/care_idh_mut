##############################
### Run ArchR peak calling analysis on CAREmut multiome ATAC data for all malignant states - NPC, OPC, Undiff, MES, AC
### Author: Kevin Johnson
### Updated: 2024.05.05
##############################

## Part 3: Peak calling on RNA-defined MALIGNANT states across IDH-mutant tumors

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
# This project comes from: `02_run_archr_caremut_all_cell_types.R`
CARE_filt_rna_all <- loadArchRProject("Save-CAREmut-All-RNA")
getAvailableMatrices(CARE_filt_rna_all) # GeneScoreMatrix and TileMatrix; 118,268 nuclei

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# CAREmut data that was processed by Seurat. Here is the associated metadata.
mut_md <- read.table("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20240212.txt", sep = "\t", header = TRUE)
mut_md_verhaak <- mut_md %>% 
  filter(lab=="Verhaak lab")

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

# There is a slight difference in ArchR and Seurat naming convention.
CARE_filt_rna_all$CellID <-  gsub("#", "-", CARE_filt_rna_all$cellNames)

# Restricting to malignant cells and the malignant state classification.
cell_class_filt <- mut_md_verhaak_state[mut_md_verhaak_state$CellID %in% CARE_filt_rna_all$CellID, ]
cell_class_filt_ord <- cell_class_filt[match(CARE_filt_rna_all$CellID, cell_class_filt$CellID), ]
all(CARE_filt_rna_all$CellID==cell_class_filt_ord$CellID)
CARE_filt_rna_all$CellStateGroup <- cell_class_filt_ord$group

# Subset entire CARE IDHmut object to only MALIGNANT cells.
idxSample <- BiocGenerics::which(CARE_filt_rna_all$CellType_final%in%"Malignant")
cellsSample <- CARE_filt_rna_all$cellNames[idxSample]

# Create a subset of the ArchR object to save malignant-only analyses.
CARE_filt_rna_malignant <- subsetArchRProject(
  ArchRProj = CARE_filt_rna_all,
  cells = cellsSample,
  outputDirectory = "Save-CAREmut-Malignant-RNA",
  force = TRUE)

# Sanity check that it was correctly loaded - otherwise it may write results to old directory.
CARE_filt_rna_malignant <- loadArchRProject("Save-CAREmut-Malignant-RNA")

# Plot frags against TSS enrichment for malignant cells only - what's the median values? 9.138 TSS and 13,920 fragments.
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

# Quickly assess whether there are any differences in QC metrics for cells that were not found to be concordant for their cluster identity.
malignant_archr_df <- data.frame(CARE_filt_rna_malignant@cellColData)
malignant_archr_df$concordance_state <- ifelse(malignant_archr_df$CellType_ATAC==malignant_archr_df$CellType_final, "concordant", "discordant")
table(malignant_archr_df$concordance_state, malignant_archr_df$idh_codel_subtype)

# Slightly lower TSS and nFrags for discordant clusters. 
ggplot(malignant_archr_df, aes(x=concordance_state, y=TSSEnrichment)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(.~idh_codel_subtype) 

ggplot(malignant_archr_df, aes(x=concordance_state, y=nFrags)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(.~idh_codel_subtype) 


# Add grade information from clinical tables
patient_md <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240227.txt", sep="\t", header = TRUE)
sample_md <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240227.txt", sep="\t", header = TRUE)

atac_md_mut <- data.frame(getCellColData(CARE_filt_rna_malignant))
atac_md_mut_annot <- atac_md_mut %>% 
  left_join(sample_md, by=c("care_id", "sample_barcode", "patient_id", "timepoint", "idh_codel_subtype")) %>% 
  mutate(atacCellNames = rownames(atac_md_mut))

# Do these cell names retain the same order?
ifelse(all(getCellNames(CARE_filt_rna_malignant)==atac_md_mut_annot$atacCellNames),
       sprintf("All cell names match. Proceed"), sprintf("Warning! Cell names do not match!"))

# Add a few RNA features to the ArchR object.
CARE_filt_rna_malignant$Grade <- paste0("G", atac_md_mut_annot$grade_num)
CARE_filt_rna_malignant$hypermutation <- atac_md_mut_annot$hypermutation
CARE_filt_rna_malignant$subtype_grade <- paste0(CARE_filt_rna_malignant$idh_codel_subtype, "_", CARE_filt_rna_malignant$Grade)

# What's the breakdown of sample information: mostly codel malignant cells (50K vs 20K)
table(CARE_filt_rna_malignant$idh_codel_subtype)
table(CARE_filt_rna_malignant$Sample, CARE_filt_rna_malignant$subtype_grade)
table(CARE_filt_rna_malignant$Sample, CARE_filt_rna_malignant$CellStateGroup)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Make pseudo bulk measurements for Cell state groups
# See github discussion on how to structure groups: https://github.com/GreenleafLab/ArchR/discussions/696
# Pseudo-bulk refers to a grouping of single cells where the data from each single cell is combined into a single pseudo-sample.
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
# Account for biases in data quality via TSSEnrichment and nFrags - these are the default settings on the tutorial and seem advisable.
markersPeaks <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "PeakMatrix", 
  groupBy = "CellStateGroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 2000
)

saveRDS(markersPeaks, "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/malignant_cell_state_markersPeaks_20240505.RDS")
markersPeaks <- readRDS("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/malignant_cell_state_markersPeaks_20240505.RDS")

# Extract the marker peaks. Get access the GRanges object via `returnGR = TRUE`.
# Throughout these analyses, I am sticking with the FDR <= 0.05 and Log2FC >= 1 so that it is simple.
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)
lengths(markerList) # AC = 33454, MES = 14250, NPC = 2309, OPC = 12355, Undifferentiated = 5829

table(markerList$Undifferentiated@seqnames)
table(markerList$OPC@seqnames)
table(markerList$NPC@seqnames)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = FALSE,
  nLabel = 1
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-mut-malignant-state", width = 8, height = 6, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

pdf(paste0(fig_dir, "caremut_malignant_state_marker_peaks.pdf"), width = 6, height = 4, useDingbats = FALSE, bg = "transparent")
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


markerGenes <- c("AQP4", "CD44", 
                 "DLL3", "PDGFRA",
                 "HOXD11", "MET",
                 "HOXB2")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)


CARE_filt_rna_malignant <- saveArchRProject(ArchRProj = CARE_filt_rna_malignant, outputDirectory = "Save-CAREmut-Malignant-RNA", load = TRUE) 

p <- plotBrowserTrack(
  ArchRProj = CARE_filt_rna_malignant, 
  groupBy = "CellStateGroup", 
  geneSymbol = c("PDGFRA"),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)["Undifferentiated"],
  upstream = 50000,
  downstream = 50000
)

grid::grid.draw(p$PDGFRA)
plotPDF(p, name = "Plot-Tracks-With-PDGFRA", width = 5, height = 5, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)
############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
## Add add motif annotations.
if("Motif" %ni% names(CARE_filt_rna_malignant@peakAnnotation)){
  CARE_filt_rna_malignant <- addMotifAnnotations(ArchRProj = CARE_filt_rna_malignant, motifSet = "cisbp", name = "Motif")
}

# For motifs amongst the open chromatin peak regions.
enrichRegions <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CARE_filt_rna_malignant,
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

plotPDF(heatmapRegions, name = "Regions-Enriched-Marker-Peak-Motif-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

heatmapATAC_df <- plotEnrichHeatmap(enrichRegions, 
                                    n = 1500, 
                                    transpose = TRUE,
                                    returnMatrix = TRUE)

df <- data.frame(TF = rownames(enrichRegions), mlog10Padj = assay(enrichRegions))


library(viridisLite)
colnames(heatmapATAC_df) <- sapply(strsplit(colnames(heatmapATAC_df), " "), "[[", 1)
colnames(heatmapATAC_df) <- sapply(strsplit(colnames(heatmapATAC_df), "_"), "[[", 1)
enrichment_of_interest <- c("TCF12", "ASCL1", "CREB5", "TAL1", "JUNB", "FOS", "NFIC","SOX9")
feature_order <- c("OPC", "NPC", "Undifferentiated", "MES", "AC")
enrichment_order <- c("TCF12", "ASCL1", "TAL1", "CREB5", "JUNB", "FOS", "NFIC","SOX9")

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

motif_dev_df <- getVarDeviations(CARE_filt_rna_malignant, name = "MotifMatrix", plot = TRUE)

## How can one extract a subset of motifs for downstream analyses? getFeatures()
motifs <- c("ASCL1", "TCF12", "JUNB", "FOS", "RFX2", "TAL1", "NFIC")

markerMotifs <- getFeatures(CARE_filt_rna_malignant, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- markerMotifs[grep("z:", markerMotifs)]


p <- plotGroups(ArchRProj = CARE_filt_rna_malignant, 
                groupBy = "CellStateGroup", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(CARE_filt_rna_malignant)
)

plotPDF(p, name = "Plot-State-Motifs-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

## Trying to extract the relevant TFs from the MotifMatrix. Try Z-scores 
motif_df <- getMatrixFromProject(
  ArchRProj = CARE_filt_rna_malignant,
  useMatrix = "MotifMatrix",
  useSeqnames = "z",
  verbose = TRUE,
  binarize = FALSE
)

motif_df_zscore <- assay(motif_df)
motif_df_zscore_out <- as.data.frame(as.matrix(t(motif_df_zscore)))
rownames(motif_df_zscore_out) <-  gsub("#", "-", rownames(motif_df_zscore_out))
motif_df_zscore_out$CellID <- rownames(motif_df_zscore_out)

# Write out the results so that they can be used downstream. 
write.table(motif_df_zscore_out, "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/archr_care_celltype_tf_motif_activity_zscore_20240505.txt", sep="\t", row.names = FALSE, col.names = TRUE)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
### Plot ENCODE enrichment
CARE_filt_rna_malignant <- addArchRAnnotations(ArchRProj = CARE_filt_rna_malignant, collection = "EncodeTFBS")

enrichEncode <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CARE_filt_rna_malignant,
  peakAnnotation = "EncodeTFBS",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

df <- data.frame(TF = rownames(enrichEncode), mlog10Padj = assay(enrichEncode)[,4])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

# These are cell type relevant ENCODE lines: NH-A = human astrocytes, hESC = human embryonic stem cell
indices <- grep("NH_A|hESC", df$TF)

# Subset to the proteins we care about
df_sub <- df[indices, ]

ggplot(df_sub, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df_sub[rev(seq_len(5)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 40, transpose = TRUE)

plotPDF(heatmapEncode, name = "EncodeTFBS-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)


df <- data.frame(TF = rownames(enrichEncode), mlog10Padj = assay(enrichEncode)[,5])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

enrichEncode_out <- assay(enrichEncode)
enrichEncode_out$TF <- rownames(enrichEncode_out)

write.table(enrichEncode_out, "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/archr_care_mut_malignant_tf_enrichment_amongst_markerpeaks_20240505.txt", sep="\t", row.names = FALSE, col.names = TRUE)



#enrichEncode_out <- assay(enrichEncode)[1]
# Use unadjusted for now
enrichEncode_out <- enrichEncode@assays@data[[1]]
enrichEncode_out$TF <- rownames(enrichEncode_out)

indices <- grep("NH_A|hESC", enrichEncode_out$TF)
enrichEncode_out_sub <- enrichEncode_out[indices, ]
enrichEncode_out_sub$TF <- gsub("\\.\\.\\.", "", enrichEncode_out_sub$T)

enrichEncode_out_sub$TF_listed <- sapply(strsplit(enrichEncode_out_sub$TF, "^\\d+\\."), "[[", 2)
enrichEncode_out_sub$TF_listed <- gsub("_39-|-H1_|_39-H1_", "_", enrichEncode_out_sub$TF_listed)
enrichment_of_interest <- c("EZH2_hESC", "EZH2_NH_A", "TCF12_hESC", "SUZ12_hESC")
feature_order <- c("OPC", "NPC", "Undifferentiated", "MES", "AC")
enrichment_order <- c("EZH2_NH_A", "EZH2_hESC", "SUZ12_hESC", "TCF12_hESC")

enrichEncode_out_sub_filtered <- enrichEncode_out_sub[enrichEncode_out_sub$TF_listed%in%enrichment_of_interest,]
rownames(enrichEncode_out_sub_filtered) <- enrichEncode_out_sub_filtered$TF_listed
enrichEncode_final <- enrichEncode_out_sub_filtered[,1:5]
enrichEncode_final_ordered <- enrichEncode_final[enrichment_order, feature_order]

library(ComplexHeatmap)
encode_hmap <- ComplexHeatmap::Heatmap(enrichEncode_final_ordered,                    
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


pdf(paste0(fig_dir, "malignant_states_archr_differential_peak_enrichment_encode_heatmap.pdf"), width = 5, height = 5, useDingbats = FALSE, bg = "transparent")
draw(encode_hmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
## Compute the per-cell deviations across all of our motif annotations.
## This function has an optional parameter called matrixName that allows us to define the name of deviations.
## The option below creates a deviation matrix in each of the Arrow files called "DevMatrix" or MotifMatrix. Force indicates whether the matrix listed should be overwritten.
CARE_filt_rna_malignant <- addDeviationsMatrix(
  ArchRProj = CARE_filt_rna_malignant, 
  peakAnnotation = "EncodeTFBS",
  matrixName = "EncodeMatrix",
  threads = getArchRThreads(),
  force = TRUE
)

encode_dev_df <- getVarDeviations(CARE_filt_rna_malignant, name = "EncodeMatrix", plot = TRUE)

key_tfs <- c("EZH2", "SUZ12", "JUN", "FOS")

encode_tf_df_zscore <- assay(encode_tf_df)

markerEncode <- getFeatures(CARE_filt_rna_malignant, select = paste(key_tfs, collapse="|"), useMatrix = "EncodeMatrix")
markerEncode <- markerEncode[grep("z:", markerEncode)]


p <- plotGroups(ArchRProj = CARE_filt_rna_malignant, 
                groupBy = "CellStateGroup", 
                colorBy = "EncodeMatrix", 
                name = markerEncode,
                imputeWeights = getImputeWeights(CARE_filt_rna_malignant)
)

plotPDF(p, name = "Plot-Groups-Encode-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

## Trying to extract the relevant TFs from the MotifMatrix. Try Z-scores 
encode_tf_df <- getMatrixFromProject(
  ArchRProj = CARE_filt_rna_malignant,
  useMatrix = "EncodeMatrix",
  useSeqnames = "z",
  verbose = TRUE,
  binarize = FALSE
)

encode_tf_df_zscore <- assay(encode_tf_df)

encode_tf_df_zscore_out <- as.data.frame(as.matrix(t(encode_tf_df_zscore)))
rownames(encode_tf_df_zscore_out) <-  gsub("#", "-", rownames(encode_tf_df_zscore_out))
encode_tf_df_zscore_out$CellID <- rownames(encode_tf_df_zscore_out)

cell_class_encode_tf <- cell_class %>% 
  inner_join(encode_tf_df_zscore_out, by="CellID")

write.table(cell_class_encode_tf, "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/archr_care_mut_malignant_tf_encode_activity_zscore_20240505.txt", sep="\t", row.names = FALSE, col.names = TRUE)

# Due to the enrichment of EZH2, check whether the atac module gene activity score for EZH2 target genes is consistent
tmp <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/reference/archr/BENPORATH_PRC2_TARGETS.v2023.2.Hs.tsv", sep="\t", header = TRUE)
prc2_targets <- tmp[17, 2]

all_genes <- getFeatures(CARE_filt_rna_malignant)
prc2_targets <- as.character(unlist(strsplit(prc2_targets, ",")))
prc2_targets_filt <- prc2_targets[prc2_targets%in%all_genes]
prc2_features <- list(
  PRC2score = prc2_targets_filt
)

CARE_filt_rna_malignant <- addModuleScore(CARE_filt_rna_malignant, features = prc2_features, useMatrix = "GeneScoreMatrix")

df_prc2 <- data.frame(CARE_filt_rna_malignant@cellColData)

# Higher levels between the 
ggplot(df_prc2, aes(x=CellStateGroup, y=Module.PRC2score)) +
  geom_boxplot() +
  theme_bw() +
  facet_grid(.~idh_codel_subtype)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
### Plot ATAC enrichment 
CARE_filt_rna_malignant <- addArchRAnnotations(ArchRProj = CARE_filt_rna_malignant, collection = "ATAC")

enrichATAC <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CARE_filt_rna_malignant,
  peakAnnotation = "ATAC",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 5, cutOff = 3.5, transpose = TRUE)

plotPDF(heatmapATAC, name = "ATAC-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

df <- data.frame(TF = rownames(enrichATAC), mlog10Padj = assay(enrichATAC)[,3])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

# Interestingly, NPC cells were enriched for Excitatory neuron profiles
head(df)

enrichATAC_out <- assay(enrichATAC)
enrichATAC_out$Bulk <- rownames(enrichATAC_out)

write.table(enrichATAC_out, "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/archr_care_mut_malignant_bulk_atac_enrichment_amongst_markerpeaks_20240505.txt", sep="\t", row.names = FALSE, col.names = TRUE)


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Vierstra  produced "motif archetypes" which represent clustered motifs that have been essentially deduplicated based on similarity.
# Original approach
motifPWMs <- readRDS("/gpfs/gibbs/pi/verhaak/kcj28/reference/archr/Vierstra-Human-Motifs.rds")
CARE_filt_rna_malignant <- addMotifAnnotations(ArchRProj = CARE_filt_rna_malignant, motifPWMs = motifPWMs, annoName = "Vierstra", force = TRUE)
# Natively supported approach
CARE_filt_rna_malignant <- addMotifAnnotations(ArchRProj = CARE_filt_rna_malignant, motifSet="vierstra", collection="archetype", annoName = "Vierstra_v2", force = TRUE)

# I didn't notice a large difference in these results versus the ones from regular motif
enrichViestra <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CARE_filt_rna_malignant,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

heatmapViestra <- plotEnrichHeatmap(enrichViestra, n = 10, cutOff = 10, transpose = TRUE)

plotPDF(heatmapViestra, name = "Viestra-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_malignant, addDOC = FALSE)

df <- data.frame(TF = rownames(enrichViestra), mlog10Padj = assay(enrichViestra)[,5])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Pairwise differential expression. In the snRNA data, there was a clear anti-correlation for OPC-like and AC-like cells. Assessing that difference here.
markerOPC <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "PeakMatrix",
  groupBy = "CellStateGroup",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "OPC",
  bgdGroups = "AC"
)

# Approximately equal up- and down-regulated. Total peaks: 284568 (6798 up and 8642 down)
volcano_opc_v_ac <- plotMarkers(seMarker = markerOPC, name = "OPC", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "Volcano")
volcano_opc_v_ac

markerNPC_v_OPC <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "PeakMatrix",
  groupBy = "CellStateGroup",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "NPC",
  bgdGroups = "OPC"
)

# Limited differences. Total peaks: 284568 (69 up and 285 down)
volcano_npc_v_opc <- plotMarkers(seMarker = markerNPC_v_OPC, name = "NPC", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "Volcano")
volcano_npc_v_opc

markerMES_v_AC <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "PeakMatrix",
  groupBy = "CellStateGroup",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "MES",
  bgdGroups = "AC"
)

# Approximately equal up- and down-regulated. Total peaks: 284568 ( up and  down)
volcano_mes_v_ac <- plotMarkers(seMarker = markerMES_v_AC, name = "MES", cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "Volcano")
volcano_mes_v_ac


MesMotifsUp <- peakAnnoEnrichment(
  seMarker = markerMES_v_AC,
  ArchRProj = CARE_filt_rna_malignant,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

df <- data.frame(TF = rownames(MesMotifsUp), mlog10Padj = assay(MesMotifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))


# Save output where new images will be deposited.
paste0("Memory Size = ", round(object.size(CARE_filt_rna_malignant) / 10^6, 3), " MB")
CARE_filt_rna_malignant <- saveArchRProject(ArchRProj = CARE_filt_rna_malignant, outputDirectory = "Save-CAREmut-Malignant-RNA", load = TRUE) 

malignant_archr_df <- data.frame(CARE_filt_rna_malignant@cellColData)
write.table(malignant_archr_df, "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/archr_care_mut_malignant_state_atac_metadata_20240505.txt", sep="\t", row.names = FALSE, col.names = TRUE)



CARE_filt_rna_malignant <- loadArchRProject("Save-CAREmut-Malignant-RNA")

### END ###