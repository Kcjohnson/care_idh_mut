##############################
### Examine the differences in malignant cell states across IDH-WT and IDH-MUT glioma
### Author: Kevin Johnson
### Updated: 2025.10.14
##############################

## Part 10: Examine whether IDH-WT tumors adopt a more open chromatin confirmation

## ArchR creates several directories automatically when creating arrow files and ArchR projects.
workdir <- "/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/"
setwd(workdir)

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
cell_class <- read.table("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/processed_data/rna/classification/caremut_all_select_state_assignment_20240416.txt",  sep = "\t", header = TRUE)
care_mut_mal_md <- cell_class %>% 
  mutate(idh_codel_subtype = recode(idh_codel_subtype, `IDHmut-codel` = "IDH-O",
                                    `IDHmut-noncodel` = "IDH-A"),
         State = recode(State, `MP_AC1_MUT` = "AC-like",
                        `MP_OPC_MUT` = "OPC-like",
                        `MP_NPC_MUT` = "NPC-like",
                        `MP_MES_MUT` = "MES-like",
                        `MP_AC2_MUT` = "AC-like",
                        "Undifferentiated" = "Undifferentiated"),
         idh_status = "IDHmut") %>% 
  dplyr::select(CellID, care_id, idh_codel_subtype, idh_status, State)

# Read in the RNA samples from Synapse
care_wt_mal_md <-readRDS("/vast/palmer/pi/verhaak/kcj28/care_wt/github/data/malignant_meta_data_2025_01_08.RDS")

care_wt_mal_md <- care_wt_mal_md %>% 
  mutate(idh_codel_subtype = "IDHwt",
         idh_status = "IDHwt", 
         CellID = sub("^([A-Z]{2}\\d{2})(0-)", "\\1-\\2", CellID)) %>% 
  dplyr::select(CellID, care_id = ID, idh_codel_subtype, idh_status, State)

all_care_mal_md <- care_wt_mal_md %>% 
  bind_rows(care_mut_mal_md)

## Create a data.frame that contains the essential information to be merged and inspected.
atac_df <- data.frame(projMONITOR_filt$cellNames, projMONITOR_filt@cellColData$Sample, projMONITOR_filt@cellColData$TSSEnrichment)
table(atac_df$projMONITOR_filt.cellColData.Sample)
# The cell names are a little different between RNA (Seurat) and ATAC (ArchR). Need to create a common linker.
atac_df$CellID <- gsub("#", "-", atac_df$projMONITOR_filt.cellNames)

# 81,527 cells present in ATACseq data that pass doublet removal and are also found in snRNAseq data that passed QC for IDH-mutant.
sum(atac_df$CellID%in%all_care_mal_md$CellID)

# There's only 10,439 malignant cells in both the ATAC data and passed QC for IDH-WT CARE.
sum(atac_df$CellID%in%care_mut_mal_md$CellID)
sum(atac_df$CellID%in%care_wt_mal_md$CellID)

care_wt_mal_md_atac <- care_wt_mal_md %>% 
  filter(CellID%in%atac_df$CellID)

# Combine the two data.frames by adding on the RNA data and filtering out what's leftover.
atac_df_filt_rna <- atac_df %>% 
  inner_join(all_care_mal_md, by="CellID") 


## Identify which cells to keep in the analysis. 
tmp <- getCellNames(projMONITOR_filt)
rna_cells = tmp[which(tmp%in%atac_df_filt_rna$projMONITOR_filt.cellNames)]  

## Subset the cells to only those that also have RNA.
CARE_all_rna_malignant <- subsetCells(ArchRProj = projMONITOR_filt, cellNames = rna_cells)


all(atac_df_filt_rna$projMONITOR_filt.cellNames==CARE_all_rna_malignant$cellNames)

CARE_all_rna_malignant$CellID <- atac_df_filt_rna$CellID
CARE_all_rna_malignant$care_id <- atac_df_filt_rna$care_id
CARE_all_rna_malignant$idh_codel_subtype <- atac_df_filt_rna$idh_codel_subtype
CARE_all_rna_malignant$idh_status <- atac_df_filt_rna$idh_status
CARE_all_rna_malignant$State <- atac_df_filt_rna$State
CARE_all_rna_malignant$type_state <- paste0(atac_df_filt_rna$idh_status, "-", atac_df_filt_rna$State)

markerState_wt <- getMarkerFeatures(
  ArchRProj = CARE_all_rna_malignant, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "idh_status",
  testMethod = "wilcoxon",
  maxCells = 5000,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "IDHwt", 
  bgdGroups = "IDHmut"
)

upregulated <- getMarkers(markerState_wt, cutOff = "FDR <= 0.05 & Log2FC >= 1")
downregulated <- getMarkers(markerState_wt, cutOff = "FDR <= 0.05 & Log2FC <= -1")
upregulated_genes <- as.data.frame(upregulated$`IDHwt`)
downregulated_genes <- as.data.frame(downregulated$`IDHwt`)
dag_genes_wt <- bind_rows(upregulated_genes, downregulated_genes) %>% 
  mutate(state = "IDHwt", 
         direction = ifelse(Log2FC>1, "upregulated", "downregulated"))

table(dag_genes_wt$direction, dag_genes_wt$seqnames)


markersGS <- getMarkerFeatures(
  ArchRProj = CARE_all_rna_malignant, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "idh_codel_subtype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

wt <- markerList$IDHwt
noncodel <- markerList$`IDH-A`
codel <- markerList$`IDH-O`$name

markerGenes  <- c(
  "EGFR", # IDH-WT
  "HOXA1",
  "HOXA4",
  "HOXA5",
  "HIST1H3E",
  "SOX10", # IDH-MUT
  "SOX1")

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

markerState_idh <- getMarkerFeatures(
  ArchRProj = CARE_all_rna_malignant, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "idh_status",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmap_idh <- plotMarkerHeatmap(
  seMarker = markerState_idh, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE,
  plotLog2FC = TRUE
)

ComplexHeatmap::draw(heatmap_idh, heatmap_legend_side = "bot", annotation_legend_side = "bot")


marker_idh_state <- getMarkerFeatures(
  ArchRProj = CARE_all_rna_malignant, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "type_state",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(marker_idh_state, cutOff = "FDR <= 0.05 & Log2FC >= 1")

ac_wt <- markerList$`IDHwt-AC`
noncodel <- markerList$`IDH-A`
codel <- markerList$`IDH-O`$name

markerGenes  <- c(
  "EGFR", # IDH-WT
  "HOXA1",
  "HOXA4",
  "HOXA5",
  "HIST1H3E",
  "SOX10", # IDH-MUT
  "SOX1")


heatmap_idh <- plotMarkerHeatmap(
  seMarker = marker_idh_state, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  #labelMarkers = markerGenes,
  transpose = TRUE,
  plotLog2FC = TRUE
)

ComplexHeatmap::draw(heatmap_idh, heatmap_legend_side = "bot", annotation_legend_side = "bot")

