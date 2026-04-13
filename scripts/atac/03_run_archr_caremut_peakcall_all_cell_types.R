##############################
### Run ArchR peak calling analysis on CAREmut multiome ATAC data for all cell types
### Author: Kevin Johnson
### Updated: 2026.04.08
##############################

## Part 3: Peak calling on all RNA-defined CELL TYPES across IDH-mutant tumors

# Specify directories
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

# Load the ArchR object for IDHmut final analysis set. Doublets have been removed, only cells that intersect with passed qc for RNA, and RNA annotated cell states.
CARE_filt_rna_all <- loadArchRProject("Save-CAREmut-All-RNA")

# Remove "Unresolved" (RNA-based) cells from further analysis since it may cause issues. Unresolved are mostly malignant/astrocytes, which seem to be most prone to misclassification in glioma.
atac_df <- data.frame(CARE_filt_rna_all@cellColData)

atac_df_filt <- atac_df %>% 
  filter(CellType_final!="Unresolved")
current_rna_atac_cells <- getCellNames(CARE_filt_rna_all)
rna_cells_to_keep = current_rna_atac_cells[which(current_rna_atac_cells%in%rownames(atac_df_filt))]  

# Subset the cells to only those that also have RNA.
CARE_filt_rna_all_filt <- subsetCells(ArchRProj = CARE_filt_rna_all, cellNames = rna_cells_to_keep)

# Add grade information
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt", sep="\t", header = TRUE)
sample_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240608.txt", sep="\t", header = TRUE)

atac_md_mut <- data.frame(getCellColData(CARE_filt_rna_all_filt))
atac_md_mut_annot <- atac_md_mut %>% 
  left_join(sample_md, by=c("care_id", "sample_barcode", "patient_id", "timepoint", "idh_codel_subtype")) %>% 
  mutate(atacCellNames = rownames(atac_md_mut))

## Do these cell names retain the same order?
ifelse(all(getCellNames(CARE_filt_rna_all_filt)==atac_md_mut_annot$atacCellNames),
       sprintf("All cell names match. Proceed"), sprintf("Warning! Cell names do not match!"))

## Add a few RNA features to the ArchR object.
CARE_filt_rna_all_filt$Grade <- paste0("G", atac_md_mut_annot$grade_num)
CARE_filt_rna_all_filt$hypermutation <- atac_md_mut_annot$hypermutation
CARE_filt_rna_all_filt$subtype_grade <- paste0(CARE_filt_rna_all_filt$idh_codel_subtype, "_", CARE_filt_rna_all_filt$Grade)

# What's the breakdown of sample information:83K (Oligo.) vs 34K (Astro.)
table(CARE_filt_rna_all_filt$idh_codel_subtype)
table(CARE_filt_rna_all_filt$sample_barcode, CARE_filt_rna_all_filt$timepoint)
table(CARE_filt_rna_all_filt$patient_id, CARE_filt_rna_all_filt$timepoint, CARE_filt_rna_all_filt$idh_codel_subtype)
table(CARE_filt_rna_all_filt$patient_id, CARE_filt_rna_all_filt$timepoint, CARE_filt_rna_all_filt$idh_codel_subtype)
table(CARE_filt_rna_all_filt$idh_codel_subtype, CARE_filt_rna_all_filt$Grade)
table(CARE_filt_rna_all_filt$patient_id, CARE_filt_rna_all_filt$Grade)
table(CARE_filt_rna_all_filt$CellType_final)

# 117,173 cells
CARE_filt_rna_all_filt <- saveArchRProject(ArchRProj = CARE_filt_rna_all_filt, 
                 outputDirectory = "Save-CAREmut-All-RNA-Filtered", 
                 load = TRUE, 
                 dropCells = TRUE,
                 overwrite = TRUE) 


############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Assess differential gene activity scores across the different RNA-based cell types.
#  devtools::install_github('immunogenomics/presto', repos = BiocManager::repositories())

## Define markers based on RNA-cell type assignment.
markersGS <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_all_filt, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "CellType_final",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  # Note we are increasing the default number of maxCells to 2000 (some cell types do not have this ammount: lymphocytes (304), endothelial (1666), mural (19990).
  # We could increase to higher number, but then there is an imbalance with the rarer TME cell types.
  maxCells = 2000,
)

## Extract the marker list. These are the default thresholds for ArchR.
# Throughout this analysis, I set a cutoff criteria of FDR < 0.05 and Log2FC >= 1. The ArchR tutorial uses variable cut-offs and I could not find a good explanation for why.
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
saveRDS(markerList, "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/archr_all_celltypes_markerlist.RDS")

# Examine the distribution of differentially accessible genes
lapply(markerList, nrow)
markerList$Astrocyte
markerList$Malignant # PTPRZ1 and SOX6 towards the top of the list.
markerList$Oligodendrocyte
markerList$Mural
markerList$Endothelial
markerList$Lymphocyte
markerList$Myeloid
markerList$InhNeuron
markerList$ExcNeuron

# Selected RNA marker genes
markerGenes <- c("PTPRZ1", "SOX6", "AQP4", "SLC1A2", "MBP", "PLP1", "RBFOX1", "SNAP25", "GAD1", "GAD2", "VWF", "FLT1", "COL1A2", "DCN", "CSF1R", "MSR1", "CD247", "CLEC2D")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", 
  labelMarkers = markerGenes,
  transpose = FALSE)

heatmap_out <- ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap-All-Cells", width = 5, height = 5, ArchRProj = CARE_filt_rna_all_filt, addDOC = FALSE)

pdf(paste0(fig_dir, "caremut_tme_state_atac_marker_genes.pdf"), width = 5, height = 5, useDingbats = FALSE, bg = "transparent")
heatmap_out
dev.off()

## Code to create a bubble plot from ArchR. Adapted from: https://github.com/NoemieL/bubble-plot-ArchR/blob/main/Script.
# This takes a bit of time to run if done interactively. Can skip as it does not provide a clear visual
gene_score <- getMatrixFromProject(CARE_filt_rna_all_filt, useMatrix="GeneScoreMatrix")
dense_matrix <- as.matrix(assay(gene_score))
gene_score_data <- as.data.frame(dense_matrix)
genedata <- rowData(gene_score)
gen_score_info <- cbind(genedata, gene_score_data)
gen_score_info2 <- as.data.frame(t(as.data.frame(gen_score_info[,7:length(colnames(gen_score_info))])))
colnames(gen_score_info2) <- genedata$name
rownames(gen_score_info2) <- colnames(gen_score_info)[7:length(colnames(gen_score_info))]
# gen_score_info2 is a large data.frame
for(i in rownames(gen_score_info2)){
  gen_score_info2[i,"CellType_final"]=CARE_filt_rna_all_filt$CellType_final[which(CARE_filt_rna_all_filt$cellNames==i)]
  gen_score_info2[i,"Samples"]=CARE_filt_rna_all_filt$Sample[which(CARE_filt_rna_all_filt$cellNames==i)]
}


# Genes for which to calculate the average accessibility score across the different cell types
gene_list <- c("PTPRZ1", "SOX6", "GPC5", "SLC1A2", "MBP", "ST18", "RBFOX1", "SYT1", "GAD1", "GRIK1", "ABCB1", "FLT1", "COL1A2", "DCN", "CD74", "DOCK8", "SKAP1", "TOX")

bubble_plot_info = data.frame()
for(i in gene_list){
  for(k in 1:length(unique(gen_score_info2$CellType_final))){
    a=nrow(bubble_plot_info)
    l=unique(gen_score_info2$CellType_final)[k]
    bubble_plot_info[a+1,"gene_name"]=i
    bubble_plot_info[a+1,"CellType_final"]=l
    eval(parse(text=(paste("bubble_plot_info[",a,"+1,'pct_exp']=(length(gen_score_info2[(gen_score_info2$",i,">0 & gen_score_info2$CellType_final=='",l,"'),'",i,"'])/nrow(gen_score_info2[gen_score_info2$CellType_final=='",l,"',]))*100", sep=""))))
    eval(parse(text=(paste("bubble_plot_info[",a,"+1,'avg_exp']=mean(gen_score_info2[gen_score_info2$",i,">0 & gen_score_info2$CellType_final=='",l,"','",i,"'])", sep=""))))
  }
}

# Save this information in case the plot needs to be remade/resized:
saveRDS(bubble_plot_info, file="/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac/archr_marker_gene_score_accessibility_all_celltypes.RDS")
library(scales)
library(RColorBrewer)

bubble_plot_info$CellType_final <- factor(bubble_plot_info$CellType_final, levels=c("Lymphocyte", "Myeloid", "Mural", "Endothelial", "InhNeuron", "ExcNeuron", "Oligodendrocyte", "Astrocyte", "Malignant"))
bubble_plot_info$gene_name <- factor(bubble_plot_info$gene_name, levels = gene_list)

pdf(paste0(fig_dir, "care_mut_rna_all_marker_gene_bubble_plot.pdf"), width=8, height=6, useDingbats = FALSE)
ggplot(data = bubble_plot_info, mapping = aes_string(x = 'gene_name', y = 'CellType_final')) +
  geom_point(mapping = aes_string(size = 'pct_exp', color = "avg_exp")) +
  scale_color_viridis_c(limits = c(0, 3.5), oob = scales::squish) +
  scale_size(range = c(0, 9), breaks = c(25, 50, 75), limits = c(0, 100)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(size = guide_legend(title = 'Percent\nopen chromatin'),
         color = guide_colorbar(title = "Average accessibility")) +
  labs(
    x = 'Gene',
    y = 'Cell type'
  )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
dev.off()

# Clean up these large objects.
rm(gene_score)
rm(gene_score_data)
rm(gen_score_info)
rm(gen_score_info2)
rm(dense_matrix)
gc()

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
## Make pseudo bulk measurements.
## There's a known issue with this on HPC systems: https://github.com/GreenleafLab/ArchR/issues/248
## Might be solved by setting threads = 1
# Pseudo-bulk refers to a grouping of single cells where the data from each single sample is combined into a single pseudo-sample.
set.seed(123)
CARE_filt_rna_all_filt <- addGroupCoverages(ArchRProj = CARE_filt_rna_all_filt, 
                                            groupBy = "CellType_final", 
                                            minCells = 100,
                                            maxCells = 1000,
                                            threads = getArchRThreads(),
                                            # Overwite the data in the ArchRProject object if the pseudo-bulk replicate information already exists
                                            force = TRUE)

## Is macs2 in the path variable?
pathToMacs2 <- findMacs2()

# Iterative overlap peak merging procedure
CARE_filt_rna_all_filt <- addReproduciblePeakSet(
  ArchRProj = CARE_filt_rna_all_filt, 
  groupBy = "CellType_final", 
  pathToMacs2 = pathToMacs2,
  threads = getArchRThreads(),
)

## Identifying marker peaks with ArchR.
CARE_filt_rna_all_filt <- addPeakMatrix(CARE_filt_rna_all_filt)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################

## Identifying marker PEAKS - features that are unique to a specific cell grouping.
##  ArchR takes into account for biases in data quality via TSSEnrichment and nFrags.
markersPeaks <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_all_filt, 
  useMatrix = "PeakMatrix", 
  groupBy = "CellType_final",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 2000
)


## Extract the marker peaks. Get access the GRanges object via `returnGR = TRUE`.
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")
# Tens of thousands of peaks for most cell types
lapply(markerList, nrow)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = TRUE,
  labelMarkers = NULL
)
heatmapPeaks
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-All-CellTypes", width = 8, height = 6, ArchRProj = CARE_filt_rna_all_filt, addDOC = FALSE)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
### Plot ATAC enrichment
CARE_filt_rna_all_filt <- addArchRAnnotations(ArchRProj = CARE_filt_rna_all_filt, collection = "ATAC")

enrichATAC <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CARE_filt_rna_all_filt,
  peakAnnotation = "ATAC",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)


heatmapATAC <- plotEnrichHeatmap(enrichATAC, 
                                 n = 10, 
                                 transpose = TRUE)
ComplexHeatmap::draw(heatmapATAC, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# The number in parentheses appears to be the max -log10(adj P value)
plotPDF(heatmapATAC, name = "ATAC-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_all_filt, addDOC = FALSE)

heatmapATAC_df <- plotEnrichHeatmap(enrichATAC, 
                                    n = 10, 
                                    transpose = TRUE,
                                    returnMatrix = TRUE)
library(viridisLite)
colnames(heatmapATAC_df) <- sapply(strsplit(colnames(heatmapATAC_df), " "), "[[", 1)
enrichment_of_interest <- c("Heme_CD8", "Heme_CD4", "IAtlas_T_MemoryTreg", "Heme_Mono", "Brain_Microglia", "IAtlas_Monocyte_Bulk",
                         "Cancer_HNSC", "Cancer_MESO", "Cancer_GBM", "Cancer_LGG",
                         "Brain_Opcs", "Brain_Opcs", "Brain_Excitatory_neurons", "Brain_Inhibitory_neurons", "Brain_Astrocytes", "Brain_Oligodendrocytes")
feature_order <- c("Lymphocyte", "Myeloid", "Mural", "Endothelial", "ExcNeuron", "InhNeuron","Astrocyte", "Malignant", "Oligodendrocyte")
heatmapATAC_df_filtered <- t(heatmapATAC_df[, colnames(heatmapATAC_df)%in%enrichment_of_interest])
heatmapATAC_df_ordered <- heatmapATAC_df_filtered[,feature_order]
manual_hmap <- ComplexHeatmap::Heatmap(heatmapATAC_df_ordered,                    
                                       show_row_dend = FALSE,
                                       cluster_columns = FALSE,
                                       show_column_dend = FALSE,
                                       col=viridis(100),
                                       name = "Norm. Enrichment -log10(P-adj) [0-Max]",
                                       heatmap_legend_param = list(
                                         legend_direction = "horizontal",
                                         legend_width = unit(5, "cm")
                                         ))

# Highly consistent enrichment across expected features for these cell types (malignant = TCGA GBM/LGG, Oligodendrocytes for normal oligodendrocytes etc.)
pdf(paste0(fig_dir, "all_celltypes_archr_differential_peak_enrichment_public_atac_heatmap.pdf"), width = 5, height = 5, useDingbats = FALSE, bg = "transparent")
draw(manual_hmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Add motif annotations.
if("Motif" %ni% names(CARE_filt_rna_all_filt@peakAnnotation)){
  CARE_filt_rna_all_filt <- addMotifAnnotations(ArchRProj = CARE_filt_rna_all_filt, motifSet = "cisbp", name = "Motif")
}

# Performing an enrichment of these motifs among the differentially accessible peaks.
enrichRegions <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = CARE_filt_rna_all_filt,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.05 & Log2FC >= 1"
)

heatmapRegions <- plotEnrichHeatmap(enrichRegions, 
                                    transpose = TRUE, 
                                    n = 5,
                                    clusterCols= FALSE)

plotPDF(heatmapRegions, name = "Regions-Enriched-Marker-Peak-Heatmap", width = 8, height = 6, ArchRProj = CARE_filt_rna_all_filt, addDOC = FALSE)

getAvailableMatrices(CARE_filt_rna_all_filt) # "GeneScoreMatrix" "PeakMatrix" "TileMatrix" 
names(CARE_filt_rna_all_filt@peakAnnotation) # "ATAC" "Motif"

# Save the project so that the annotations can be more easily accessed in the future. 
saveArchRProject(ArchRProj = CARE_filt_rna_all_filt, 
                 outputDirectory = "Save-CAREmut-All-RNA-Filtered", 
                 load = FALSE, 
                 dropCells = FALSE,
                 overwrite = TRUE) 

### END ###