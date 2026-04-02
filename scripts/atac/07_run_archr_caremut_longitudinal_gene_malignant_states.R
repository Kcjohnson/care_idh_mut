##############################
### Run ArchR differentially peak activity analysis on CAREmut multiome ATAC data longitudinal analysis of malignant cell states
### Author: Kevin Johnson
### Updated: 2024.06.08
##############################

## Part 7: Longitudinal per-patient analysis of differential gene chromatin accessibility per cell state

## ArchR creates several directories automatically when creating arrow files and ArchR projects.
workdir <- "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/"
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
fig_dir <- "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/figures/archr/"

# Load sample-specific metadata
sample_md <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/data/metadata/clinical_samples_genomic_md_20240608.txt", sep="\t", header = TRUE)


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
mut_md <- read.table("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20240212.txt", sep = "\t", header = TRUE)

# For this analysis, only ATAC profiles were generated for samples processed by the Verhaak lab and defined by RNA as malignant cells
mut_md_verhaak <- mut_md %>% 
  filter(lab=="Verhaak lab", CellType_final=="Malignant")

## Create a data.frame that contains the essential information to be merged and inspected.
atac_df <- data.frame(projMONITOR_filt$cellNames, projMONITOR_filt@cellColData$Sample, projMONITOR_filt@cellColData$TSSEnrichment)

# The cell names are a little different between RNA (Seurat) and ATAC (ArchR). Need to create a common linker.
atac_df$CellID <- gsub("#", "-", atac_df$projMONITOR_filt.cellNames)
# 71,088 cells present in ATACseq data that pass doublet removal and are also found in snRNAseq data that passed QC for IDH-mutant.
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


all(CARE_filt_rna_malignant$CellID==cell_class_filt_ord$CellID)
CARE_filt_rna_malignant$care_id <- cell_class_filt_ord$care_id
CARE_filt_rna_malignant$idh_codel_subtype <- cell_class_filt_ord$idh_codel_subtype


cell_class_filt_ord_annot <- cell_class_filt_ord %>% 
  left_join(sample_md, by=c("sample_barcode", "case_barcode", "idh_codel_subtype", "care_id", "patient_id", "timepoint"))
all(CARE_filt_rna_malignant$CellID==cell_class_filt_ord_annot$CellID)
CARE_filt_rna_malignant$grade <- paste0("G", cell_class_filt_ord_annot$grade_num)


# CRITICAL STEP!
# See following github discussion for treatment vs. control style anylases: https://github.com/GreenleafLab/ArchR/discussions/696
CARE_filt_rna_malignant <- addCellColData(ArchRProj = CARE_filt_rna_malignant, data = paste0(CARE_filt_rna_malignant@cellColData$CellStateGroup,"_x_",CARE_filt_rna_malignant@cellColData$Sample), name = "StateBySample", cells = getCellNames(CARE_filt_rna_malignant), force = TRUE)
table(CARE_filt_rna_malignant$StateBySample)

CARE_filt_rna_malignant <- addCellColData(ArchRProj = CARE_filt_rna_malignant, data = paste0(CARE_filt_rna_malignant@cellColData$CellStateGroup,"_x_",CARE_filt_rna_malignant@cellColData$idh_codel_subtype), name = "StateBySubtype", cells = getCellNames(CARE_filt_rna_malignant), force = TRUE)
table(CARE_filt_rna_malignant$StateBySubtype)

CARE_filt_rna_malignant <- addCellColData(ArchRProj = CARE_filt_rna_malignant, data = paste0(CARE_filt_rna_malignant@cellColData$CellStateGroup,"_x_",CARE_filt_rna_malignant@cellColData$idh_codel_subtype, "_x_",CARE_filt_rna_malignant@cellColData$grade), name = "StateByGrade", cells = getCellNames(CARE_filt_rna_malignant), force = TRUE)
table(CARE_filt_rna_malignant$StateByGrade)

# Compare across subtypes. Note that the maximum cells allowed to be used for a single group are 500 by default. Let's increase that for some of these analyses.

markerState <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "StateBySubtype",
  testMethod = "wilcoxon",
  maxCells = 5000,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "AC_x_IDHmut-noncodel", 
  bgdGroups = "AC_x_IDHmut-codel"
)

upregulated <- getMarkers(markerState, cutOff = "FDR <= 0.05 & Log2FC >= 1")
downregulated <- getMarkers(markerState, cutOff = "FDR <= 0.05 & Log2FC <= -1")
upregulated_genes <- as.data.frame(upregulated$`AC_x_IDHmut-noncodel`)
downregulated_genes <- as.data.frame(downregulated$`AC_x_IDHmut-noncodel`)
dag_genes_ac <- bind_rows(upregulated_genes, downregulated_genes) %>% 
  mutate(state = "AC-like")

markerState_opc <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "StateBySubtype",
  testMethod = "wilcoxon",
  maxCells = 5000,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "OPC_x_IDHmut-noncodel", 
  bgdGroups = "OPC_x_IDHmut-codel"
)

upregulated_opc <- getMarkers(markerState_opc, cutOff = "FDR <= 0.05 & Log2FC >= 1")
downregulated_opc <- getMarkers(markerState_opc, cutOff = "FDR <= 0.05 & Log2FC <= -1")
upregulated_genes_opc <- as.data.frame(upregulated_opc$`OPC_x_IDHmut-noncodel`)
downregulated_genes_opc <- as.data.frame(downregulated_opc$`OPC_x_IDHmut-noncodel`)
dag_genes_opc <- bind_rows(upregulated_genes_opc, downregulated_genes_opc) %>% 
  mutate(state = "OPC-like")

markerState_npc <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "StateBySubtype",
  testMethod = "wilcoxon",
  maxCells = 5000,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "NPC_x_IDHmut-noncodel", 
  bgdGroups = "NPC_x_IDHmut-codel"
)

upregulated_npc <- getMarkers(markerState_npc, cutOff = "FDR <= 0.05 & Log2FC >= 1")
downregulated_npc <- getMarkers(markerState_npc, cutOff = "FDR <= 0.05 & Log2FC <= -1")
upregulated_genes_npc <- as.data.frame(upregulated_npc$`NPC_x_IDHmut-noncodel`)
downregulated_genes_npc <- as.data.frame(downregulated_npc$`NPC_x_IDHmut-noncodel`)
dag_genes_npc <- bind_rows(upregulated_genes_npc, downregulated_genes_npc) %>% 
  mutate(state = "NPC-like")

markerState_mes <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "StateBySubtype",
  testMethod = "wilcoxon",
  maxCells = 5000,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "MES_x_IDHmut-noncodel", 
  bgdGroups = "MES_x_IDHmut-codel"
)

upregulated_mes <- getMarkers(markerState_mes, cutOff = "FDR <= 0.05 & Log2FC >= 1")
downregulated_mes <- getMarkers(markerState_mes, cutOff = "FDR <= 0.05 & Log2FC <= -1")
upregulated_genes_mes <- as.data.frame(upregulated_mes$`MES_x_IDHmut-noncodel`)
downregulated_genes_mes <- as.data.frame(downregulated_mes$`MES_x_IDHmut-noncodel`)
dag_genes_mes <- bind_rows(upregulated_genes_mes, downregulated_genes_mes) %>% 
  mutate(state = "MES-like")

markerState_undiff <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "StateBySubtype",
  testMethod = "wilcoxon",
  maxCells = 5000,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Undifferentiated_x_IDHmut-noncodel", 
  bgdGroups = "Undifferentiated_x_IDHmut-codel"
)

upregulated_undiff <- getMarkers(markerState_undiff, cutOff = "FDR <= 0.05 & Log2FC >= 1")
downregulated_undiff <- getMarkers(markerState_undiff, cutOff = "FDR <= 0.05 & Log2FC <= -1")
upregulated_genes_undiff <- as.data.frame(upregulated_undiff$`Undifferentiated_x_IDHmut-noncodel`)
downregulated_genes_undiff <- as.data.frame(downregulated_undiff$`Undifferentiated_x_IDHmut-noncodel`)
dag_genes_undiff <- bind_rows(upregulated_genes_undiff, downregulated_genes_undiff) %>% 
  mutate(state = "Undifferentiated")


state_dag_genes <- bind_rows(dag_genes_ac,
                             dag_genes_opc,
                             dag_genes_npc,
                             dag_genes_mes,
                             dag_genes_undiff) %>% 
  mutate(direction = ifelse(Log2FC>=1, "IDH-A upregulated", "IDH-A downregulated"))

write.table(state_dag_genes, "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/state_controlled_dag/state_controlled_subtype_dag.txt", sep="\t", row.names = FALSE, col.names = TRUE)

top_genes <- state_dag_genes %>% 
  dplyr::group_by(name, seqnames) %>% 
  dplyr::summarise(gene_hits = n()) %>% 
  filter(gene_hits > 4)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# ATAC AddModuleScores for RNA-based metaprograms
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## Define markers based on cell state. Use 2000 cells as max to keep groups equal sized
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
lapply(markerListGS, nrow)
undiff_atac <-markerListGS$Undifferentiated
undiff_atac_sig_hits <- undiff_atac %>% 
  as.data.frame() %>% 
  arrange(Log2FC)
undiff_genes <- list(rev(undiff_atac_sig_hits$name))
names(undiff_genes) <- "undifferentiated_atac"
undiff_genes_top50 <- undiff_genes[[1:50]]

mut_mp <- read.table("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/nmf_res/metaprograms/malignant_selected/care_mut_selected_malignant_metaprograms.csv", sep = ",", header = T, row.names = 1)
mut_mp_list <- lapply(names(mut_mp), function(col_name) mut_mp[[col_name]])
names(mut_mp_list) <-paste0(colnames(mut_mp), "_MUT")

# Need to filter to those genes present in the dataset.
all_genes <- getFeatures(CARE_filt_rna_malignant)

filtered_mp_list <- lapply(mut_mp_list, function(genes) {
  genes[genes %in% all_genes]
})

CARE_filt_rna_malignant <- addModuleScore(CARE_filt_rna_malignant, features = filtered_mp_list, useMatrix = "GeneScoreMatrix")
CARE_filt_rna_malignant <- addModuleScore(CARE_filt_rna_malignant, features = undiff_genes, useMatrix = "GeneScoreMatrix")
CARE_filt_rna_malignant <- addModuleScore(CARE_filt_rna_malignant, features = undiff_genes_top50, useMatrix = "GeneScoreMatrix")

df <- getCellColData(CARE_filt_rna_malignant, select = c("CellID", "CellStateGroup",  "Module.MP_AC1_MUT", "Module.MP_OPC_MUT", "Module.MP_NPC_MUT", "Module.MP_CC_MUT", "Module.undifferentiated_atac"))
ggplot(df %>% data.frame(), aes(x=CellStateGroup, y=Module.MP_OPC_MUT)) +
  geom_boxplot()
ggplot(df %>% data.frame(), aes(x=CellStateGroup, y=Module.MP_AC1_MUT)) +
  geom_boxplot()
ggplot(df %>% data.frame(), aes(x=CellStateGroup, y=Module.MP_NPC_MUT)) +
  geom_boxplot()
ggplot(df %>% data.frame(), aes(x=CellStateGroup, y=Module.undifferentiated_atac)) +
  geom_boxplot()

df_out <- df %>% data.frame()
write.table(df_out, "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/atac_malignant_module_scores.txt", sep="\t", row.names = FALSE, col.names = TRUE)

## Examine some of the differences in RNA metaprogram accessibility across time points
metaprogram_atac <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/atac_malignant_module_scores.txt", header = TRUE)
patient_md <- read.delim("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/clinical_patient_genomic_md_20240608.txt", sep="\t", header = TRUE)


metaprogram_atac <- metaprogram_atac %>% 
  mutate(timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2)))

cell_scores <- metaprogram_atac %>% 
  dplyr::select(Module.MP_AC1_MUT:Module.MP_NPC_MUT)

# Simple metric to define RNA metaprogram accessibility differences - Standard deviation 
metaprogram_atac$module_standard_dev <- apply(cell_scores, 1, function(x) {
  sd(x)   # higher SD = more committed
})

metaprogram_atac$entropy_based <- apply(cell_scores, 1, function(x) {
  x_prob <- x - min(x)          # shift so no negative values
  if(sum(x_prob) == 0) x_prob <- rep(1, length(x)) # handle zero vector
  x_prob <- x_prob / sum(x_prob) # normalize to probabilities
  -sum(x_prob * log2(x_prob + 1e-9)) # Shannon entropy
})

metaprogram_atac$entropy_scaled <- (metaprogram_atac$entropy_based - min(metaprogram_atac$entropy_based)) /
  (max(metaprogram_atac$entropy_based) - min(metaprogram_atac$entropy_based))

metaprogram_summary <- metaprogram_atac %>% 
  dplyr::group_by(care_id, CellStateGroup) %>% 
  dplyr::summarise(mean_module_sd = mean(module_standard_dev),
                   mean_entropy = mean(entropy_scaled),
                   cells = n(), 
                   avg_FRIP = mean(FRIP),
                   avg_TSSEnrichment = mean(TSSEnrichment),
                   avg_nFrags = mean(nFrags)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(timepoint = paste0("T", sapply(strsplit(care_id, "T"), "[[", 2))) %>% 
  dplyr::mutate(patient = sapply(strsplit(care_id, "T"), "[[", 1))

cases_to_keep <- metaprogram_atac %>% 
  dplyr::select(care_id, timepoint) %>% 
  distinct() %>% 
  filter(timepoint!="T3") %>% 
  dplyr::mutate(patient = sapply(strsplit(care_id, "T"), "[[", 1)) %>% 
  dplyr::group_by(patient) %>%
  dplyr::summarise(counts = n()) %>% 
  filter(counts > 1)

sample_states_to_keep <- metaprogram_atac %>% 
  filter(timepoint!="T3") %>% 
  mutate(care_id_state = paste0(care_id, "-", CellStateGroup)) %>% 
  dplyr::group_by(care_id_state) %>%
  dplyr::summarise(counts = n()) %>% 
  filter(counts > 25) 

metaprogram_longitudinal_df <- metaprogram_summary %>% 
  dplyr::mutate(care_id_state = paste0(care_id, "-", CellStateGroup)) %>% 
  filter(timepoint!="T3", patient%in%cases_to_keep$patient, care_id_state%in%sample_states_to_keep$care_id_state) %>% 
  inner_join(patient_md, by=c("patient"="patient_id"))
  
source("/vast/palmer/pi/verhaak/kcj28/care_mut/scripts/misc/plot_theme.R")

pdf("/vast/palmer/pi/verhaak/kcj28/care_mut/results/figures/archr/atac_module_score_variability.pdf", height = 4, width = 5)
ggplot(metaprogram_longitudinal_df %>% 
         mutate(timepoint = recode(timepoint, `T1` = "Sample 1",
                                   `T2` = "Sample 2")), aes(x=timepoint, y=mean_module_sd)) +
  geom_boxplot() +
  geom_point() +
  #geom_line(aes(group=patient), color="gray70", linetype=2) +
  plot_theme +
  labs(x = "Time point", y = "Mean RNA metaprogram\naccessility score standard deviation") +
  facet_grid(.~CellStateGroup) +
  stat_compare_means(method="wilcox", label = "p.format") +
  theme(axis.text.x  = element_text(angle = 45, hjust=1))
dev.off()

ggplot(metaprogram_longitudinal_df, aes(x=timepoint, y=mean_module_sd)) +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group=patient), color="gray70", linetype=2) +
  plot_theme +
  labs(x = "Sample time point", y = "Mean RNA metaprogram\naccessility score standard deviation") +
  facet_grid(acquired_genetic_alt_t1t2~CellStateGroup) 

# Entropy appears to decrease.
ggplot(metaprogram_longitudinal_df, aes(x=timepoint, y=mean_entropy)) +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group=patient), color="gray70", linetype=2) +
  plot_theme +
  labs(x = "Time point", y = "Mean metaprogram entropy") +
  facet_grid(.~CellStateGroup) 

ggplot(metaprogram_longitudinal_df, aes(x=timepoint, y=avg_FRIP)) +
  geom_boxplot() +
  geom_point() +
  geom_line(aes(group=patient), color="gray70", linetype=2) +
  geom_line() + 
  plot_theme  +
  facet_grid(.~CellStateGroup)+ 
  labs(x = "Time point", y = "Mean Fraction of Reads in Promoters")
  

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#Patient specific longitudinal changes
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

## Objective: For each patient and cell state loop through longitudinal comparisons so
# that differential peaks can be determined within patient while controlling for cell state.

# Set up for loop run
archr_project_name <- "CARE_filt_rna_malignant"
patient_id <- unique(sapply(strsplit(CARE_filt_rna_malignant@cellColData$Sample, "-"), "[[", 1))
cell_state = unique(CARE_filt_rna_malignant@cellColData$CellStateGroup)

# Create empty vectors and a list to store the results.
up_genes <- down_genes <- use_cells <- bgd_cells <- c()
outlist <- list()

atac_md_filt <- as.data.frame(CARE_filt_rna_malignant@cellColData)
table(atac_md_filt$CellStateGroup, atac_md_filt$Sample)

# Check to see whether there will be any cases where there are 3 tumor samples per patient. Manually remove any 3 time point samples with low numbers 
table(CARE_filt_rna_malignant@cellColData$Sample)

# Determine the longitudinal difference in cell states. Since we are interested in more subtle changes, let's use a lower threshold for significant difference.
# Log2FC of 1 would represent approximately ~100% increase in accessibility.

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
    if(use_cells < 50) next
    if(bgd_cells < 50) next
    
    out <- paste0(sprintf("Analyzing patient %s for the comparison: ", patient_id[i]),  sprintf("Use group: %s", state_use), " vs ", sprintf("Background group: %s", state_bgd))
    print(out)
    
    start.time <- Sys.time()
    # Note that the maximum cells allowed to be used for a single group are 500 by default. This is okay for these analyses because there are some samples 
    # with many more cells at a given time point.
    markerState <- getMarkerFeatures(
      ArchRProj = CARE_filt_rna_malignant, 
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
    write.table(dag_df, paste0("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/state_controlled_dag/", comparison, "_longitudinal_dag.txt"), sep="\t", row.names = FALSE, col.names = TRUE)
    
    results <- data.frame(patient_id[i], cell_state[j], state_use, use_cells, state_bgd, bgd_cells, up_genes, down_genes)
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
  mutate(SampleID = sapply(strsplit(state_use, "_x_"), "[[", 2)) %>% 
  inner_join(sample_filt, by="SampleID") 

patient_list_results_filt 

# Write out the results so that they can be used downstream. 
write.table(patient_list_results_filt, "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/archr_care_cell_type_longitudinal_gene_score_min40cells_perstate_20240608.txt", sep="\t", row.names = FALSE, col.names = TRUE)
patient_list_results_filt <- read.delim("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/archr_care_cell_type_longitudinal_gene_score_min40cells_perstate_20240608.txt", sep="\t", header = TRUE)

res <- patient_list_results_filt %>% 
  left_join(patient_md, by="patient_id")
res$down_genes <- ifelse(res$down_genes!=0, res$down_genes*-1, res$down_genes)

library(tidyverse)
library(ggpubr)
library(EnvStats)
source("/gpfs/gibbs/pi/verhaak/kcj28/care_mut/scripts/misc/plot_theme.R")

malignant_res <- res %>% 
  dplyr::select(patient_id = patient_id.i., cell_state = cell_state.j., grade_change_t1t2, acquired_genetic_alt_t1t2, idh_codel_subtype, up_genes, down_genes) %>% 
  pivot_longer(cols= c(up_genes:down_genes),
               names_to = "Type",
               values_to = "Values") 
malignant_res$acquired_genetic_alt_t1t2 <- ifelse(is.na(malignant_res$acquired_genetic_alt_t1t2), "No acquired alt.", "Acquired genetic alt.")
malignant_res$Type <- factor(malignant_res$Type, levels=c("up_genes", "down_genes"))
malignant_res$Type <- ifelse(malignant_res$Type=="up_genes", "Increased accessibility", "Decreased accessibility")


png(paste0(fig_dir, "longitudinal_state_controll_atac_malignant_gene_scores_codel.png"), width = 9, height = 5, units = 'in', res = 300, bg = "transparent")
ggplot(malignant_res %>% filter(cell_state%in%c("AC", "OPC"), idh_codel_subtype=="IDHmut-codel"), aes(x = patient_id, y = Values, fill = Type)) +
  geom_bar(stat = 'identity', position = 'identity', color = 'black') +
  labs(x = 'Patients', y = 'Longitudinal differential accessibility\nstate-conrolled (IDH-O)', fill="Longitudinal change") +
  scale_fill_manual(values = c('blue','red')) +
  plot_theme +
  facet_grid(cell_state~acquired_genetic_alt_t1t2, scales="free") +
  ylim(-4500,2000) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

patients_to_include <- patient_list_results_filt %>% 
  filter(use_cells>100, bgd_cells>100)

ggplot(malignant_res %>% filter(cell_state%in%c("AC", "OPC", "Undifferentiated"),
                                patient_id%in%patients_to_include$patient_id.i.), aes(x = patient_id, y = Values, fill = Type)) +
  geom_bar(stat = 'identity', position = 'identity', color = 'black') +
  labs(x = 'Patients', y = 'Longitudinal differential accessibility\nstate-conrolled (IDH-O)', fill="Longitudinal change") +
  scale_fill_manual(values = c('blue', 'red')) +
  plot_theme +
  facet_grid(cell_state~acquired_genetic_alt_t1t2, scales="free") +
  ylim(-4500,2000) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 

ggplot(malignant_res %>% filter(cell_state%in%c("AC", "OPC", "Undifferentiated"), patient_id%in%c("SJ04","SJ06", "SJ07", "SJ13", "SJ15")), aes(x = patient_id, y = Values, fill = Type)) +
  geom_bar(stat = 'identity', position = 'identity', color = 'black') +
  labs(x = 'Patients', y = 'Up or down longitudinal ATAC peaks', fill="Longitudinal change") +
  scale_fill_manual(values = c('blue', 'red')) +
  plot_theme +
  facet_grid(cell_state~., scales="free") 

# Create a facet box plot for each cell state comparing the total number of differentially accessible genes.
patients_to_include <- patient_list_results_filt %>% 
  filter(use_cells>20, bgd_cells>20)

malignant_res_compare <- malignant_res %>% filter(cell_state%in%c("AC", "OPC", "Undifferentiated"),
                         patient_id%in%patients_to_include$patient_id.i.) %>% 
  mutate(all_dag = ifelse(Type=="Decreased accessibility", -1*Values, Values)) %>% 
  group_by(patient_id, cell_state, acquired_genetic_alt_t1t2) %>% 
  summarise(total_dag = sum(all_dag))

ggplot(malignant_res_compare, aes(x=acquired_genetic_alt_t1t2, y=total_dag)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color=cell_state)) +
  stat_compare_means(method = "wilcox", paired = FALSE, size =4, label="p.format") +
  #facet_grid(.~cell_state, scales="free") +
  plot_theme + 
  stat_n_text()

ggplot(malignant_res, aes(x = patient_id, y = Values, fill = Type)) +
  geom_bar(stat = 'identity', position = 'identity', color = 'black') +
  labs(title = 'Longitudinal differentially accessible peaks (malignant)',
       x = 'Patients', y = 'Up or down longitudinal ATAC peaks', fill="Longitudinal change") +
  scale_fill_manual(values = c('blue', 'red')) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(.~grade_change, scales="free", space="free") 
dev.off()

## Testing whether all T1 vs all T2 could be used.
CARE_filt_rna_malignant$timepoint <- paste0("T", sapply(strsplit(as.character(CARE_filt_rna_malignant@cellColData$care_id), "T"), "[[", 2))
CARE_filt_rna_malignant$patient_id <- sapply(strsplit(as.character(CARE_filt_rna_malignant@cellColData$care_id), "T"), "[[", 1)
CARE_filt_rna_malignant$state_timepoint <- paste0(CARE_filt_rna_malignant$CellStateGroup, "_", CARE_filt_rna_malignant$timepoint )

markerState_all <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "state_timepoint",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "OPC_T2", 
  bgdGroups = "OPC_T1"
)

markerList_t1t2 <- getMarkers(markerState_all, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1")
dag_genes_opc <-  as.data.frame(markerList_t1t2[["OPC_T2"]])

markerState_all_ac <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "state_timepoint",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "AC_T2", 
  bgdGroups = "AC_T1"
)

markerList_t1t2_ac <- getMarkers(markerState_all_ac, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1")
dag_genes_ac <-  as.data.frame(markerList_t1t2_ac[["AC_T2"]])

markerState_all_undiff <- getMarkerFeatures(
  ArchRProj = CARE_filt_rna_malignant, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "state_timepoint",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Undifferentiated_T2", 
  bgdGroups = "Undifferentiated_T1"
)

markerList_t1t2_undiff <- getMarkers(markerState_all_undiff, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1")
dag_genes_undiff <-  as.data.frame(markerList_t1t2_undiff[["Undifferentiated_T2"]])

#
# Get a list of files with the specified suffix
file_list <- list.files(path = "/gpfs/gibbs/pi/verhaak/kcj28/care_mut/results/archr/state_controlled_dag/", pattern = "_longitudinal_dag.txt$", full.names = TRUE)

# Initialize an empty list to store data frames
data_frames <- list()

# Loop through each file and process
for (file in file_list) {
  # Read the file into a data frame
  df <- read_delim(file, delim = "\t", col_types = cols())
  
  # Check if the data frame has only a header (i.e., no rows)
  if (nrow(df) == 0) {
    warning(paste("File", file, "has only a header and will be skipped."))
    next
  }
  
  # Extract the file name without the suffix
  file_name <- basename(file)
  file_name <- sub("_longitudinal_dag.txt$", "", file_name)
  
  # Add the file name as a new column in the data frame
  df <- df %>%
    mutate(file_name = file_name)
  
  # Append the data frame to the list
  data_frames[[file_name]] <- df
}

data_frames_w_hits <- data_frames[1:41]
all_dag_longitudinal <- bind_rows(data_frames_w_hits) %>% 
  mutate(case = sapply(strsplit(file_name, "_"), "[[", 1),
         state = sapply(strsplit(file_name, "_"), "[[", 2))

all_dag_longitudinal_select <- all_dag_longitudinal %>% 
  filter(case=="SJ04", Log2FC>0)

upregulated_freq <- all_dag_longitudinal_select %>%
  as.data.frame() %>% 
  dplyr::group_by(name) %>% 
  dplyr::summarise(obs = n()) %>% 
  filter(obs > 2)

all_dag_longitudinal_select <- all_dag_longitudinal %>% 
  filter(case=="SJ08", Log2FC>0)

upregulated_freq <- all_dag_longitudinal %>%
  as.data.frame() %>% 
  dplyr::group_by(name) %>% 
  dplyr::summarise(obs = n()) %>% 
  filter(obs >= 10)


### END ###