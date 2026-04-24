##################################
# Assign malignant cells to CARE IDHmut metaprogram/state based on metaprogram scoring and p-value 
# Author: Kevin Johnson
##################################

library(tidyverse)
library(RColorBrewer)
library(viridis)
library(Seurat)
library(Matrix)
library(openxlsx)

proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
fig_dir     <- file.path(proj_dir, "figures/")
out_data_dir <- file.path(proj_dir, "processed_data/rna/")
script_dir  <- file.path(proj_dir, "scripts")

setwd(proj_dir)

# Avishay Spitzer provided these helper scripts. Use the score_within_samples*() functions for malignant cells.
source(file.path(script_dir, "utils", "plot_theme.R"))
source(file.path(script_dir, "utils", "caremut_utils.R"))

### ### ### ### ### ### ### ###
# Set-up
### ### ### ### ### ### ### ###
# Load in the CARE IDH-mutant metadata and count matrices.
md <- read.table(paste0(proj_dir, "/processed_data/rna/care_mut_cleaned_snrna_metadata_n75_20260320.txt"), sep = "\t", row.names = 1, header = TRUE)
md_trim <- md %>% 
  dplyr::select(SampleID, case_barcode, idh_codel_subtype, care_id) %>% 
  distinct()
rownames(md_trim) <- NULL

# Load in the CAREmut UMI produced at the beginning of the project.
umi_data_all <- readRDS("data/snrna/care_mut_umi_data_all_20230729.RDS")
names(umi_data_all)


# Prepare all relevant signatures in order to score CARE-MUT malignant cells

# Load CARE IDH-mutant metaprograms
mut_mp <- read.table("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/metaprograms/care_mut_selected_malignant_metaprograms.csv", sep = ",", header = T, row.names = 1)
mut_mp_list <- lapply(names(mut_mp), function(col_name) mut_mp[[col_name]])
names(mut_mp_list) <-paste0(colnames(mut_mp), "_MUT")

# Use top 50 genes from each signature, when possible, to be consistent with features used

# Public IDH-A and IDH-O signatures from Tirosh Nature 2016 and Venteicher 2017 Science papers (SmartSeq2 single cell data)
# Venteicher Astrocytoma
venteicher <- readWorkbook("data/misc/venteicher_table_s3.xlsx", startRow = 5, colNames = TRUE)
venteicher_signatures <- venteicher %>% 
  dplyr::select("Venteicher_OC_scRNA"=`Oligo-program.(Fig..2C)`, "Venteicher_AC_scRNA"=`Astro-program.(Fig..2C)`,
                "Venteicher_Stemness_scRNA"=`Stemness.program.(Fig..3C)`)
venteicher_signatures$Venteicher_OC_scRNA <- trimws(venteicher_signatures$Venteicher_OC_scRNA)
venteicher_signatures$Venteicher_AC_scRNA <- trimws(venteicher_signatures$Venteicher_AC_scRNA)
venteicher_signatures$Venteicher_Stemness_scRNA <- trimws(venteicher_signatures$Venteicher_Stemness_scRNA)
venteicher_sig_list <- as.list(venteicher_signatures[1:50,1:3])
venteicher_sig_list <- lapply(venteicher_sig_list, function(x) x[!is.na(x)])

# Tirosh Oligodendroglioma
tirosh <- readWorkbook("data/misc/tirosh_nature_2016_supplementary_table_1.xlsx", startRow = 9, colNames = TRUE)
tirosh_signatures <- tirosh %>% 
  dplyr::select("Tirosh_OC_scRNA"=`OC.(PCA-only)`,
                "Tirosh_AC_scRNA"=`AC.(PCA-only)`,
                "Tirosh_Stemness_scRNA"=`stemness`)
tirosh_signatures$Tirosh_OC_scRNA <- trimws(tirosh_signatures$Tirosh_OC_scRNA)
tirosh_signatures$Tirosh_AC_scRNA <- trimws(tirosh_signatures$Tirosh_AC_scRNA)
tirosh_signatures$Tirosh_Stemness_scRNA <- trimws(tirosh_signatures$Tirosh_Stemness_scRNA)
tirosh_signatures_list <- as.list(tirosh_signatures[1:50,1:3])

# Neftel et al Cell 2019 signatures 
neftel_signatures <- read.delim("data/misc/neftel_metamodule_genelists.csv", header = TRUE, sep = ",")
colnames(neftel_signatures) <- paste0(colnames(neftel_signatures), "_Neftel2019")
neftel_signatures_list <- as.list(neftel_signatures)

# CARE IDH-wildtype metaprograms
care_wt_mps <- read.delim("data/misc/carewt_malignant_metaprograms.txt", sep = "\t", header = T)
colnames(care_wt_mps) <- gsub("MP_", "WT_MP", colnames(care_wt_mps))
colnames(care_wt_mps) <- gsub("[0-9]_", "_", colnames(care_wt_mps))
colnames(care_wt_mps) <- gsub("MP1", "MP", colnames(care_wt_mps))
wt_mp_list <- lapply(names(care_wt_mps), function(col_name) care_wt_mps[[col_name]])
names(wt_mp_list) <-colnames(care_wt_mps)


# Public brain development and injury signatures
# Read in markers from Liu et al Cell 2023 - developing brain cells
liu_markers <- read.csv("data/misc/hNSPC_marker_genes.csv", header = T, stringsAsFactors = F) 

liu_markers <- apply(liu_markers, 2, function(x) {
  sapply(strsplit(x, split = "_"), function(y) y[[2]][1]) %>%
    head(50)
})

liu_markers_list <- setNames(lapply(1:ncol(liu_markers), function(i) liu_markers[, i]), colnames(liu_markers))
names(liu_markers_list) <- paste0(names(liu_markers_list), "_Liu2023")
names(liu_markers_list)[4] <- "GPC_Liu2023"
names(liu_markers_list)[10] <- "oRG_Liu2023"
names(liu_markers_list)[13] <- "cycling.RG_Liu2023"
names(liu_markers_list)[17] <- "vRG_Liu2023"

# Read in the data from Sadick et al Neuron 2022.
sadick_astrocytes <- readWorkbook("data/misc/sadick_neuron_2022_astrocyte_modules.xlsx", sheet = 2)
# Need to rename these Clusters based on function in the paper.
sadick_astrocytes_list <- lapply(sadick_astrocytes, function(x) head(x, 50))
names(sadick_astrocytes_list) <- c("C0_protective", "C1_reactive", "C2_protective", "C3_reactive", "C4_protective", "C5_protective",
                                   "C6_reactive", "C7_reactive", "C8_protective")
names(sadick_astrocytes_list) <- paste0(names(sadick_astrocytes_list), "-sadick2022")
reactive_astrocytes_list <- sadick_astrocytes_list[4]
names(reactive_astrocytes_list) <- "React.Astro_Sadick2022"


# Metaprograms from Gavish et al Nature 2023
gavish_signatures <- readWorkbook("data/misc/gavish_nature_2023_signatures.xlsx", sheet = 1)
colnames(gavish_signatures) <- paste0(colnames(gavish_signatures), "_Gavish2023")
gavish_signatures_select <- gavish_signatures[c(6, 12, 13, 16, 25, 26, 27, 28, 29, 38, 39)]
gavish_signatures_select_list <- lapply(gavish_signatures_select, function(x) head(x, 50))

sigs_list <- c(mut_mp_list, venteicher_sig_list, tirosh_signatures_list, neftel_signatures_list, wt_mp_list, liu_markers_list, reactive_astrocytes_list, gavish_signatures_select_list)

#########################################################################################
##### Assign within-sample metaprogram score
#########################################################################################

# Extract the malignant cells. I am removing SJ02-3, which had only a few malignant cells post-QC and should not be considered for downstream analyses.
md_malignant <- md %>% 
  filter(CellType_final=="Malignant", SampleID!="SJ02-3") 

# Define a function that subsets columns based on CellID variable
subset_columns <- function(mat, cell_ids) {
  col_idx <- which(colnames(mat) %in% cell_ids)
  mat[, col_idx, drop = FALSE]
}

# Apply the function to each element of the list
umi_data_all_malignant <- map(umi_data_all, subset_columns, cell_ids = md_malignant$CellID)
names(umi_data_all_malignant)

# Perform the scoring within a single sample. Note that AddModuleScore results differ depending on cell/sample dataset.
set.seed(43)
mp_scores <- score_within_samples_caremut(umi_data_all_malignant, md = md_malignant, sigs = sigs_list)

# Warning messages:
# 1: In simpleLoess(y, x, w, span, degree = degree, parametric = parametric,  ... :
#                    pseudoinverse used at -2.0523

# Save output so that it can be used for downstream analyses:
write.table(mp_scores, file = paste0(proj_dir, "/results/scoring/malignant_signature_scores_final_metaprograms.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# Reading backing in once complete.
# mp_scores <- read.table(paste0(proj_dir, "/results/scoring/malignant_signature_scores_final_metaprograms.txt"), sep = "\t", row.names = 1, header = TRUE)

mp_scores_md <- mp_scores %>% 
  dplyr::select(CellID, SampleID, MP_AC1_MUT:MP_CC_MUT)


#####################################################################
### Cell state scores - null distribution
#####################################################################
# mdata is an object (tibble) that contains the meta-data for the cells that should be classified (Sample, Patient, Timepoint etc.).
# The object also must include the CellID variable that identifies the cells.
# Due to the shuffling each variable in the tibble is theoretically normally distributed (or is at least close to being ND).
mdata <- mp_scores_md %>% 
  dplyr::select(CellID, SampleID, ends_with("MUT")) %>% 
  as_tibble()

## The approach requires variables to be named in a certain way:
sigs <- sigs_list[1:6]

umi_data_all_malignant_clean <- Filter(function(x) nnzero(x) > 0, umi_data_all_malignant)
names(umi_data_all_malignant_clean)

# We call this function to generate a NULL distribution to facilitate classification. According to the configured parameters
# it will sample 5000 cells from the pool of cells, shuffle the expression values while maintaining the mean expression of
# each gene and score the artificial cells for the meta-programs. It will repeat the process 20 times to generate a NULL 
# distribution of ~100K cells. It returns a tibble of ~100K x n (where n is the number of meta-programs)
set.seed(43)
permuted_data_all_mut <- generate_null_dist(umi_data_list = umi_data_all_malignant_clean,
                                    md = mdata,
                                    sigs = sigs,
                                    n_iter = 20, n_cells = 5000, verbose = T)

# Save output of permuted data
saveRDS(permuted_data_all_mut, paste0(out_data_dir, "caremut_signatures_permuted_data_all.RDS"))
#permuted_data_all_mut <- readRDS(paste0(out_data_dir, "caremut_signatures_permuted_data_all.RDS"))

permuted_data <- permuted_data_all_mut

state_programs <- sigs

vars <- names(state_programs)

# We fit a normal distribution to each of the variables of the permuted data
scores_nd <- lapply(colnames(permuted_data), function(mp) {
  
  x <- permuted_data %>% pull(mp)
  
  fit <- MASS::fitdistr(x, "normal")
  class(fit)
  
  para <- fit$estimate
  
  tibble(MP = mp, Mean = para[1], SD = para[2])  
})

# This provides the mean and standard deviation per signature/metaprogram
scores_nd <- do.call(rbind, scores_nd)
# Create vectors for these data
mean_vec <- setNames(scores_nd$Mean, scores_nd$MP)
sd_vec <- setNames(scores_nd$SD, scores_nd$MP)

####################################################################################################################################
# Plot the actual scores vs. the NULL distribution for each MP
####################################################################################################################################
library(reshape2)
library(ggdist)
set.seed(43)
norm_fit <- lapply(colnames(permuted_data), function(mp) tibble(MP = mp,
                                                                Sig = rnorm(n = 100000,
                                                                            mean = mean_vec[mp],
                                                                            sd = sd_vec[mp])))
norm_fit <- do.call(rbind, norm_fit)

dm <- rbind(melt(data = permuted_data,
                 measure.vars = vars) %>% mutate(DataType = "Permuted"),
            melt(data = mdata %>% select(ends_with("_MUT")),
                 measure.vars = vars) %>% mutate(DataType = "Actual"),
            norm_fit %>%
              rename(variable = MP, value = Sig) %>%
              mutate(DataType = "Classifier")) %>%
  mutate(DataType = factor(DataType, c("Permuted", "Actual", "Classifier")))

dm_stats <- dm %>%
  group_by(variable) %>%
  filter(DataType == "Classifier") %>%
  summarise(Q95 = quantile(value, .95), Q99 = quantile(value, .99))

ggplot(data = dm, aes(x = value, y = after_stat(ncount), color = DataType, linetype = DataType)) +
  facet_wrap(~variable, scales = "free_x", nrow = 2) +
  geom_freqpoly(bins = 100, size = 1, show.legend = c(color = T, linetype = F)) +
  scale_color_manual(name = "", values = c("Permuted" = "dodgerblue", "Actual" = "red", "Classifier" = "black")) +
  scale_linetype_manual(values = c("Permuted" = "solid", "Actual" = "solid", "Classifier" = "dashed")) +
  scale_fill_discrete(name = "Data distribution") +
  geom_vline(data = dm_stats, mapping = aes(xintercept = Q99), linetype = "dashed", size = 1) +
  xlab("Program score") +
  ylab("Count (scaled to 1)") + 
  plot_theme


####################################################################################################################################
# Cell state classification
####################################################################################################################################
# Names of all MPs to-be-classified should be included in this vector (can be used to exclude MPs that reflect artifact/low quality etc.)
state_vars <- names(state_programs)

# Melt the data to make the computation easier
state_data <- melt(data = mdata %>%
                     select(CellID, all_of(state_vars)),
                   id.vars = "CellID",
                   variable.name = "Program",
                   value.name = "Score",
                   measure.vars = state_vars)
state_data <- as_tibble(state_data)
state_data$Program <- as.character(state_data$Program)
table(state_data$Program)

# Generate a Z-score for each MP (using the mean and SD of the NULL distribution)
state_data <- state_data %>%
  mutate(Score_z = (Score - mean_vec[Program]) / sd_vec[Program])

# Compute a p-value for each (CellID, MP) pair. This is a one-sided test with the hypothesis that
# the actual score is not greater than expected by chance
state_data <- state_data %>%
  mutate(p.val = pnorm(Score_z, lower.tail = F))

# Correct for multiple testing **within each cell** using the Holm method
state_data <- state_data %>%
  group_by(CellID) %>%
  mutate(p.adj = p.adjust(p.val, "holm"))

# We consider the p-value to be significant if the adjusted p-value is less than 0.05. 
# However, we can also use nominal p-value depending on threshold/dataset
state_data$p.sig <- state_data$p.adj < .05

# Compute the classification statistics for each MP/gene set 
state_stats <- state_data %>%
  group_by(Program) %>%
  summarise(n = sum(p.sig == T), N = n(), Freq = n / N)


# This is the actual classification step. We filter out the statistically insignificant scores and classify the cell
# to the MP with the maximal signal. Classify CC separate from the other states.
# The classification will do the following:
# 1. Separate out cell cycle since cells can be classified as cycling independent of AC/MES/OPC/NPC state
# 2. Restrict only those programs with adj. p-value < 0.05
# 3. Rank by descending score such multiple significant scores will be assigned to the one with the highest Score (AddModuleScore).
state_data_classify <- state_data %>% 
  filter(!Program%in%c("MP_CC_MUT")) %>% 
  group_by(CellID) %>%
  filter(p.sig == T) %>%
  arrange(desc(Score)) %>%
  filter(!duplicated(CellID)) %>%
  ungroup()

# Perform this step separately for Cell Cycle:
state_data_classify_cc <- state_data %>% 
  filter(Program=="MP_CC_MUT") %>% 
  group_by(CellID) %>%
  filter(p.sig == T) %>%
  arrange(desc(Score)) %>%
  filter(!duplicated(CellID)) %>%
  ungroup()


####################################################################################################################################
# Assign the state
####################################################################################################################################
# If you have a large proportion of cells that are "Undifferentiated" (i.e. they did not achieve a significant adjusted p-value for any MP)
# then you can adjust the multiple-adjustment method or classification threshold to less stringent values.
state_vec <- setNames(rep("Undifferentiated", nrow(mdata)), mdata$CellID)
state_vec[state_data_classify$CellID] <- state_data_classify$Program
table(state_vec)
table(state_vec) / length(state_vec)

mdata$State <- state_vec[mdata$CellID]
table(mdata$State)

state_complexity <- mdata %>% 
  inner_join(md, by=c("CellID", "SampleID")) 

state_complexity %>% 
  mutate(State = factor(State, c(names(mut_mp_list), "Undifferentiated"))) %>%
  ggplot(aes(x = State, y = nFeature_RNA)) + 
  ggdist::stat_halfeye(adjust = .5, width = .75, justification = -.2, .width = 0, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.color = NA) +
  coord_cartesian(xlim = c(1.2, NA)) +
  xlab("") +
  ylab("Complexity") +
  geom_hline(yintercept = median(state_complexity$nFeature_RNA), linetype = "dashed", size = 1, color = "red") +
  scale_y_continuous(breaks = seq(0, 7000, 1000)) +
  theme_bw() +
  theme(panel.grid.major = element_line())

## 1. Cell cycle is not considered a cellular state but rather a feature (since cells can have a clear identity such as OPC or NPC and still be cycling) 
mdata <- mdata %>% 
  inner_join(md_trim, by="SampleID") %>% 
  mutate(isCC = ifelse(mdata$CellID%in%state_data_classify_cc$CellID, TRUE, FALSE))

table(mdata$State, mdata$isCC)

# Write out the current cell state classification approach.
write.table(mdata, file = paste0(proj_dir, "/results/scoring/caremut_malignant_cell_state_assignment.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)


### END ###