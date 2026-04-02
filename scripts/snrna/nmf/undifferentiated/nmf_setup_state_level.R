##################################
# Preparing the snRNAseq data from CAREmut dataset for NMF in the undifferentiated malignant population
# Author: Kevin Johnson
# Date: 2026.04.01
##################################

library(tidyverse)
library(Seurat)
library(Matrix)


## Specify directories:
proj_dir    <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut"
out_data_dir <- file.path(proj_dir, "results/nmf_res_caremut/undifferentiated/")

setwd(proj_dir)

##############################################################
#### Load dataset
##############################################################
umi_data_all <- readRDS(paste0(proj_dir, "/data/snrna/care_mut_umi_data_all_20230729.RDS"))

# Malignant state assignment  
malignant_md <- read.delim(paste0(proj_dir, "/results/scoring/malignant_first_round_cell_state_assignment.txt"), sep="\t", header = TRUE)
malignant_md <- malignant_md %>% 
  mutate(CellType = recode(State, `MP1_OPC_MUT` = "OPC",
                             `MP2_AC_MUT` = "AC",
                             `MP7_AC2_MUT` = "AC",
                             `MP5_MES_MUT` = "MES",
                             `Undifferentiated` = "Undifferentiated")) 

# Restrict to only Undifferentiated malignant cells.
md <- malignant_md %>% 
  filter(CellType == "Undifferentiated")
####################################################################################################################################
####################################################################################################################################
# Prepare expression matrices for NMF computation
#
# Prepare an expression matrix PER SAMPLE FOR EACH CELLTYPE by centering the expression values and setting negative
# values to zero. We the save the matrices as RDS files to enable running the NMF algorithm using HPC.
#
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# md is the meta-data table that contains for each cell the CellID, Sample and CellType variables
# umi_data_all is a named list containing the UMI counts matrix for each sample (the list is named by sample names)

celltypes <- unique(md$CellType)
length(celltypes)

md$Sample <- md$SampleID
samples <- unique(md$Sample)
length(samples)

# The matrices will be saved to this folder
nmf_matrices_path <- file.path(proj_dir, "processed_data/nmf_input/")

####################################################################################################################################
# Helper function that normalizes UMI counts to per-million
####################################################################################################################################
umi2upm <- function(m) {
  count_sum <- colSums(m)
  upm_data <- (t(t(m)/count_sum)) * 1e+06
  upm_data
}


####################################################################################################################################
# Compute the expression matrices for NMF per sample/cell type
####################################################################################################################################
start_time <- Sys.time()
for(ct in celltypes) {
  
  save_dir <- paste0(nmf_matrices_path, ct, "/")
  
  if(!dir.exists(paths = save_dir)) {
    print(paste0("Creating directory for ", ct))
    dir.create(path = save_dir)
  }
  
  print("Preparing expression matrices for NMF")
  
  for(i in 1:length(samples)) {
    
    gc()
    
    sname <- samples[i]
    
    print(paste0("******************************* ", sname, " - start, i = ", i, " *******************************"))
    
    d <- md %>%
      filter(Sample == sname, CellType == ct)
    
    # If there are less than 10 cells we do not generate the matrix  
    if(nrow(d) < 10) {
      print(paste0("Less than 10 ", ct, " cells found for sample ", sname))
      next
    }
    
    # select the UMI counts matrix for the current sample, subset for the specific cell type and normalize counts
    m <- umi_data_all[[sname]]
    m <- m[, d$CellID]
    m <- umi2upm(m)
    dim(m)
    
    print(paste0("Found ", ncol(m), " cells"))
    
    # Select the highly-expressed genes for NMF (should be ~8000 genes)
    rm <- log2(rowMeans(m) + 1)
    
    genes <- names(rm[rm > 4])
    length(genes)
    
    print(paste0("Found ", length(genes), " highly expressed genes"))
    
    # log-transform and center    
    genes <- genes[genes %in% rownames(m)]
    m <- m[genes, ]
    m <- log2(m / 10 + 1)
    m <- apply(m, 1, function(x) x - mean(x))
    m <- t(m)
    dim(m)
    
    # Zero negative values    
    m[m < 0] <- 0
    dim(m)
    
    # Remove all-zero rows/columns (if any) to avoid NMF errors    
    m <- m[rowSums(m) > 0, ]
    m <- m[, colSums(m) > 0]
    dim(m)
    
    print("Saving")
    st <- Sys.time()
    saveRDS(object = m, file = paste0(save_dir, sname, ".RDS"))
    et <- Sys.time()
    print(et - st)
    
    print(paste0("******************************* ", sname, " - end, i = ", i, " *******************************"))
  }
}
end_time <- Sys.time()
end_time - start_time

### END ###