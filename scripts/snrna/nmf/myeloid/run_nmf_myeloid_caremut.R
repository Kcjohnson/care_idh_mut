#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

print("************************************ Starting ************************************")

library(Matrix)

## Test parameters:
# in_path <- "/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/nmf_2026/Myeloid/"
# mat_index <- 1
# out_path <- "/vast/palmer/pi/verhaak/kcj28/care_mut/results/nmf_res_2026/myeloid_n73_no_parallel/"
# rank_lb <- 3
# rank_ub <- 10
# nrun <- 10

in_path <- as.character(args[[1]])   ### The name of a list, where each element is a matrix on which we want to run NMF
mat_index <- as.numeric(args[[2]])    ### The index of the current matrix we want to run here  
out_path <- as.character(args[[3]])   ### The name of a list, where each element is a matrix on which we want to run NMF
rank_lb <- as.numeric(args[[4]])      ### lower rank of NMF
rank_ub <- as.numeric(args[[5]])      ### upper rank of NMF
nrun <- as.numeric(args[[6]])         ### number of runs 

print(paste0("Loading packages"))

library(NMF)
library(BiocParallel)
set.seed(1)

files_list <- list.files(in_path)

filename <- files_list[mat_index]

print(paste0("Openning index ", mat_index, " in directory ", in_path, ". File name is ", filename))

filepath <- paste0(in_path, filename)

nmf_mat  <- readRDS(file = filepath)  ### List with all matrices (samples) 

mat_name <- gsub("\\.RDS", "", filename) ### name of current sample 

nmf_mat <- as.matrix(nmf_mat)     ### current sample we want to run 

print(paste0("Loaded matrix with ", nrow(nmf_mat), " rows and ", ncol(nmf_mat), " columns"))

print(paste0("Proceeding with a matrix with ", nrow(nmf_mat), " rows and ", ncol(nmf_mat), " columns"))

print(paste0("Running NMF for sample ", mat_name, ", ", "LB=", rank_lb, ", UB=", rank_ub, ", nrun=", nrun, ", method=snmf/r"))

# This runs the NMF algorithm. Due to an issue that we encountered on our HPC infrastructure we cannot use parallelization but
# this might not be the case on your system so check it (you'll need to set the .pbackend parameter)
nmf_res <- nmf(x = nmf_mat, rank = rank_lb:rank_ub, method = "snmf/r", nrun = nrun , .opt = "v", .pbackend = NA)

print(paste0("NMF completed successfully"))

out_filename <- paste0(out_path, mat_name, "_rank", rank_lb, "_", rank_ub, "_nrun", nrun, ".RDS")

print(paste0("Saving result to ", out_filename))

### Would be good to save using this nomenclature so that the code for finding robust NMFs and defining MPs runs smoothly: 
saveRDS(nmf_res, file = out_filename)  

print("************************************ Ended ************************************")


