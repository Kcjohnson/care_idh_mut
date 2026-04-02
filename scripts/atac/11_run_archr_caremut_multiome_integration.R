##############################
### Run ArchR RNA and ATAC Multiome integration
### Author: Kevin Johnson
### Updated: 2025.10.20
##############################

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

CARE_filt_rna_malignant <- loadArchRProject("Save-CAREmut-Malignant-RNA")

## Specify the RNA files to be loaded in - these will be the same as the ATAC libraries but selected for RNA
rna_files <- c("/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL01-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL01-2/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL02-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL02-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL02-2/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL03-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL03-2/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL03-3/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL04-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL04-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL05-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL05-2/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL06-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL06-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL07-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL07-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL08-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL08-2/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL09-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL09-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL10-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/NL10-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN01-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN01-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN02-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN02-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN03-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN03-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN04-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN04-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN05-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN05-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN06-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN06-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN07-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN07-2/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN08-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN08-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN10-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN10-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN13-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN13-2/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN15-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN15-2/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN16-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN16-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN17-0/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch1_arc_2/SN17-1/outs/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200891-SC2200892_T0443/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200893-SC2200894_T0902/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200895-SC2200896_T0500/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200897-SC2200898_T0897/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200915-SC2200916_T0816/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200917-SC2200918_T0907/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201049-SC2201050_T0615/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201051-SC2201052_T0772/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201053-SC2201054_T0784/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201055-SC2201056_T1070/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201078-SC2201079_T0509/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201080-SC2201081_T1033/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201082-SC2201083_T0824/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201084-SC2201085_T1069/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC22001086R-SC2201087_T0871/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201088-SC2201089_T0887/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201090-SC2201091_T0973/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201092-SC2201093_T1020/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201094-SC2201095_T0870/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201096-SC2201097_T0925/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201098-SC2201099_T0797/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201100-SC2201101_T0893/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201102-SC2201103_T1059/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201104-SC2201105_T1087/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201106-SC2201107_T0645/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201108-SC2201109_T0848/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201118-SC2201119_T0841/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201120-SC2201121_T0861/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201110-SC2201111_T0609/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201112-SC2201113_T0899/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201122-SC2201123_T0856/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC22001124-SC2201125R_T0896/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201114-SC2201115_T0978/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201116-SC2201117_T0989/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200919-SC2200920_T1012/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200921-SC2200922_T1041/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200833-SC2200834_GT21-03642/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200791-SC2200792_GT21-03649/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200795-SC2200796_GT21-03651/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200793-SC2200794_GT21-03650/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200837-SC2200838_GT21-03653/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200839-SC2200840_GT21-03654/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200867-SC2200868_GT21-03657/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200869-SC2200870_GT21-03659/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200871-SC2200872_GT21-03660/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201015-SC2201016_GT21-03662/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201017-SC2201018_GT21-03663/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200924-SC2200925_GT21-03664/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200926-SC2200927_GT21-03667/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200928-SC2200929_GT21-03668/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201019-SC2201020_GT21-03673/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201021-SC2201022_GT21-03674/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201023-SC2201024_GT21-03679/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201025-SC2201026_GT21-03680/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201027-SC2201028_GT21-03682/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201029-SC2201030_GT21-03683/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201031-SC2201032_GT21-03688/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201033-SC2201034_GT21-03689/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201035-SC2201036_GT21-03694/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2201037-SC2201038_GT21-03695/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200930-SC2200931_GT21-03703/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200932-SC2200933_GT21-03704/cellranger-arc/filtered_feature_bc_matrix.h5",
                 "/gpfs/gibbs/pi/verhaak/shared/monitor/other_proc_data/batch2/SC2200934-SC2200935_GT21-03705/cellranger-arc/filtered_feature_bc_matrix.h5")


sample_names <- c("NL01-1",
                  "NL01-2",
                  "NL02-0",
                  "NL02-1",
                  "NL02-2",
                  "NL03-1",
                  "NL03-2",
                  "NL03-3",
                  "NL04-0",
                  "NL04-1",
                  "NL05-1",
                  "NL05-2",
                  "NL06-0",
                  "NL06-1",
                  "NL07-0",
                  "NL07-1",
                  "NL08-1",
                  "NL08-2",
                  "NL09-0",
                  "NL09-1",
                  "NL10-0",
                  "NL10-1",
                  "SN01-0",
                  "SN01-1",
                  "SN02-0",
                  "SN02-1",
                  "SN03-0",
                  "SN03-1",
                  "SN04-0",
                  "SN04-1",
                  "SN05-0",
                  "SN05-1",
                  "SN06-0",
                  "SN06-1",
                  "SN07-1",
                  "SN07-2",
                  "SN08-0",
                  "SN08-1",
                  "SN10-0",
                  "SN10-1",
                  "SN13-0",
                  "SN13-2",
                  "SN15-0",
                  "SN15-2",
                  "SN16-0",
                  "SN16-1",
                  "SN17-0",
                  "SN17-1",
                  "NL11-0",
                  "NL11-1",
                  "NL12-0",
                  "NL12-1",
                  "NL13-0",
                  "NL13-1",
                  "NL14-0",
                  "NL14-1",
                  "NL15-0",
                  "NL15-1",
                  "NL16-0",
                  "NL16-1",
                  "NL17-0",
                  "NL17-1",
                  "NL18-0",
                  "NL18-1",
                  "NL19-0",
                  "NL19-1",
                  "NL20-0",
                  "NL20-1",
                  "NL21-0",
                  "NL21-1",
                  "NL22-0",
                  "NL22-1",
                  "NL23-0",
                  "NL23-1",
                  "NL24-0",
                  "NL24-1",
                  "NL25-0",
                  "NL25-1",
                  "NL26-0",
                  "NL26-1",
                  "NL27-0",
                  "NL27-1",
                  "NL28-0",
                  "NL28-1",
                  "SJ01-1",
                  "SJ03-0",
                  "SJ03-1",
                  "SJ03-2",
                  "SJ04-1",
                  "SJ04-2",
                  "SJ05-1",
                  "SJ06-1",
                  "SJ06-2",
                  "SJ07-1",
                  "SJ07-2",
                  "SJ08-0",
                  "SJ08-2",
                  "SJ08-3",
                  "SJ10-1",
                  "SJ10-2",
                  "SJ12-1",
                  "SJ12-2",
                  "SJ13-1",
                  "SJ13-2",
                  "SJ15-1",
                  "SJ15-2",
                  "SJ17-1",
                  "SJ17-2",
                  "SJ20-1",
                  "SJ20-2",
                  "SJ20-3")

names(rna_files) <- sample_names


# 1. Import each RNA feature matrix
rna_list <- lapply(names(rna_files), function(smp) {
  import10xFeatureMatrix(
    input = rna_files[[smp]],
    names = smp,
    strictMatch = FALSE
  )
})
names(rna_list) <- names(rna_files)

# 2. Combine them (optional — see below)
# You can either add one at a time, or merge into one combined object
# If you want a single seRNA object:
common_genes <- Reduce(intersect, lapply(rna_list, rownames))
rna_counts <- lapply(rna_list, function(se) assay(se, "counts")[common_genes, , drop = FALSE])
rna_combined <- do.call(cbind, rna_counts)

# Wrap back into a SummarizedExperiment - this ran into issues
#rna_combined_se <- rna_list[[1]][common_genes, ]
#assay(rna_combined_se, "counts") <- rna_combined

# Use shared gene order and metadata
gene_metadata <- rowData(rna_list[[1]])[common_genes, , drop = FALSE]

# Construct a clean SummarizedExperiment
rna_combined_se <- SummarizedExperiment(
  assays = list(counts = rna_combined),
  rowData = gene_metadata
)

# Quick sanity check
dim(rna_combined_se)
length(rownames(rna_combined_se))
length(colnames(rna_combined_se))

# Use gene ranges from one of the imported 10x objects
rowRanges(rna_combined_se) <- rowRanges(rna_list[[1]])[rownames(rna_combined_se)]

# Most of the cells are there from RNA (70,556 out of 71,088).
cellsToKeep <- which(getCellNames(CARE_filt_rna_malignant) %in% colnames(rna_combined_se))
cellsSample <- getCellNames(CARE_filt_rna_malignant)[cellsToKeep]

# The cells seem to be distributed acros many samples and the missing cells are therefore not a mapping error.
cellsMissing <- which(!getCellNames(CARE_filt_rna_malignant) %in% colnames(rna_combined_se))
getCellNames(CARE_filt_rna_malignant)[cellsMissing]

# Create a subset of the ArchR object to save malignant-only analyses.
projMulti <- subsetArchRProject(ArchRProj = CARE_filt_rna_malignant, cells = cellsSample, outputDirectory = "Save-ArchR-Multiome-Analysis", force = TRUE)

# Load the ArchR project restricted to malignant cells that passed the ARC pipeline
projMulti <- loadArchRProject("Save-ArchR-Multiome-Analysis")


# Add the gene expression matrix to this ArchR project
projMulti2 <- addGeneExpressionMatrix(
  input = projMulti,               
  seRNA = rna_combined_se,
  force = TRUE
)

# Confirming that the cell state grouping is properly recorded.
table(projMulti2$CellStateGroup)

# Check to make sure that the expected matrices have already been determined
getAvailableMatrices(projMulti2) 

# Not sure what happened with copying the arrow files to this new project but since this returns "Not All Seqnames Identical", I am moving to "GeneScoreMatrix" approach.
projMulti2 <- addIterativeLSI(
  ArchRProj = projMulti2, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC"
)

# Repeat for RNA - may go to a smaller number of features since it likely picks up sample specific chromatin/expression features at a point

projMulti2 <- addIterativeLSI(
  ArchRProj = projMulti2,
  useMatrix = "GeneScoreMatrix", 
  name = "LSI_ATAC", 
  iterations = 2, 
  clusterParams = list( 
    resolution = c(0.2), 
    sampleCells = 20000, 
    n.start = 10
  ), 
  varFeatures = 2500, 
  dimsToUse = 1:20, 
  force = TRUE
)


projMulti2 <- addIterativeLSI(
  ArchRProj = projMulti2, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)


projMulti2 <- addCombinedDims(projMulti2, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.4, force = TRUE)
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 0.4, force = TRUE)
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_Combined", name = "Clusters_Combined", resolution = 0.4, force = TRUE)

cols <- c("#AA2756", "#F77D58", "#7fbf7b", "#E8F5A3", "gray90")
names(cols)  <- names(table(projMulti2$CellStateGroup)) 

p1 <- plotEmbedding(projMulti2, name = "CellStateGroup", embedding = "UMAP_ATAC", size = 1, labelAsFactors=F, labelMeans=F, pal = cols)
p2 <- plotEmbedding(projMulti2, name = "CellStateGroup", embedding = "UMAP_RNA", size = 1, labelAsFactors=F, labelMeans=F, pal = cols)
p3 <- plotEmbedding(projMulti2, name = "CellStateGroup", embedding = "UMAP_Combined", size = 1, labelAsFactors=F, labelMeans=F, pal = cols)

p <- lapply(list(p1,p2,p3), function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

do.call(cowplot::plot_grid, c(list(ncol = 3),p))



p_state_rna <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "CellStateGroup", embedding = "UMAP_RNA", size = 1, pal = cols, labelMeans=F) 
p_state_rna

# Highly dependent on sample
p_sample_rna <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "Sample", embedding = "UMAP_RNA", size = 1, labelMeans=F) 
p_sample_rna

p_state <- plotEmbedding(ArchRProj = projMulti2, colorBy = "cellColData", name = "CellStateGroup", embedding = "UMAP_ATAC", size = 1, pal = cols, labelMeans=F) 
p_state


# We don't need to re-run these steps because we've already run these steps in a prior script.
# pathToMacs2 <- findMacs2()
# projMulti2 <- addGroupCoverages(ArchRProj = projMulti2, groupBy = "CellStateGroup", verbose = FALSE, force=TRUE)
# projMulti2 <- addReproduciblePeakSet(ArchRProj = projMulti2, groupBy = "CellStateGroup", pathToMacs2 = pathToMacs2, force=TRUE)
# projMulti2 <- addPeakMatrix(ArchRProj = projMulti2, force=TRUE)

# It says "unused force=TRUE", which is odd because I thought this was used in that last one.
# projMulti2 <- addPeak2GeneLinks(ArchRProj = projMulti2, reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix", force=TRUE)
projMulti2 <- addPeak2GeneLinks(ArchRProj = projMulti2, reducedDims = "LSI_Combined", useMatrix = "GeneExpressionMatrix")

# Inspect some of the markers:
se <- getMarkerFeatures(ArchRProj = projMulti2,
                        groupBy = "CellStateGroup",
                        bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"))

heatmap_gex <- plotMarkerHeatmap(
  seMarker = se, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  nLabel = 4,
  transpose = TRUE
)

# Extract some peak2gene links
p2g <- getPeak2GeneLinks(
  ArchRProj = projMulti2,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = TRUE
)

p2g[[1]]

markerGenes  <- c(
  "AQP4", # AC-like
  "PDGFRA", # OPC-like
  "OLIG1",
  "HOXD11",
  "CD44",
  "VIM",
  "ANXA2",
  "DLL3")

p <- plotBrowserTrack(
  ArchRProj = projMulti2, 
  groupBy = "CellStateGroup", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(projMulti2)
)
grid::grid.newpage()
grid::grid.draw(p$AQP4)

grid::grid.newpage()
grid::grid.draw(p$OLIG1)

grid::grid.newpage()
grid::grid.draw(p$PDGFRA)

grid::grid.newpage()
grid::grid.draw(p$CD44)

grid::grid.newpage()
grid::grid.draw(p$VIM)

grid::grid.newpage()
grid::grid.draw(p$ANXA2)

grid::grid.newpage()
grid::grid.draw(p$DLL3)

# Print these to the ArchR project
plotPDF(plotList = p, 
        name = "Plot-Tracks-Key-Genes-with-Peak2GeneLinks.pdf", 
        ArchRProj = projMulti2, 
        addDOC = FALSE, width = 5, height = 5)


# Re-run with returnLoops since this causes errors for some reason downstream
p2g <- getPeak2GeneLinks(
  ArchRProj = projMulti2,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

# Getting both the gene and peak name
p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
p2g$peakName <- (metadata(p2g)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g$idxATAC]
# Confirm that it's working as expected
p2g

metadata(p2g)$seRNA

# This seems to be about the best one could hope for in terms of separation
p2g_heat <- plotPeak2GeneHeatmap(ArchRProj = projMulti2, 
                               groupBy = "CellStateGroup",
                               palGroup=cols,
                               k = 10,
                               corCutOff = 0.45,             
                               varCutOffATAC = 0.25,
                               varCutOffRNA = 0.25,
                               nPlot = 25000)

p2g_heat

p2g_heat_patient <- plotPeak2GeneHeatmap(ArchRProj = projMulti2, 
                               groupBy = "patient_id",
                               k = 10,
                               corCutOff = 0.45,             
                               varCutOffATAC = 0.25,
                               varCutOffRNA = 0.25,
                               nPlot = 20000)

p2g_heat_patient

plotPDF(p2g_heat, name = "Plot-Peak2GeneLinks-Heatmap.pdf", ArchRProj = projMulti2, addDOC = FALSE, width = 5, height = 5)
plotPDF(p2g_heat_patient, name = "Plot-Peak2GeneLinks-Heatmap-Patient.pdf", ArchRProj = projMulti2, addDOC = FALSE, width = 5, height = 5)


# Check point to save
projMulti2 <- saveArchRProject(ArchRProj = projMulti2, outputDirectory = "Save-ArchR-Multiome-Analysis", overwrite = TRUE, load = TRUE)


# The analysis below aims to determine peak2gene relationships for each cell state.
# Need to store the results following each cell state and timepoint.
cell_states <- unique(projMulti2$CellStateGroup)
timepoints  <- c("T1", "T2")

# Keep per-subset results
p2g_links_by_subset <- list()

for (state in cell_states) {
  for (tp in timepoints) {
    
    message("P2G for ", state, " @ ", tp)
    
    cells_subset <- getCellNames(projMulti2)[
      projMulti2$CellStateGroup == state & projMulti2$timepoint == tp
    ]
    
    if (length(cells_subset) < 100) {
      message("Skipping ", state, " ", tp, " (", length(cells_subset), " cells)")
      next
    }
    
    # compute links within this subset
    projMulti2 <- addPeak2GeneLinks(
      ArchRProj  = projMulti2,
      reducedDims = "LSI_Combined", 
      useMatrix = "GeneExpressionMatrix",
      cellsToUse  = cells_subset,
      verbose     = TRUE,
      logFile     = paste0("ArchRLogs/P2G_", state, "_", tp, ".log")
    )
    
    
    # immediately retrieve and stash the links produced by the call above
    p2g_now <- getPeak2GeneLinks(
      ArchRProj     = projMulti2,
      corCutOff     = 0.0,        # pull all; filter later
      varCutOffATAC = 0.0,
      varCutOffRNA  = 0.0,
      returnLoops   = FALSE
    )
    
    p2g_now$GeneSymbol <- mcols(metadata(p2g_now)$geneSet)$name[p2g_now$idxRNA]
    p2g_now$Peak <- (metadata(p2g_now)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2g_now$idxATAC]
    
    if (!is.null(p2g_now) && nrow(as.data.frame(p2g_now)) > 0) {
      df <- as.data.frame(p2g_now) %>%
        transmute(
          peakName = Peak,
          geneName = GeneSymbol,
          Correlation = Correlation,
          FDR = FDR
        ) %>%
        mutate(cell_state = state, timepoint = tp)
      key <- paste(state, tp, sep = "_")
      p2g_links_by_subset[[key]] <- df
    }
  }
}

### Compare T2 vs T1 within each cell state ### 
p2g_diff_all <- bind_rows(p2g_links_by_subset)

p2g_diff_list <- lapply(cell_states, function(state) {
  t1 <- p2g_links_by_subset[[paste0(state, "_T1")]]
  t2 <- p2g_links_by_subset[[paste0(state, "_T2")]]
  if (is.null(t1) || is.null(t2)) return(NULL)
  
  full_join(t1, t2, by = c("peakName", "geneName"),
            suffix = c("_T1", "_T2")) %>%
    mutate(
      corChange = Correlation_T2 - Correlation_T1,
      direction = case_when(
        # Peak present only at T2 → new / gained accessibility
        is.na(Correlation_T1) & !is.na(Correlation_T2) ~ "New Peak / Gain",
        
        # Peak present only at T1 → lost accessibility
        !is.na(Correlation_T1) & is.na(Correlation_T2) ~ "Lost Peak / Loss",
        
        # Peak present in both, correlation increased
        !is.na(corChange) & corChange >  0.25 ~ "Increased",
        
        # Peak present in both, correlation decreased
        !is.na(corChange) & corChange < -0.25 ~ "Decreased",
        
        # Otherwise no meaningful change
        TRUE ~ "Stable"
      ),
      cell_state = state
    )
})

p2g_diff_all <- bind_rows(p2g_diff_list)

# Store output
saveRDS(p2g_links_by_subset, "/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/longitudinal_malignant_p2g_links_by_subset.RDS")
saveRDS(p2g_diff_all, "/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/longitudinal_malignant_p2g_differenences.RDS")

# Load back in since this took quite some time.
p2g_links_by_subset <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/longitudinal_malignant_p2g_links_by_subset.RDS")
p2g_diff_all <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/results/archr/longitudinal_malignant_p2g_differenences.RDS")

# Setting a higher bar for correlation strength since we want to garner information about only the most relevant chromatin-gene relationships.
p2g_diff_all_filtered <- p2g_diff_all %>% 
  filter(Correlation_T1>0.5 | Correlation_T2 > 0.5) 

# Produce a summary for plotting
p2g_summary <- p2g_diff_all_filtered %>%
  dplyr::group_by(cell_state, direction) %>%
  dplyr::summarise(n_links = n(), .groups = "drop") %>%
  dplyr::group_by(cell_state) %>%
  dplyr::mutate(pct_links = 100 * n_links / sum(n_links))

p2g_summary

ggplot(p2g_summary, aes(x=cell_state, y=pct_links, fill=direction)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=c(
    "Increased"="steelblue3",
    "Decreased"="firebrick3",
    "New Peak / Gain"="#2171b5",
    "Lost Peak / Loss"="#a50f15",
    "Stable"="grey80"
  )) +
  labs(x="Cell State", y="Fraction of P2G Links", fill="Direction")

### ### ### ### ### ### ### ### ### ###
# Read in the differentially expressed genes from longitudinal pseudobulk analyses
### ### ### ### ### ### ### ### ### ###

library(dplyr)
library(readr)
library(purrr)
library(stringr)

# Define the files where pseudobulk results were stored
files <- c(
  "/vast/palmer/pi/verhaak/kcj28/care_mut/results/pseudobulk/pseudobulk_ac_t1t2_deg.txt",
  "/vast/palmer/pi/verhaak/kcj28/care_mut/results/pseudobulk/pseudobulk_opc_t1t2_deg.txt",
  "/vast/palmer/pi/verhaak/kcj28/care_mut/results/pseudobulk/pseudobulk_npc_t1t2_deg.txt",
  "/vast/palmer/pi/verhaak/kcj28/care_mut/results/pseudobulk/pseudobulk_mes_t1t2_deg.txt",
  "/vast/palmer/pi/verhaak/kcj28/care_mut/results/pseudobulk/pseudobulk_undiff_t1t2_deg.txt"
)

# Helper to extract cell type from filename
get_celltype <- function(path) {
  str_match(basename(path), "pseudobulk_(.*?)_t1t2_deg")[,2]
}

# Define the mapping  between the different naming schema
celltype_map <- c(
  "ac" = "AC",
  "opc" = "OPC",
  "npc" = "NPC",
  "mes" = "MES",
  "undiff" = "Undifferentiated"
)

# Read and combine all pseudobulk DEG tables
pseudobulk_all <- map_dfr(files, function(f) {
  raw_type <- get_celltype(f)
  read.delim(f, check.names = FALSE) %>%
    mutate(
      cell_type_raw = raw_type,
      cell_type = celltype_map[raw_type],
      source_file = basename(f)
    )
})

# Quick sanity check
table(pseudobulk_all$cell_type)

### Examine overlap and test for enrichment ###

up_genes_by_state <- pseudobulk_all %>%
  filter(log2FoldChange > 0, padj < 0.1) %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(
    up_genes = list(unique(feature)),
    .groups = "drop"
  )

up_genes_by_state

p2g_collapsed <- p2g_diff_all_filtered %>%
  dplyr::group_by(cell_state, geneName) %>%
  dplyr::summarise(
    direction = names(sort(table(direction), decreasing = TRUE))[1],
    .groups = "drop"
  )


overlap_summary <- p2g_collapsed %>%
  dplyr::group_by(cell_state, direction) %>%
  dplyr::summarise(
    n_total = n(),
    n_overlap = sum(geneName %in% unlist(
      up_genes_by_state$up_genes[up_genes_by_state$cell_type == unique(cell_state)]
    )),
    prop_overlap = n_overlap / n_total,
    .groups = "drop"
  ) %>%
  arrange(cell_state, prop_overlap)

ggplot(overlap_summary, aes(x = reorder(direction, prop_overlap), y = prop_overlap, fill = direction)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~cell_state, scales = "free_y") +
  coord_flip() +
  labs(
    x = "P2G grouping",
    y = "Proportion overlapping with T2 upregulated genes",
    title = "Overlap between chromatin change groupings and T2 upregulated DEGs"
  ) +
  theme_minimal()


deg_long <- pseudobulk_all %>%
  filter(!is.na(log2FoldChange), !is.na(padj), group=="T2") %>%
  mutate(deg_direction = ifelse(log2FoldChange > 0, "Up_T2", "Down_T2"),
         state_gene = paste0(cell_type, "_", feature),
         deg_status = dplyr::case_when(
           deg_direction == "Up_T2" & padj<0.1 ~ "Upregulated",
           deg_direction == "Down_T2" & padj<0.1 ~ "Downregulated",
             padj > 0.1 ~ "Stable"))

p2g_collapsed_merge <- p2g_collapsed %>% 
  mutate(state_gene = paste0(cell_state, "_", geneName)) 

deg_long_p2g <- deg_long %>% 
  left_join(p2g_collapsed_merge, by="state_gene") %>% 
  mutate(p2g_direction = ifelse(is.na(direction), "no_peak_to_gene", direction))

cor.test(df_tierB$Correlation_T1, df_tierB$log2FoldChange, method = "spearman")


table(deg_long_p2g$deg_direction, deg_long_p2g$p2g_direction)

p2g_summary <- deg_long_p2g %>%
  dplyr::group_by(deg_status, p2g_direction) %>%
  dplyr::summarise(counts = n()) %>% 
  dplyr::mutate(prop = counts / sum(counts)) %>%
  ungroup()


ggplot(p2g_summary, aes(x = deg_status, y = prop, fill = p2g_direction)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(
    "Increased" = "#1b9e77",
    "Decreased" = "#d95f02",
    "Stable" = "#7570b3",
    "Lost Peak / Loss" = "#e7298a",
    "New Peak / Gain" = "#66a61e",
    "no_peak_to_gene" = "grey70"
  )) +
  labs(
    x = "Differential Expression Category",
    y = "Proportion of genes",
    fill = "Peak-to-Gene Change",
    title = "Association between differential expression and chromatin linkage categories"
  ) +
  theme_minimal(base_size = 13) 

tab <- p2g_summary %>%
  select(deg_status, p2g_direction, counts) %>%
  tidyr::pivot_wider(
    names_from = p2g_direction,
    values_from = counts,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("deg_status") %>%
  as.matrix()

chisq.test(tab)


fisher_results <- p2g_summary %>%
  dplyr::group_by(p2g_direction) %>%
  dplyr::group_modify(~{
    mat <- matrix(
      c(
        .x$counts[.x$deg_status == "Upregulated"],
        sum(.x$counts[.x$deg_status == "Upregulated"]),
        .x$counts[.x$deg_status == "Stable"],
        sum(.x$counts[.x$deg_status == "Stable"])
      ),
      nrow = 2
    )
    test <- fisher.test(mat)
    tibble(
      direction = unique(.x$p2g_direction),
      pval = test$p.value,
      odds_ratio = test$estimate
    )
  }) %>%
  mutate(padj = p.adjust(pval, method = "BH"))

### END ###