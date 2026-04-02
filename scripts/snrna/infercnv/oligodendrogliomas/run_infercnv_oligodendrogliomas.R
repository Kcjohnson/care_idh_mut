# run_inferCNV.R - Runs inferCNV for RNAseq data pre-processed by Seurat

library(infercnv)
library(optparse)

# Input arguments
option_list <- list(
  make_option(c("-i","--SAMPLE_ID"), action="store", default=NULL, type="character",
              help="SAMPLE_ID (e.g. SJ02-1)"),
  make_option(c("-t","--TEMP_DIR"), action="store", default=NULL, type="character",
              help="Temporary directory."),
  make_option(c("-o","--OUTDIR"), action="store", default=NULL, type="character",
              help="Output directory."),
  make_option(c("-n","--num_cores"), action="store", default=NULL, type="integer",
              help="Number of cores used for multithreading-compatible functions.")
)

args <- parse_args(OptionParser(option_list=option_list))
sample_id <- args$SAMPLE_ID
temp_dir <- args$TEMP_DIR
outdir <- args$OUTDIR
num_cores <- args$num_cores

# Test placeholders
# sample_id <- "SJ02-1"
# outdir <- "/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/infercnv_outputs/"

##### Input arguments #####
# Sample-level counts matrix.
counts <- readRDS(paste0("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/infercnv_inputs/oligodendrogliomas/",sample_id,"_counts.RDS"))

# Sample-level sample annotation file post QC and cell state assignment.
anno_file <- paste0("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/infercnv_inputs/oligodendrogliomas/",sample_id,"_sample_annotation_file.txt")

# Load annotation file.
anno <- read.delim(anno_file, header = FALSE, stringsAsFactors = FALSE)

## Filter the counts object so that only one sample per patient is analyzed.
counts_filt <- counts[,colnames(counts)%in%anno$V1]

print("Do all cell IDs match between annotation and filtered counts?")
all(colnames(counts_filt)==anno$V1)

# Extract the cell types present.
cell_types <- unique(anno$V2)

## Tabulate the number of cells per state
print("Breakdown of cell types:")
table(anno$V2)

# Extract the reference (nontumor) cell types present
cell_types.reference <- names(table(anno$V2))[names(table(anno$V2))%in%c("Oligodendrocyte", "Myeloid")]
sprintf("Reference population(s) will be: %s", cell_types.reference)

###########################
# Create the infercnv object.
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_filt),
                                    annotations_file=anno_file,
                                    delim="\t",
                                    gene_order_file="/vast/palmer/pi/verhaak/kcj28/reference/infercnv/hg38_gencode_v27.txt",
                                    ref_group_names=cell_types.reference)


# Perform infercnv operations to reveal cnv signal.
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 0.1 for 10x-genomics
                             out_dir=paste0(outdir,"infercnv_",sample_id),  # dir is auto-created for storing outputs
                             num_ref_groups=length(cell_types.reference),
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
                             num_threads = num_cores
)


# For further denoising, you can apply median filtering, but this might remove true signal from the plot.
outdir_median <- paste0(outdir,"infercnv_",sample_id, "/median_filtered")
infercnv_obj_median_filtered = infercnv::apply_median_filtering(infercnv_obj)

print("Plotting median filtered heatmap")
infercnv::plot_cnv(infercnv_obj_median_filtered,
                   out_dir = outdir_median,
                   output_filename = 'infercnv.median_filtered',
                   x.range = "auto",
                   x.center = 1,
                   title = "infercnv_median_filtered",
                   color_safe_pal = FALSE)

### END ####