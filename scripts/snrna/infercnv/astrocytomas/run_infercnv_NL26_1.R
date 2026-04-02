# NL26-1 did not have any non-malignant cells identified so I will be using the reference cells from its initial tumor

library(infercnv)


# Sample-level counts matrix for both initial and recurrent tumors
counts_initial <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/infercnv_inputs/astrocytomas/NL26-0_counts.RDS")
counts_recurrence <- readRDS("/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/infercnv_inputs/astrocytomas/NL26-1_counts.RDS")

# Sample-level sample annotation file post QC and cell state assignment.
anno_file_initial <- "/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/infercnv_inputs/astrocytomas/NL26-0_sample_annotation_file.txt"
anno_file_recurrence <- "/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/infercnv_inputs/astrocytomas/NL26-1_sample_annotation_file.txt"

# Load in annotation, restrict to reference cells in the initial and combine with recurrence.
anno_initial <- read.delim(anno_file_initial, header = FALSE, stringsAsFactors = FALSE)
anno_recurrence <- read.delim(anno_file_recurrence, header = FALSE, stringsAsFactors = FALSE)

anno_initial_reference <- anno_initial[anno_initial$V2%in%c("Myeloid", "Oligodendrocyte"), ]
counts_initial_filt <- counts_initial[,colnames(counts_initial)%in%anno_initial_reference$V1]

# These will be the new inputs.
all_annotation <- rbind(anno_initial_reference, anno_recurrence)
counts <- cbind(counts_initial_filt, counts_recurrence)

print("Do all cell IDs match between annotation and input counts?")
all(colnames(counts)==all_annotation$V1)

# Write out the annotation to be loaded back in.
write.table(all_annotation, 
            file = "/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/infercnv_inputs/astrocytomas/NL26_1_plus_ref_sample_annotation_file.txt",
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Combined annotation file
anno_file <- "/vast/palmer/pi/verhaak/kcj28/care_mut/processed_data/rna/infercnv_inputs/astrocytomas/NL26_1_plus_ref_sample_annotation_file.txt"

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
                             out_dir="/vast/palmer/pi/verhaak/kcj28/care_mut/results/infercnv/astrocytoma_samples/infercnv_NL26-1",  # dir is auto-created for storing outputs
                             num_ref_groups=length(cell_types.reference),
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
                             num_threads = 4
)


# For further denoising, you can apply median filtering, but this might remove true signal from the plot.
outdir_median <- paste0("/vast/palmer/pi/verhaak/kcj28/care_mut/results/infercnv/astrocytoma_samples/infercnv_NL26-1", "/median_filtered")
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