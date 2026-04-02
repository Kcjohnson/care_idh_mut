##############################
### Run ArchR analyses on CARE IDH-mutant multiome data (48 samples)
### Author: Kevin Johnson
### Updated: 2026.03.30
##############################

# ArchR creates several directories automatically when creating arrow files and ArchR projects.
workdir <- "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/results/atac"
setwd(workdir)

# Load necessary packages
library(ArchR)
library(parallel)
library(pheatmap)
library(chromVARmotifs)


#### Set-up #####

# Check available cores. ArchR recommends setting total number of cores 1/2 to 3/4 of available cores,
ncores <- detectCores() 
nthreads <- ncores/2
addArchRThreads(threads = nthreads) 
# Each R session requires that the genome is also specified and must match alignment.
addArchRGenome("hg38")

# These samples were processed across two batches
batch1_dirs <- list.files("/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1")
batch2_dirs <- list.files("/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2")
all_dirs <- c(batch1_dirs, batch2_dirs)

message("Number of unique samples found: ", length(unique(all_dirs)))

## Create Arrow files using the desired input. It can be fragments or bam files. Choosing fragments here.
  # NL01
input_files <- c("/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL01-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL01-2/outs/atac_fragments.tsv.gz",
	# NL03
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL03-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL03-2/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL03-3/outs/atac_fragments.tsv.gz",
	# NL04
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL04-0/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL04-1/outs/atac_fragments.tsv.gz",
	# NL05
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL05-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/NL05-2/outs/atac_fragments.tsv.gz",
	# SN05
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/SN05-0/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/SN05-1/outs/atac_fragments.tsv.gz",
  # SN07
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/SN07-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/SN07-2/outs/atac_fragments.tsv.gz",
	# SN17
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/SN17-0/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch1/SN17-1/outs/atac_fragments.tsv.gz",
  # NL11
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL11-0/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL11-1/outs/atac_fragments.tsv.gz",
	# NL12
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL12-0/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL12-1/outs/atac_fragments.tsv.gz",
	# NL23
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL23-0/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL23-1/outs/atac_fragments.tsv.gz",
	# NL26
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL26-0/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/NL26-1/outs/atac_fragments.tsv.gz",
	# SJ03
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ03-0/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ03-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ03-2/outs/atac_fragments.tsv.gz",
	# SJ04
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ04-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ04-2/outs/atac_fragments.tsv.gz",
	# SJ06
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ06-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ06-2/outs/atac_fragments.tsv.gz",
	# SJ07
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ07-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ07-2/outs/atac_fragments.tsv.gz",
	# SJ08
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ08-0/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ08-2/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ08-3/outs/atac_fragments.tsv.gz",
	# SJ10
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ10-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ10-2/outs/atac_fragments.tsv.gz",
	# SJ12
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ12-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ12-2/outs/atac_fragments.tsv.gz",
	# SJ13
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ13-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ13-2/outs/atac_fragments.tsv.gz",
	# SJ15
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ15-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ15-2/outs/atac_fragments.tsv.gz",
	# SJ17
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ17-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ17-2/outs/atac_fragments.tsv.gz",
	# SJ20
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ20-1/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ20-2/outs/atac_fragments.tsv.gz",
	"/vast/palmer/pi/verhaak/kcj28/cellranger_arc/batch2/SJ20-3/outs/atac_fragments.tsv.gz")

message("Number of unique files found: ", length(unique(input_files)))

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) == 0) {
  message("PASS: All input files exist")
} else {
  message("FAIL: ", length(missing_files), " missing files:")
  message(paste(missing_files, collapse = "\n"))
}

# Use the shortened multiome IDs for the sample names.
sample_names <- basename(dirname(dirname(input_files)))

if (all(sample_names == all_dirs)) {
  message("PASS: All sample names match input directories")
} else {
  mismatches <- which(sample_names != all_dirs)
  message("FAIL: Mismatches found at indices: ", paste(mismatches, collapse = ", "))
  message("Expected: ", paste(all_dirs[mismatches], collapse = ", "))
  message("Got:      ", paste(sample_names[mismatches], collapse = ", "))
}

# Assign these as the names that ArchR will read in
names(input_files) <- sample_names

# Set.seed for any random processes that might be run by ArchR.
set.seed(1)

message("Creating arrow files")

# May want to play around with the `minTSS` variable if it excludes too many cells. Setting to default = 4.
ArrowFiles <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = names(input_files),
  minTSS = 4, # Don't set this too high because you can always increase later. Will differ between cell lines, human samples, and different datasets.
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

message("Calculating Doublet Scores")

# Doublet inference usually detects ~5% or more doublets.
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)


##### Analysis #########
## Creation of an ArchR project that can combine multiple "arrow" files 
## There are different ways to access this ArchR project, see: https://www.archrproject.com/bookdown/manipulating-an-archrproject.html
projCARE <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Save-AllSamples-2026",
  copyArrows = TRUE # This is recommended so that if you modify the Arrow files you have an original copy for later usage.
)

## Save ArchR project:
saveArchRProject(ArchRProj = projCARE, outputDirectory = "Save-AllSamples-2026", load = FALSE)

### END ####