# ----------------------------------------------------------------------------------------------------
# Load packages  
# ----------------------------------------------------------------------------------------------------

library(NMF)
library(ggplot2)
library(scales)
# upload this function
source("/vast/palmer/pi/verhaak/kcj28/care_idh_mut/scripts/utils/custom_magma.R")

# -------
## Use the same variable names that appear below:
#  Genes_nmf_w_basis  - list you produced with gene matrices (all NMFs stacked). Each sample is an entry in the list
#  There are parameters that can be defined using the 'robust_nmf_programs' function, and also after that when finding the clusters  
#  use the same naming for samples in Genes_nmf_w_basis for this to run smoothly (i.e. end each name with '_rank4_9_nruns10.RDS' and for the colnames in each
# matrix '_rank4_9_nrun10.RDS.4.1' where here 4 is the rank and 1 is the program within the rank) 
# -------

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# NMF module public functions
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

load_nmf_results <- function(path, samples = NULL) {
  
  nmf_files <- list.files(path)
  length(nmf_files)
  
  if(!is.null(samples)) {
    snames <- sapply(strsplit(nmf_files, split = "_"), function(x) x[[1]][1])
    nmf_files <- nmf_files[snames %in% samples]
  }
  
  print("Loading...")
  
  nmf_res_list <- lapply(nmf_files, function(x) readRDS(file = paste0(path, x)))
  names(nmf_res_list) <- nmf_files
  
  snames <- sapply(strsplit(names(nmf_res_list), split = "_"), function(x) x[[1]][1])

  Genes_nmf_w_basis <- lapply(1:length(nmf_res_list), function(i) get_nmf_programs(nmf_res_list[[i]]$fit, names(nmf_res_list)[i]))
  names(Genes_nmf_w_basis) <- nmf_files
  length(Genes_nmf_w_basis)
  
  rm(nmf_res_list)
  gc()
  
  return(Genes_nmf_w_basis)
}

# ----------------------------------------------------------------------------------------------------
# Find robust NMFs
# ----------------------------------------------------------------------------------------------------

derive_NMF_metaprograms <- function(Genes_nmf_w_basis, n_genes = 50, save_path = "/vast/palmer/pi/verhaak/kcj28/care_idh_mut/processed_data/nmf_res/", population ="AllCells", plot = T, save = T, verbose = F) {
 
  if(verbose == T)
    print("Starting robust NMF algorithm")
  
  # get gene programs (top 50 genes by NMF score)
  nmf_programs_sig <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:n_genes]))
  nmf_programs_sig <- lapply(nmf_programs_sig, toupper) ## convert all genes to uppercase 
  
  if(verbose == T)
    print("Calling robust_nmf_programs (may take a few minutes)")
  
  # for each sample, select robust NMF programs (i.e. observed using different ranks in the same sample),
  # remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other samples. 
  nmf_filter_ccle <- robust_nmf_programs(nmf_programs_sig, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10) # this can take a few minutes to run 
  nmf_programs_sig <- lapply(nmf_programs_sig, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
  nmf_programs_sig <- do.call(cbind, nmf_programs_sig)
  
  # calculate similarity between programs
  nmf_intersect <- apply(nmf_programs_sig , 2, function(x) apply(nmf_programs_sig , 2, function(y) length(intersect(x,y)))) 
  
  # hierarchical clustering of the similarity matrix 
  nmf_intersect_hc_ccle <- hclust(as.dist(n_genes - nmf_intersect), method="average") 
  nmf_intersect_hc_ccle <- reorder(as.dendrogram(nmf_intersect_hc_ccle), colMeans(nmf_intersect))
  nmf_intersect         <- nmf_intersect[order.dendrogram(nmf_intersect_hc_ccle), order.dendrogram(nmf_intersect_hc_ccle)]
  
  print("Contents of nmf_intersect")
  print(nmf_intersect)
  ### Should save your output
  saveRDS(nmf_intersect, file = paste0(save_path, "nmf_intersect.RDS"))
  saveRDS(nmf_programs_sig, file = paste0(save_path, "nmf_programs_sig.RDS"))
  
  if(verbose == T)
    print("Clustering robust NMF programs")
  
  ### use a clustering approach that updates MPs in each iteration (see SI figure in Pan cancer MP paper)
  
  nmf_intersect_KEEP    <- nmf_intersect
  nmf_programs_sig_KEEP <- nmf_programs_sig
  
  ### Parameters (later change to function form):
  Min_intersect_initial <- 10    # the minimal intersection cutoff for defining the Founder NMF program of a cluster
  Min_intersect_cluster <- 10    # the minimal intersection cuttof for adding a new NMF to the forming cluster 
  Min_group_size        <- 5     # the minimal group size to consider for defining the Founder_NMF of a MP

  Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)
  
  Cluster_list <- list()   ### Every entry contains the NMFs of a chosec cluster
  k <- 1
  Curr_cluster <- c()
  MP_list      <- list()
  
  while (Sorted_intersection[1]>Min_group_size) {   ### CHECK!
    
    Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
    
    ### intersection between all remaining NMFs and Genes in MP 
    Genes_MP                   <- nmf_programs_sig[,names(Sorted_intersection[1])] # initial genes are those in the first NMF. Genes_MP always has only 50 genes consisting of the current MP
    nmf_programs_sig           <- nmf_programs_sig[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs_sig))]  # remove selected NMF
    Intersection_with_Genes_MP <- sort(apply(nmf_programs_sig, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
    NMF_history                <- Genes_MP  # has all genes in all NMFs in the current cluster, for newly defining Genes_MP after adding a new NMF 
    
    ### Create gene list - composed of intersecting genes in descending order + genes with highest NMF scores to add up to 50 genes. Update Curr_cluster each time
    
    while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {   ### Define current cluster 
      
      Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
      
      Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs_sig[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all NMFs in the current cluster 
      Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[n_genes])]   ### genes with overlap equal to the 50th gene
      
      if (length(Genes_at_border)>1){
        ### Sort last genes in Genes_at_border according to maximal NMF gene scores
        ### Run over all NMF programs in Curr_cluster and extract NMF scores for each gene
        Genes_curr_NMF_score <- c()
        for (i in Curr_cluster) {
          curr_study           <- paste( strsplit(i , "[.]")[[1]][1 : which(strsplit(i , "[.]")[[1]] == "RDS")]   , collapse = "."  )
          print(curr_study)
          Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i] 
          # print(length(Q))
          names(Q)             <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])  ### sometimes when adding genes the names do not appear 
          Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
        }
        Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
        Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   ## take only the maximal score of each gene - which is the first entry after sorting
        
        Genes_MP_temp <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[n_genes])]) , names(Genes_curr_NMF_score_sort))
        
      } else {
        Genes_MP_temp <- names(Genes_MP_temp)[1:n_genes] 
      }
      
      NMF_history   <- c(NMF_history , nmf_programs_sig[,names(Intersection_with_Genes_MP)[1]]) 
      Genes_MP <- Genes_MP_temp[1:n_genes]
      
      nmf_programs_sig      <- nmf_programs_sig[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs_sig))]  # remove selected NMF
      
      Intersection_with_Genes_MP <- sort(apply(nmf_programs_sig, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
      
    }
    
    Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
    MP_list[[paste0("MP_",k)]]           <- Genes_MP
    k <- k+1
    
    nmf_intersect             <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]  # remove current chosen cluster
    
    Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   ### Sort intersection of remaining NMFs not included in any of the previous clusters
    
    Curr_cluster <- c()
    # print(dim(nmf_intersect)[2])
  }
  
  #### *****  Sort Jaccard similarity plot according to new clusters:
  
  inds_sorted <- c()
  
  for (j in 1:length(Cluster_list)){
    
    inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_KEEP)))
    
  }
  inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_KEEP)[2],inds_sorted)))) ### combine inds in clusters with inds without clusters
  
  saveRDS(nmf_intersect_KEEP, file = paste0(save_path, "nmf_intersect_KEEP.RDS"))
  saveRDS(inds_new, file = paste0(save_path, "inds_new.RDS"))
  
  # plot re-ordered similarity matrix heatmap     
  nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_KEEP[inds_new,inds_new]) 
  
  if(plot == T) {
    print("Saving plot input")
    
    saveRDS(nmf_intersect_meltI_NEW, file = paste0(save_path, "nmf_intersect_meltI_NEW.RDS"))
    
    if(verbose == T)
      print("Plotting")
    
    p <-
      ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
      geom_tile() + 
      scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
      scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
      theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
      guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

    #print(p)    
    pdf(paste0(save_path, "similarity_heatmap_mp_care_mut_", population, ".pdf"), width=6, height=4)
    print(p)
    dev.off()
    
    # plot similarity matrix cancer type annotation
    print("Loading metadata file from: /vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/caremut_sample_identifier_linker.txt")
    meta_care <- read.table("/vast/palmer/pi/verhaak/kcj28/care_mut/data/metadata/caremut_sample_identifier_linker.txt", sep = "\t", row.names = 1, header = TRUE)
    rownames(meta_care) <- meta_care$SampleID
    sprintf("Available metadata column: %s", colnames(meta_care))

    nmf_annot_type_tumor  <- data.frame("sample"=sub("\\_.*", "", colnames(nmf_intersect_KEEP[inds_new,inds_new])))
    nmf_annot_type_tumor$idh_codel_subtype <- meta_care$idh_codel_subtype[match(nmf_annot_type_tumor$sample, rownames(meta_care))]
    nmf_annot_type_tumor$idh_codel_subtype <- factor(nmf_annot_type_tumor$idh_codel_subtype, levels=sort(unique(nmf_annot_type_tumor$idh_codel_subtype)))
    nmf_annot_type_tumor$timepoint <- meta_care$timepoint[match(nmf_annot_type_tumor$sample, rownames(meta_care))]
    nmf_annot_type_tumor$timepoint <- factor(nmf_annot_type_tumor$timepoint, levels=sort(unique(nmf_annot_type_tumor$timepoint)))
    nmf_annot_type_tumor$lab <- meta_care$lab[match(nmf_annot_type_tumor$sample, rownames(meta_care))]
    nmf_annot_type_tumor$lab <- factor(nmf_annot_type_tumor$lab, levels=sort(unique(nmf_annot_type_tumor$lab)))
    
    saveRDS(nmf_annot_type_tumor, file = paste0(save_path, "nmf_annot_type_tumor.RDS"))
    
    p_annotation_bottom <- ggplot(nmf_annot_type_tumor, aes(y="", x=1:nrow(nmf_annot_type_tumor), fill=idh_codel_subtype, color=idh_codel_subtype)) +
      geom_tile() +
      scale_color_manual(values =c( "IDHmut-codel"= "#2c7bb6", "IDHmut-noncodel"="#d7191c"), name="") +
      scale_fill_manual(values =c( "IDHmut-codel"= "#2c7bb6", "IDHmut-noncodel"="#d7191c"), name="") +                                                    
      theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10),  axis.title = element_text(size=8),  legend.text.align = 0, legend.key.size = unit(0.4, "cm"), legend.key=element_blank(), legend.position=c(1.3,-4), plot.margin = unit(c(0.5,3,-0.6,0.5), "cm")) +
      scale_x_continuous(expand = c(0,0)) +
      labs(x="", y="") +
      guides(fill = "none", color="none")
    
    p_annotation_middle <- ggplot(nmf_annot_type_tumor, aes(y="", x=1:nrow(nmf_annot_type_tumor), fill=lab, color=lab)) +
      geom_tile() +
      scale_color_manual(values =c( "Iavarone lab"= "#7570b3", "Suva lab"="#d95f02", "Verhaak lab"="#1b9e77"), name="") +
      scale_fill_manual(values =c( "Iavarone lab"= "#7570b3", "Suva lab"="#d95f02", "Verhaak lab"="#1b9e77"), name="") +                                                    
      theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10),  axis.title = element_text(size=8),  legend.text.align = 0, legend.key.size = unit(0.4, "cm"), legend.key=element_blank(), legend.position=c(1.3,-4), plot.margin = unit(c(0.5,3,-0.6,0.5), "cm")) +
      scale_x_continuous(expand = c(0,0)) +
      labs(x="", y="") +
      guides(fill = "none", color="none")
    
    p_annotation_top <- ggplot(nmf_annot_type_tumor, aes(y="", x=1:nrow(nmf_annot_type_tumor), fill=timepoint, color=timepoint)) +
      geom_tile() +
      scale_color_manual(values =c( "T1"= "#f1eef6", "T2"="#2b8cbe", "T3"="#045a8d"), name="") +
      scale_fill_manual(values =c( "T1"= "#f1eef6", "T2"="#2b8cbe", "T3"="#045a8d"), name="") +                                                    
      theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 10),  axis.title = element_text(size=8),  legend.text.align = 0, legend.key.size = unit(0.4, "cm"), legend.key=element_blank(), legend.position=c(1.3,-4), plot.margin = unit(c(0.5,3,-0.6,0.5), "cm")) +
      scale_x_continuous(expand = c(0,0)) +
      labs(x="", y="") +
      guides(fill = "none", color="none")
    
    # Combine plots
    pdf(paste0(save_path, "similarity_heatmap_mp_care_mut_annot_", population, ".pdf"), width=6, height=4, onefile=FALSE)
    egg::ggarrange(p_annotation_top, p_annotation_middle, p_annotation_bottom, p, nrow = 4, heights = c(1, 1, 1, 8))
    dev.off()
    
    top_legend_only_plot <- ggplot(nmf_annot_type_tumor, aes(y="", x=1:nrow(nmf_annot_type_tumor), fill=timepoint, color=timepoint)) +
      geom_tile() +
      scale_color_manual(values =c( "T1"= "#f1eef6", "T2"="#2b8cbe", "T3"="#045a8d"), name="") +
      scale_fill_manual(values =c( "T1"= "#f1eef6", "T2"="#2b8cbe", "T3"="#045a8d"), name="") +                                                    
      theme(panel.background = element_rect(fill = "white"),  axis.line = element_blank()) +
      scale_x_continuous(expand = c(0,0)) +
      labs(x="", y="")
    
    middle_legend_only_plot <- ggplot(nmf_annot_type_tumor, aes(y="", x=1:nrow(nmf_annot_type_tumor), fill=lab, color=lab)) +
      geom_tile() +
      scale_color_manual(values =c( "Iavarone lab"= "#7570b3", "Suva lab"="#d95f02", "Verhaak lab"="#1b9e77"), name="") +
      scale_fill_manual(values =c( "Iavarone lab"= "#7570b3", "Suva lab"="#d95f02", "Verhaak lab"="#1b9e77"), name="") +                                                    
      theme(panel.background = element_rect(fill = "white"),  axis.line = element_blank()) +
      scale_x_continuous(expand = c(0,0)) +
      labs(x="", y="")
    
    bottom_legend_only_plot <- ggplot(nmf_annot_type_tumor, aes(y="", x=1:nrow(nmf_annot_type_tumor), fill=idh_codel_subtype, color=idh_codel_subtype)) +
      geom_tile() +
      scale_color_manual(values =c( "IDHmut-codel"= "#2c7bb6", "IDHmut-noncodel"="#d7191c"), name="") +
      scale_fill_manual(values =c( "IDHmut-codel"= "#2c7bb6", "IDHmut-noncodel"="#d7191c"), name="") +                                                    
      theme(panel.background = element_rect(fill = "white"),  axis.line = element_blank()) +
      scale_x_continuous(expand = c(0,0)) +
      labs(x="", y="")
    
    # Save the plot legends 
    ggsave(paste0(save_path, "top_legend_only.pdf"), plot = top_legend_only_plot)
    ggsave(paste0(save_path, "middle_legend_only.pdf"), plot = middle_legend_only_plot)
    ggsave(paste0(save_path, "bottom_legend_only.pdf"), plot = bottom_legend_only_plot)
    
  }
  
  #### SAVE MPs : 
  
  MP <-  do.call(cbind, MP_list)
  
  # Get the current date in the format YYYY-MM-DD
    current_date <- format(Sys.Date(), "%Y-%m-%d")
    sample_num <- length(Genes_nmf_w_basis)
  
  write.csv(MP, file = paste0(save_path, "meta_programs_generated_for_", population, "_n", sample_num, "_min_groupsize", Min_group_size, "_", current_date, ".csv"))
  
  res <- list(clusters = Cluster_list,
              MP_list = MP_list,
              MP_tbl = MP)

  if(verbose == T)
    print("Done!")
  
  return(res)
}

metaprograms_enrichment <- function(mp_list, pathways = NULL, min_gs_size = 40, max_gs_size = 500) {
  
  if(is.null(pathways)) {

    pathways <- default_pathways()
    length(pathways)
  }
  
  pathways <- pathways[lengths(pathways) >= min_gs_size & lengths(pathways) <= max_gs_size]
  length(pathways)
  
  path_en <- lapply(1:length(mp_list), function(i) {
    mp_name <- names(mp_list)[i]
    mp <- mp_list[[i]]
    
    res <- lapply(1:length(pathways), function(j) {
      p <- pathways[[j]]
      tibble(MP = mp_name,
             Pathway = names(pathways)[j],
             N_mp = length(mp),
             N_p = length(p),
             MP_Int = length(intersect(mp, p)),
             MP_Freq = MP_Int / N_mp,
             Jaccard = MP_Int / length(unique(c(mp, p))))
    })
    res <- do.call(rbind, res) 
    return(res)
  })
  path_en <- do.call(rbind, path_en)
  
  path_en$int_genes <- sapply(1:nrow(path_en), function(i) intersect(pathways[[path_en$Pathway[i]]], mp_list[[path_en$MP[i]]]))

  return(path_en)  
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# NMF module utility functions
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

W <- function(x) {
  return(x@fit@W)
}

get_top_genes <- function(fit, rank, sname, n = 50) {
  w <- W(fit)  
  # res <- lapply(1:ncol(w), function(i) sort(w[, i], decreasing = T) %>% head(n = n) %>% names())
  res <- lapply(1:ncol(w), function(i) sort(w[, i], decreasing = T) %>% names())
  res <- do.call(cbind, res)
  # colnames(res) <- paste0("P_", rank, "_", 1:ncol(w))
  colnames(res) <- paste0(sname, ".", rank, ".", 1:ncol(w))
  return(res)
}

get_nmf_basis <- function(fit, rank, sname) {
  w <- W(fit)
  colnames(w) <- paste0(sname, ".", rank, ".", 1:ncol(w))
  w
}

get_nmf_programs <- function(fit, sname, n = 50) {
  
  res <- lapply(names(fit), function(fname) {
    # print(fname)
    f <- fit[[fname]]
    # get_top_genes(f, fname, sname, n = n)
    get_nmf_basis(f, fname, sname)
  })
  res <- do.call(cbind, res)
  # colnames(res) <- paste0(sname, "_", colnames(res))
  return(res)
}

# -------------------------------------------------------------------------------------------
# Function for selecting robust nonnegative matrix factorization (NMF) programs
# ------------------------------------------------------------------------------------------- 

# - nmf_programs = a list; each element contains a matrix with NMF programs (top 50 genes) generated for a specific cell line using different NMF factorization ranks. 
# - intra_min = minimum overlap with a program from the same cell line (for selecting robust programs)
# - intra_max = maximum overlap with a program from the same cell line (for removing redundant programs)
# - inter_filter = logical; indicates whether programs should be filtered based on their similarity to programs of other cell lines
# - inter_min = minimum overlap with a program from another cell line 

# Returns a character vector with the NMF programs selected

robust_nmf_programs <- function(nmf_programs, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10) {
  
  # Select NMF programs based on the minimum overlap with other NMF programs from the same cell line
  intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y))))) 
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))             
  nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min]) 
  names(nmf_sel) <- names(nmf_programs)
  
  # Select NMF programs based on i) the maximum overlap with other NMF programs from the same cell line and
  # ii) the minimum overlap with programs from another cell line
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y)))) ## calculating intersection between all programs
  
  final_filter <- NULL 
  for(i in names(nmf_sel)) {
    a <- inter_intersect[grep(i, colnames(inter_intersect), invert = T),grep(i, colnames(inter_intersect))]
    b <- sort(apply(a, 2, max), decreasing = T) # for each cell line, ranks programs based on their maximum overlap with programs of other cell lines
    if(inter_filter==T) b <- b[b>=inter_min] # selects programs with a maximum intersection of at least 10
    if(length(b) > 1) {
      c <- names(b[1]) 
      for(y in 2:length(b)) {
        if(max(inter_intersect[c,names(b[y])]) <= intra_max) c <- c(c,names(b[y])) # selects programs iteratively from top-down. Only selects programs that have a intersection smaller than 10 with a previously selected programs
      }
      final_filter <- c(final_filter, c)
    } else {
      final_filter <- c(final_filter, names(b))
    }
  }
  return(final_filter)                                                      
}

