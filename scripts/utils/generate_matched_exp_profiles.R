# Script from Avishay Spitzer
# mdata - the metadata table that contains for each cell the CellID, Patient, Timepoint, Sample, State etc
# umi_data_all - a named list per sample of UMI count matrices

euclidean_dist <- function(x, y) {
  sqrt(sum((x - y)^2))
}

# Helper function that normalizes UMI counts to per-million
umi2upm <- function(m) {
  count_sum <- colSums(m)
  upm_data <- (t(t(m)/count_sum)) * 1e+06
  upm_data
}

generate_matched_exp_profiles <- function(meta_data, exp_data, patients, states, min_cells = 25, min_exp = 1, genes_subset = NULL, valid_genes = NULL) {
  
  meta_data <- meta_data %>%
    filter(State %in% states, Patient %in% patients)
  
  print("Generating groups")
  
  start_time <- Sys.time()
  groups <- lapply(patients, function(pt) {
    d1 <- meta_data %>%
      filter(Patient == pt, Timepoint == "T1")
    dim(d1)
    
    d2 <- meta_data %>%
      filter(Patient == pt, Timepoint == "T2")
    dim(d2)
    
    d <- rbind(d1, d2)
    dim(d)
    
    d_stats <- d %>%
      group_by(Timepoint, State) %>%
      summarise(n = n())
    
    d_stats <- d_stats %>%
      group_by(State) %>%
      filter(min(n) >= min_cells) %>%
      ungroup()
    
    d_stats <- d_stats %>%
      group_by(State) %>%
      filter("T1" %in% Timepoint & "T2" %in% Timepoint)
    
    if(!("T1" %in% d_stats$Timepoint) | !("T2" %in% d_stats$Timepoint))
      return(NULL)
    
    d <- d %>%
      filter(State %in% d_stats$State)
    table(d$State, d$Timepoint)
    
    d <- d %>%
      group_by(Timepoint, State) %>%
      sample_n(min_cells) %>%
      ungroup()
    table(d$State, d$Timepoint)
    
    d <- d %>%
      group_by(State) %>%
      filter("T1" %in% Timepoint & "T2" %in% Timepoint)
    
    return(d)
  })
  end_time <- Sys.time()
  end_time - start_time
  
  groups <- do.call(rbind, groups)
  
  if(is.null(groups))
    return(NULL)
  
  groups_stats <- groups %>%
    group_by(Sample, Patient, Timepoint, State) %>%
    summarise(n = n())
  
  if(is.null(valid_genes))
    valid_genes <- rownames(exp_data[[1]])
  
  print("Generating profiles")
  
  start_time <- Sys.time()
  profiles <- lapply(patients, function(pt) {
    
    d <- groups %>%
      filter(Patient == pt)
    
    if(nrow(d) == 0) {
      print(paste0("No group found for patient ", pt))
      
      return(NULL)
    }
    
    s1 <- unique(d$Sample[d$Timepoint == "T1"])
    s2 <- unique(d$Sample[d$Timepoint == "T2"])
    
    # print(s1); print(s2)
    
    m1 <- exp_data[[s1]]
    m2 <- exp_data[[s2]]
    
    m1 <- m1[valid_genes, d$CellID[d$Timepoint == "T1"]]
    m2 <- m2[valid_genes, d$CellID[d$Timepoint == "T2"]]
    
    m1 <- umi2upm(m1)
    m2 <- umi2upm(m2)
    
    m1 <- rowMeans(m1)
    m2 <- rowMeans(m2)
    
    m1 <- log2(m1 + 1)
    m2 <- log2(m2 + 1)
    
    genes <- genes_subset
    
    if(is.null(genes_subset)) {
      m <- (m1 + m2) / 2
      genes <- names(m[m >= min_exp])
    
      print(paste0("Detected ", length(genes), " genes with log2(mean exp) > ", min_exp, " in both timepoints"))
    }
    
    m1 <- m1[genes]
    m2 <- m2[genes]
    
    diff <- m2 - m1
    
    res <- tibble(Patient = pt,
                  Gene = names(diff),
                  T1 = m1,
                  T2 = m2,
                  Diff = diff)
    
    return(res)
  })
  end_time <- Sys.time()
  end_time - start_time
  
  lengths(profiles)
  
  profiles <- do.call(rbind, profiles)
  dim(profiles)
  
  res <- list(groups = groups, group_stats = groups_stats, profiles = profiles)
  
  return(res)
}
