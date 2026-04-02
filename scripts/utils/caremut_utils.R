###################################################################
#### CAREmut study utils - scoring cells and generating null distributions 
###################################################################
generate_null_dist <- function(umi_data_list, md, sigs, genes_subset = NULL, n_iter = 20, n_cells = 10000, verbose = F) {
  
  if(verbose == T)
    print("Generating NULL distribution for state/cell-type classification")
  
  start_time <- Sys.time()
  permuted_data <- lapply(1:n_iter, function(i) {
    
    if(verbose == T)
      print(paste0("Iteration ", i))
    
    cellids <- md %>% sample_n(n_cells) %>% pull(CellID)
    
    if(!is.null(genes_subset))
      genes <- genes_subset
    else
      genes <- rownames(umi_data_list[[1]])
    
    if(verbose == T)
      print("Aggregating UMI data")
    
    umi_data <- lapply(names(umi_data_list), function(sname) {
      d <- md %>% filter(SampleID == sname, CellID %in% cellids)
      if(nrow(d) < 2)
        return(NULL)
      m <- umi_data_all[[sname]]
      m <- m[genes, d$CellID]
      m
    })
    umi_data <- do.call(cbind, umi_data)
    gc()
    
    cellids <- cellids[cellids %in% colnames(umi_data)]
    
    m <- umi_data[, cellids]
    dim(m)
    
    if(verbose == T)
      print("Permuting matrix")
    
    m <- t(apply(m, 1, gtools::permute))
    colnames(m) <- cellids; rownames(m) <- genes
    
    if(verbose == T)
      print("Creating object")
    
    tmp <- CreateSeuratObject(counts = m,
                              project = paste0("Perm", i),
                              min.cells = 3, min.features = 200)
    
    if(verbose == T)
      print("Normalizing")
    
    tmp <- NormalizeData(tmp, normalization.method = "LogNormalize", scale.factor = 10000)
    
    tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
    
    if(verbose == T)
      print("Scaling")
    
    all.genes <- rownames(tmp)
    tmp <- ScaleData(tmp, features = all.genes)
    
    if(verbose == T)
      print("Scoring")
    
    tmp <- AddModuleScore(object = tmp, features = sigs)
    
    d <- as_tibble(tmp@meta.data)
    colnames(d)[grep("Cluster", colnames(d))] <- names(sigs)
    
    d[, -c(1:3)]
  })
  end_time <- Sys.time()
  end_time - start_time
  
  permuted_data <- do.call(rbind, permuted_data)
  
  return(permuted_data)
}

score_within_samples_caremut <- function(umi_data_list, md, sigs) {
  
  samples <- unique(md$SampleID)
  
  res_all <- lapply(samples, function(sname) {
    
    print(sname)
    
    d <- md %>% dplyr::filter(SampleID == sname)
    
    ud <- umi_data_list[[sname]][, d$CellID]
    
    npcs <- 100
    
    if(nrow(d) < npcs) {
      npcs <- nrow(d) - 1
      print(paste0("npcs set to ", npcs))    
    }
    

    tmp <- CreateSeuratObject(counts = ud,
                              project = sname,
                              min.cells = 0, min.features = 0)
    
    print("Normalizing")
    
    tmp <- NormalizeData(tmp, normalization.method = "LogNormalize", scale.factor = 10000)
    
    # I don't think these steps are necessary for AddModuleScore. 
    # Leaving in since they are used in the generate_null_dist() function.
    tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
    
    print("Scaling")
    
    all.genes <- rownames(tmp)
    
    tmp <- ScaleData(tmp, features = all.genes)
    
    tmp <- AddModuleScore(object = tmp, features = sigs)  
    
    scores <- tmp@meta.data[, grep("Cluster", colnames(tmp@meta.data))]
    colnames(scores) <- names(sigs)
    
    scores <- as_tibble(scores, rownames = "CellID")
    
    res <- d %>%
      select(CellID, SampleID) %>%
      left_join(scores, by = "CellID")
    
    return(res)
  })
  res_all <- do.call(rbind, res_all)
  dim(res_all)
  
  res_all <- as_tibble(res_all)
  
  return(res_all)
}

################################################################################
# This is a slightly modified function that removes the junk genes
score_within_samples_no_junk_hierarchy <- function(umi_data_list, md, sigs) {
  
  samples <- unique(md$SampleID)
  
  res_all <- lapply(samples, function(sname) {
    
    print(sname)
    
    d <- md %>% dplyr::filter(SampleID == sname)
    
    ud <- umi_data_list[[sname]][, d$CellID]
    
    # Identify and remove "junk genes", which may impact the scoring function
    junk_genes <- c(rownames(ud)[grep("\\.", rownames(ud))],
                    rownames(ud)[grep("-AS*", rownames(ud))],
                    rownames(ud)[grep("LINC", rownames(ud))],
                    rownames(ud)[grep("^RP[S|L]", rownames(ud))],
                    rownames(ud)[grep("^MT-", rownames(ud))])
    
    print(paste0("Found ", length(junk_genes), " junk genes"))
    
    valid_genes <- rownames(ud)[!rownames(ud) %in% junk_genes]
    
    print(paste0("Keeping ", length(valid_genes), " valid genes"))
    
    ud <- ud[valid_genes, ]
    
    npcs <- 100
    
    if(nrow(d) < npcs) {
      npcs <- nrow(d) - 1
      print(paste0("npcs set to ", npcs))    
    }
    
    tmp <- CreateSeuratObject(counts = ud,
                              project = sname,
                              min.cells = 0, min.features = 0)
    
    print("Normalizing")
    
    tmp <- NormalizeData(tmp, normalization.method = "LogNormalize", scale.factor = 10000)
    
    tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 2000)
    
    print("Scaling")
    
    all.genes <- rownames(tmp)
    
    tmp <- ScaleData(tmp, features = all.genes)
    
    tmp <- AddModuleScore(object = tmp, 
                          features = sigs,
                          seed = 1)  
    
    scores <- tmp@meta.data[, grep("Cluster", colnames(tmp@meta.data))]
    colnames(scores) <- names(sigs)
    
    scores <- as_tibble(scores, rownames = "CellID")
    
    res <- d %>%
      select(CellID, SampleID) %>%
      left_join(scores, by = "CellID")
    
    return(res)
  })
  res_all <- do.call(rbind, res_all)
  dim(res_all)
  
  res_all <- as_tibble(res_all)
  
  return(res_all)
}


################################################################################
# Helper function that normalizes UMI counts to per-million
umi2upm <- function(m) {
  count_sum <- colSums(m)
  upm_data <- (t(t(m)/count_sum)) * 1e+06
  upm_data
}
################################################################################

################################################################################
# Function to score samples for the IDH-mutant hierarchy plot using log2(cpm/10+1)
################################################################################
score_within_samples_log2_hierarchy <- function(umi_data_list, md, sigs) {
  
  samples <- unique(md$SampleID)
  
  res_all <- lapply(samples, function(sname) {
    
    print(sname)
    
    d <- md %>% filter(SampleID == sname)
    
    ud <- umi_data_list[[sname]][, d$CellID]
    
    # Identify and remove "junk genes", which may impact the scoring function
    junk_genes <- c(rownames(ud)[grep("\\.", rownames(ud))],
                    rownames(ud)[grep("-AS*", rownames(ud))],
                    rownames(ud)[grep("LINC", rownames(ud))],
                    rownames(ud)[grep("^RP[S|L]", rownames(ud))],
                    rownames(ud)[grep("^MT-", rownames(ud))])
    
    print(paste0("Found ", length(junk_genes), " junk genes"))
    
    valid_genes <- rownames(ud)[!rownames(ud) %in% junk_genes]
    
    print(paste0("Keeping ", length(valid_genes), " valid genes"))
    
    ud <- ud[valid_genes, ]
    
    # Calculating row means to keep the highest expressing genes.
    # Blocking this off for now because it removes some of the signature genes for Tirosh/Venteicher
    #rm <- log2(rowMeans(umi2upm(ud)) + 1)
    
    #rm <- rm[rm > 4]
    
    #genes <- names(rm)
    
    #print(paste0("Found ", length(genes), " HE genes"))
    
    #ud <- ud[genes, ]
    
    print(paste0("Converting UMI counts to UMI counts per million for", sname))
    m <- umi2upm(ud)
    
    print("log2-transforming UMI per 100K")
    m <- as.matrix(log2(m / 10 + 1))
    
    tmp <- CreateSeuratObject(counts = m, project = sname, min.cells = 0, min.features = 0, verbose = T)
    
    print("Scoring with log2 transformed data with Seurat's AddModuleScore")
    tmp <- AddModuleScore(object = tmp, 
                          features = sigs,
                          seed = 1)  
    
    scores <- tmp@meta.data[, grep("Cluster", colnames(tmp@meta.data))]
    colnames(scores) <- names(sigs)
    
    scores <- as_tibble(scores, rownames = "CellID")
    
    res <- d %>%
      select(CellID, SampleID) %>%
      left_join(scores, by = "CellID")
    
    return(res)
  })
  res_all <- do.call(rbind, res_all)
  dim(res_all)
  
  res_all <- as_tibble(res_all)
  
  return(res_all)
}


################################################################################
# Function to score samples for the IDH-mutant hierarchy plot - Modified to restrict to highly expressed genes
################################################################################
score_within_samples_log2_hierarchy_he <- function(umi_data_list, md, sigs) {
  
  samples <- unique(md$SampleID)
  
  res_all <- lapply(samples, function(sname) {
    
    print(sname)
    
    d <- md %>% filter(SampleID == sname)
    
    ud <- umi_data_list[[sname]][, d$CellID]
    
    # Identify and remove "junk genes", which may impact the scoring function
    junk_genes <- c(rownames(ud)[grep("\\.", rownames(ud))],
                    rownames(ud)[grep("-AS*", rownames(ud))],
                    rownames(ud)[grep("LINC", rownames(ud))],
                    rownames(ud)[grep("^RP[S|L]", rownames(ud))],
                    rownames(ud)[grep("^MT-", rownames(ud))])
    
    print(paste0("Found ", length(junk_genes), " junk genes"))
    
    valid_genes <- rownames(ud)[!rownames(ud) %in% junk_genes]
    
    print(paste0("Keeping ", length(valid_genes), " valid genes"))
    
    ud <- ud[valid_genes, ]
    
    # Calculating row means to keep the highest expressing genes.
    # Blocking this off for now because it removes some of the signature genes for Tirosh/Venteicher
    rm <- log2(rowMeans(umi2upm(ud)) + 1)
    
    rm <- rm[rm > 4]
    
    genes <- names(rm)
    
    print(paste0("Found ", length(genes), " HE genes"))
    
    ud <- ud[genes, ]
    
    print(paste0("Converting UMI counts to UMI counts per million for", sname))
    m <- umi2upm(ud)
    
    print("log2-transforming UMI per 100K")
    m <- as.matrix(log2(m / 10 + 1))
    
    tmp <- CreateSeuratObject(counts = m, project = sname, min.cells = 0, min.features = 0, verbose = T)
    
    print("Scoring with log2 transformed data with Seurat's AddModuleScore")
    tmp <- AddModuleScore(object = tmp, features = sigs)  
    
    scores <- tmp@meta.data[, grep("Cluster", colnames(tmp@meta.data))]
    colnames(scores) <- names(sigs)
    
    scores <- as_tibble(scores, rownames = "CellID")
    
    res <- d %>%
      select(CellID, SampleID) %>%
      left_join(scores, by = "CellID")
    
    return(res)
  })
  res_all <- do.call(rbind, res_all)
  dim(res_all)
  
  res_all <- as_tibble(res_all)
  
  return(res_all)
}


####### From scrabble: Quadrant hierarchy plot for the IDHwt single cell analyses 
hierarchy = function(m, quadrants = NULL, log.scale = T) {
  
  if (!is.null(quadrants)) {
    stopifnot(all(unlist(quadrants) %in% colnames(m)))
    dat = as.data.frame(sapply(quadrants, function(col) do.call(pmax, list(as.data.frame(m[, col])))))
  } else {
    stopifnot(ncol(m) == 4)
    dat = as.data.frame(m)
  }
  
  rows = rownames(m)
  colnames(dat) = c('bl', 'br', 'tl', 'tr')
  
  dat = dat %>%
    dplyr::mutate(bottom = pmax(bl, br),
                  top = pmax(tl, tr),
                  b.center = br - bl,
                  t.center = tr - tl,
                  x = ifelse(bottom > top, b.center, t.center), # dependent var
                  x.scaled = (sign(x) * log2(abs(x) + 1)),
                  y = top - bottom, # independent var
                  y.scaled = (sign(y) * log2(abs(y) + 1)))
  
  if (!log.scale) dat = dplyr::transmute(dat, X = x, Y = y)
  else dat = dplyr::transmute(dat, X = x.scaled, Y = y.scaled)
  rownames(dat) = rows
  class(dat) = append(class(dat), 'hierarchy')
  dat
}

## Butterfly plot
plot_hierarchy = function(X,
                          quadrant.names = c('bl', 'br', 'tl', 'tr'),
                          main = NULL,
                          xlab = 'Relative meta-module score [log2(|SC1-SC2|+1)]',
                          ylab = 'Relative meta-module score [log2(|SC1-SC2|+1)]',
                          groups = NULL,
                          group.cols = NULL, 
                          legend = T,
                          legend.pos = 'bottom',
                          legend.horiz = T) {
  
  if (is.null(groups)) col = 'darkred'
  else col = 'grey85'
  
  plot(X[,1], X[,2], pch = 20, col = col, main = main, xlab = xlab, ylab = ylab)
  
  if (is.null(groups)) legend = F
  else {
    stopifnot(!is.null(names(groups)))
    stopifnot(all(groups %in% rownames(X)))
    groups = split(groups, names(groups))
    Xgrp = sapply(groups, function(rows) X[rows,,drop = F], simplify = F)
    if (!is.null(group.cols)) colgrp = group.cols[names(groups)]
    else colgrp = rainbow(n = length(Xgrp))
    Map(points,
        x = sapply(Xgrp, `[[`, 1, simplify = F),
        y = sapply(Xgrp, `[[`, 2, simplify = F),
        col = colgrp,
        MoreArgs = list(pch = 20))
  }
  
  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)
  
  if (legend) {
    legend(legend.pos,
           fill = colgrp,
           legend = names(groups),
           horiz = legend.horiz,
           cex = 0.8,
           box.col = 'white',
           bg = 'white',
           box.lwd = 0)
  }
  
  Names = quadrant.names
  cex = 1.2
  mtext(side = 1, adj = 0, text = Names[1], cex = cex, line = cex - 1)
  mtext(side = 1, adj = 1, text = Names[2], cex = cex, line = cex - 1)
  mtext(side = 3, adj = 0, text = Names[3], cex = cex)
  mtext(side = 3, adj = 1, text = Names[4], cex = cex)
}


## To plot density of points
# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#### END ####