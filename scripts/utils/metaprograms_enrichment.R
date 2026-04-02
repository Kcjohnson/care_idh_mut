# Script to determine the overlap between a set of metaprograms and previously published genesets/metaprograms

metaprograms_enrichment <- function(mp_list, pathways = NULL, min_gs_size = 40, max_gs_size = 500) {
  
  if(is.null(pathways)) {
    
    stop("Error: Pathway list or dataframe must be provided.")
    
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
