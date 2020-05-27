# shuffle barcode family labels for cells of the same organoid
shuffle_barcodes <- function(org_lab, # organoid labels
                             clone_lab, # barcode family labels
                             num_shuffle = 100) # number of shuffling to generate
{
  res <- sapply(1:num_shuffle, function(i){
    clones_org <- split(data.frame(idx = 1:length(org_lab), organoid = org_lab, GeneBarcodeMerge = clone_lab), org_lab)
    clones_shuffled <- do.call(rbind, lapply(clones_org, function(clones){
      idx_barcoded <- which(!is.na(clones$GeneBarcodeMerge))
      clones$GeneBarcodeMerge[idx_barcoded] <- sample(clones$GeneBarcodeMerge[idx_barcoded])
      return(clones)
    }))
    clones_shuffled <- clones_shuffled[order(clones_shuffled$idx),]
    return(clones_shuffled$GeneBarcodeMerge)
  })
  rownames(res) <- names(clone_lab)
  return(res)
}

# calculate barcode family composition similarity between cell clusters
count_clone_linkage <- function(cluster_lab, # cell cluster labels
                                clone_lab, # barcode family labels
                                candidate_clones = NULL) # if not NULL, only barcode families given are considered
{
  cluster_lab <- factor(cluster_lab)
  clone_lab <- factor(clone_lab)
  if (is.null(candidate_clones))
    candidate_clones <- levels(clone_lab)
  
  sapply(levels(cluster_lab), function(cl1){
    clone_cl1 <- intersect(candidate_clones, sort(unique(clone_lab[cluster_lab == cl1])))
    freq_clones_cl1 <- table(clone_lab[cluster_lab == cl1])[clone_cl1]
    num_cell_clone <- sum(clone_lab %in% candidate_clones & cluster_lab == cl1, na.rm=T)
    
    sapply(levels(cluster_lab), function(cl2){
      clone_cl2 <- intersect(candidate_clones, sort(unique(clone_lab[cluster_lab == cl2])))
      freq_clones_cl2 <- table(clone_lab[cluster_lab == cl2])[clone_cl2]
      shared_clones <- intersect(clone_cl1, clone_cl2)
      num_links <- sum(freq_clones_cl1[shared_clones] * freq_clones_cl2[shared_clones])
      if (cl1 == cl2)
        num_links <- num_links - num_cell_clone
      return(num_links)
    })
  })
}

# calculate barcode family composition similarity and normalize to random similarities based on shuffled labels
count_clone_linkage_and_test <- function(cluster_lab, # cell cluster labels
                                         clone_lab, # barcode family labels
                                         shuffled_clone_lab, # matrix of shuffled barcode family labels
                                         candidate_clones = NULL, # if not NULL, only barcode families given are considered
                                         num_threads = 100)
{
  obs_linkage <- count_clone_linkage(cluster_lab, clone_lab, candidate_clones)
  
  library(doParallel)
  registerDoParallel(num_threads)
  shuffled_linkages <- foreach(i = 1:ncol(shuffled_clone_lab), .combine = list, .multicombine = T, .maxcombine = ncol(shuffled_clone_lab)) %dopar%{
    shuffled_clone <- shuffled_clone_lab[,i]
    shuffled_linkage <- count_clone_linkage(cluster_lab, shuffled_clone, candidate_clones)
  }
  q <- Reduce("+", lapply(lapply(lapply(shuffled_linkages, "-", obs_linkage), sign), ">", 0))/length(shuffled_linkages)
  zscores <- (obs_linkage - apply(simplify2array(shuffled_linkages), 1:2, mean, na.rm=T)) / apply(simplify2array(shuffled_linkages), 1:2, sd, na.rm=T)
  zscores_norm <- (zscores - mean(as.numeric(zscores), na.rm=T))/sd(as.numeric(zscores), na.rm=T)
  res <- list(obs_link = obs_linkage, shuffled_links = shuffled_linkages, z = zscores, znorm = zscores_norm, q = q)
  return(res)
}
