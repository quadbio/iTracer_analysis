# shuffle scar family labels for cells within the same barcode family
shuffle_scars <- function(clone_lab, # barcode family labels
                          scar_lab, # scar family labels
                          num_shuffle = 1000) # number of shuffling to generate
{
  res <- sapply(1:num_shuffle, function(i){
    idx <- which(!is.na(clone_lab))
    scars_clone <- setNames(split(data.frame(idx = idx, GeneBarcodeMerge = clone_lab[idx], ScarMerge = scar_lab[idx], stringsAsFactors=F),
                                  clone_lab[idx]),
                            NULL)
    scars_clone <- do.call(rbind, lapply(scars_clone, function(x){
      x$ScarMerge <- sample(x$ScarMerge)
      return(x)
    }))
    shuffled <- rep(NA, length(clone_lab))
    shuffled[scars_clone$idx] <- scars_clone$ScarMerge
    return(shuffled)
  })
  rownames(res) <- names(clone_lab)
  return(res)
}

# calculate scar composition differences between cell clusters
calculate_scar_prop_dist <- function(cluster_lab, # cell cluster labels
                                     clone_lab, # barcode family labels
                                     scar_lab, # scar family labels
                                     candidate_clones = NULL) # if not NULL, only barcode families given are considered
{
  cluster_lab <- factor(cluster_lab)
  clone_lab <- factor(clone_lab)
  if (is.null(candidate_clones))
    candidate_clones <- levels(clone_lab)
  
  scar_freq_cl_in_clones <- setNames(lapply(intersect(levels(clone_lab), candidate_clones), function(clone){
    idx <- which(clone_lab == clone & !is.na(scar_lab))
    freq <- table(scar_lab[idx], cluster_lab[idx])
    freq <- apply(freq+1, 2, function(x) x/sum(x))
    return(freq)
  }), intersect(levels(clone_lab), candidate_clones))
  
  dist <- as.matrix(dist(t(do.call(rbind, scar_freq_cl_in_clones))))
  rownames(dist) <- levels(cluster_lab)
  colnames(dist) <- levels(cluster_lab)
  return(dist)
}

# calculate scar composition differences and normalize based on shuffled scar family label
calculate_scar_prop_dist_and_test <- function(cluster_lab, # cell cluster labels
                                              clone_lab, # barcode family labels
                                              scar_lab, # scar family labels
                                              shuffled_scar_lab, # matrix of shuffled scar family labels
                                              candidate_clones = NULL, # if not NULL, only barcode families given are considered
                                              num_threads = 100) # number of threads to use
{
  obs_dist <- calculate_scar_prop_dist(cluster_lab, clone_lab, scar_lab, candidate_clones)
  library(doParallel)
  registerDoParallel(num_threads)
  shuffled_dists <- foreach(i = 1:ncol(shuffled_scar_lab), .combine = list, .multicombine = T, .maxcombine = ncol(shuffled_scar_lab)) %dopar%{
    shuffled_scars <- shuffled_scar_lab[,i]
    shuffled_dist <- calculate_scar_prop_dist(cluster_lab, clone_lab, shuffled_scars, candidate_clones)
  }
  q <- Reduce("+", lapply(lapply(shuffled_dists, "-", obs_dist), "<=", 0)) / length(shuffled_dists)
  zscores <- (obs_dist - apply(simplify2array(shuffled_dists), 1:2, mean, na.rm=T)) / apply(simplify2array(shuffled_dists), 1:2, sd, na.rm=T)
  zscores_norm <- (zscores - mean(as.numeric(zscores),na.rm=T))/sd(as.numeric(zscores),na.rm=T)
  res <- list(obs_dist = obs_dist, shuffled_dists = shuffled_dists, z = zscores, znorm = zscores_norm, q = q)
  return(res)
}
