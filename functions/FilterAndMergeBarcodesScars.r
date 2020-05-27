library(stringdist)
library(igraph)

# filter barcode
do_filter_barcodes <- function(df_barcode) # table output by ExtractRFPbarcodesFromBam.pl
{
  # filter 1: read coverage
  count <- melt(table(df_barcode$Reads))
  count <- merge(data.frame(count=3:max(count$Var1)),count, by.y = "Var1",by.x = "count",all = TRUE)
  count[is.na(count)] <- 0
  count$fit <- predict(loess(value~log(count), count, span = 0.1), count)
  count$fit[count$fit<1] <- 1
  turn <- count$count[min(which(count$fit[-1] - count$fit[-length(count$fit)]>0))]
  idx_high_coverage <- which(df_barcode$Reads >= turn)
  df_barcode <- df_barcode[idx_high_coverage,]
  
  # filter 2: start and end with either A or T
  idx_reasonable_seq <- which(substring(df_barcode$GeneBarcode, 1, 1) %in% c("A","T") & substring(df_barcode$GeneBarcode, nchar(df_barcode$GeneBarcode), nchar(df_barcode$GeneBarcode)) %in% c("A","T"))
  df_barcode <- df_barcode[idx_reasonable_seq,]
  
  # filter 3: if there are molecules where the same UMI for different barcodes is detected in a cell, only remain the one with higher coverage
  df_barcode_cells <- split(data.frame(idx = 1:nrow(df_barcode), df_barcode, stringsAsFactors=F), df_barcode$Cell)
  idx_max_coverage_per_umi <- sort(unlist(lapply(df_barcode_cells, function(df){
    df_each_umi <- split(df, df$UMI)
    return(unlist(lapply(df_each_umi, function(df) df$idx[which.max(df$Reads)])))
  })))
  df_barcode <- df_barcode[idx_max_coverage_per_umi,]
  
  # filter 4: filter out molecules where the same UMI and barcode is detected across multiple cells
  df_barcode_barcode_umi <- split(data.frame(idx = 1:nrow(df_barcode), df_barcode, stringsAsFactors=F), paste0(df_barcode$UMI, "_", df_barcode$GeneBarcode))
  df_barcode_barcode_umi <- df_barcode_barcode_umi[sapply(df_barcode_barcode_umi, nrow) == 1]
  idx_once_barcode_umi <- sort(unlist(lapply(df_barcode_barcode_umi, function(x) x$idx)))
  df_barcode <- df_barcode[idx_once_barcode_umi,]
  
  # filter 5: filter out molecules where the hamming distance between detected barcodes within a cell is 1, to pick the barcode that are covered by more UMI or that have a higher coverage
  df_barcode_cells <- split(data.frame(idx = 1:nrow(df_barcode), df_barcode, stringsAsFactors=F), df_barcode$Cell)
  idx_diff_enough <- sort(unlist(lapply(df_barcode_cells, function(df){
    barcodes <- sort(unique(df$GeneBarcode))
    umi_barcodes <- table(df$GeneBarcode)[barcodes]
    read_barcodes <- sapply(split(df[,c("Reads","GeneBarcode")], df$GeneBarcode), function(x) sum(x$Reads))[barcodes]
    
    hamming_dist <- stringdist::stringdistmatrix(barcodes, barcodes, method="hamming")
    hamming_adj <- hamming_dist <= 1
    rownames(hamming_adj) <- barcodes
    colnames(hamming_adj) <- barcodes
    graph <- igraph::graph_from_adjacency_matrix(hamming_adj, mode = "undirected")
    components <- setNames(igraph::components(graph)$membership, names(igraph::V(graph)))
    valid_barcodes <- sapply(tapply(barcodes, components[barcodes], list), function(x) x[order(umi_barcodes[x], read_barcodes[x], decreasing=T)[1]] )
    
    return(df$idx[df$GeneBarcode %in% valid_barcodes])
  })))
  df_barcode <- df_barcode[idx_diff_enough,]
  
  return(df_barcode)
}

# filter scars
do_filter_scars <- function(df_scar) # table output by ExtractScarsFromBam.pl
{
  # filter 1: read coverage
  count <- melt(table(df_scar$Reads))
  count <- merge(data.frame(count=3:max(count$Var1)),count, by.y = "Var1",by.x = "count",all = TRUE)
  count[is.na(count)] <- 0
  count$fit <- predict(loess(value~log(count), count, span = 0.1), count)
  count$fit[count$fit<1] <- 1
  turn <- count$count[min(which(count$fit[-1] - count$fit[-length(count$fit)]>0))]
  idx_high_coverage <- which(df_scar$Reads >= turn)
  df_scar <- df_scar[idx_high_coverage,]
  
  # filter 2: if there are molecules where the same UMI for different scars is detected in a cell, only remain the one with higher coverage
  df_scar_cells <- split(data.frame(idx = 1:nrow(df_scar), df_scar, stringsAsFactors=F), df_scar$Cell)
  idx_max_coverage_per_umi <- sort(unlist(lapply(df_scar_cells, function(df){
    df_each_umi <- split(df, df$UMI)
    return(unlist(lapply(df_each_umi, function(df) df$idx[which.max(df$Reads)])))
  })))
  df_scar <- df_scar[idx_max_coverage_per_umi,]
  
  # filter 3: filter out molecules where the same UMI and scar is detected across multiple cells
  df_scar_scar_umi <- split(data.frame(idx = 1:nrow(df_scar), df_scar, stringsAsFactors=F), paste0(df_scar$UMI, "_", df_scar$Scar))
  df_scar_scar_umi <- df_scar_scar_umi[sapply(df_scar_scar_umi, nrow) == 1]
  idx_once_barcode_umi <- sort(unlist(lapply(df_scar_scar_umi, function(x) x$idx)))
  df_scar <- df_scar[idx_once_barcode_umi,]
  
  return(df_scar)
}

# only consider UMIs with both barcode and scar detected
merge_barcode_scars <- function(df_barcode, # filtered barcode table
                                df_scar) # filtered scar table
{
  # only shared cells are considered
  cells <- intersect(df_barcode$Cell, df_scar$Cell)
  df_barcode <- df_barcode[df_barcode$Cell %in% cells,]
  df_scar <- df_scar[df_scar$Cell %in% cells,]
  
  # only shared UMIs in the same cells are considered
  df_barcode_cells <- split(df_barcode, df_barcode$Cell)
  df_scar_cells <- split(df_scar, df_scar$Cell)
  df_barcode_scar <- do.call(rbind, lapply(cells, function(cell){
    df_barcode_cell <- df_barcode_cells[[cell]]
    df_scar_cell <- df_scar_cells[[cell]]
    umis <- intersect(df_barcode_cell$UMI, df_scar_cell$UMI)
    
    barcode_family <- setNames(df_barcode_cell$GeneBarcode, df_barcode_cell$UMI)[umis]
    scar_family <- setNames(df_scar_cell$Scar, df_scar_cell$UMI)[umis]
    orders <- order(barcode_family)
    df <- unique(data.frame(barcode = barcode_family[orders], scar = scar_family[orders], stringsAsFactors = F))
    return(data.frame(cell = cell,
                      GeneBarcodeMerged = paste(df$barcode, collapse="_"),
                      ScarMerged = paste(df$scar, collapse="_"), stringsAsFactors=F))
  }))
  df_barcode_scar <- df_barcode_scar[which(df_barcode_scar$GeneBarcodeMerged!=""),]
  return(df_barcode_scar)
}
