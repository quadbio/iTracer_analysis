library(igraph)
library(scatterplot3d)

spots_trees <- lapply(1:4, function(i)
  data.frame(Idx_tree = i, read.csv(paste0("lightsheet/290520_KL_checked_tree",i,"_mamut.spot_df.csv"), stringsAsFactors=F, header=T, row.names=1), stringsAsFactors=F))
edges_trees <- lapply(1:4, function(i)
  read.csv(paste0("lightsheet/290520_KL_checked_tree",i,"_mamut.edge_df.csv"), stringsAsFactors=F, header=T, row.names=1))

id_replacement <- lapply(1:4, function(i){
  sizes_trees <- sapply(spots_trees, nrow)
  new_id <- setNames(seq(sum(sizes_trees[0:(i-1)])+1, sum(sizes_trees[0:i])), spots_trees[[i]]$ID)
  return(new_id)
})
for(i in 1:4){
  spots_trees[[i]]$ID <- id_replacement[[i]][as.character(spots_trees[[i]]$ID)]
  edges_trees[[i]]$SPOT_SOURCE_ID <- id_replacement[[i]][as.character(edges_trees[[i]]$SPOT_SOURCE_ID)]
  edges_trees[[i]]$SPOT_TARGET_ID <- id_replacement[[i]][as.character(edges_trees[[i]]$SPOT_TARGET_ID)]
}
spots <- do.call(rbind, spots_trees)
rownames(spots) <- spots$ID
edges <- do.call(rbind, edges_trees)
graph <- graph_from_data_frame(edges, directed=F)



# make dendrogram of tracked lineages
max_timeframe <- 130
layout(matrix(1:4, nrow=1))
for(idx_tree in 1:4){
  spots_tree <- spots[spots$Idx_tree == idx_tree & spots$FRAME %in% 0:max_timeframe,]
  edges_tree <- edges[edges$SPOT_SOURCE_ID %in% spots_tree$ID & edges$SPOT_TARGET_ID %in% spots_tree$ID,]
  graph_tree <- graph_from_data_frame(edges_tree, directed = T)
  degrees_tree <- degree(graph_tree, mode = "all")
  idx_leaves <- which(rownames(spots_tree) %in% names(which(degrees_tree == 1)) & spots_tree$FRAME != 0)
  
  merge_at <- sapply(idx_leaves, function(ix1){
    path1 <- names(shortest_paths(graph_tree, from = as.character(spots_tree$ID[spots_tree$FRAME == 0]), to = as.character(spots_tree$ID[ix1]))$vpath[[1]])
    sapply(idx_leaves, function(ix2){
      path2 <- names(shortest_paths(graph_tree, from = as.character(spots_tree$ID[spots_tree$FRAME == 0]), to = as.character(spots_tree$ID[ix2]))$vpath[[1]])
      merged_frame <- max(spots_tree[intersect(path1, path2),"FRAME"])
    })
  })
  diag(merge_at) <- NA
  rownames(merge_at) <- spots_tree$ID[idx_leaves]
  colnames(merge_at) <- spots_tree$ID[idx_leaves]
  hcl <- hclust(as.dist((max_timeframe-merge_at)/2))
  plot(hcl, main = paste0("Tree ",idx_tree," Progeny"), xlab = "", ylab = "Traced-back time", labels = FALSE, hang = -0.1)
}


# nuclei number along times
spots_t1 <- spots[spots$Idx_tree == 1,]
plot(seq(0,200)/2, as.numeric(table(spots_t1$FRAME[spots_t1$FRAME <= 200])), pch = 16, col = "#303030", cex = 1.5, bty="n", xlab = "Hours", ylab = "Nuclei number", cex.lab = 1.2, xaxt = "n")
axis(side = 1, at = seq(0,100, 25))
df <- data.frame(frame=seq(0,200)/2, count=as.numeric(table(spots_t1$FRAME[spots_t1$FRAME <= 200])))
m <- nls(count ~ 2^(frame * a), data = df, start = list(a = 0.1))
lines(seq(0,200)/2, predict(m, newdata = data.frame(frame = (0:200)/2)))



# make 3d scatter plot of L1 nuclei along the 100h of tracking
spots_t1 <- spots[spots$Idx_tree == 1 & spots$FRAME <= 200,]
edges_t1 <- edges[edges$SPOT_SOURCE_ID %in% spots_t1$ID & edges$SPOT_TARGET_ID %in% spots_t1$ID,]
idx <- which(rownames(spots_t1) %in% names(which(degree(graph) > 2)))
idx <- which(rownames(spots_t1) %in% as.character(edges_t1$SPOT_TARGET_ID[edges_t1$SPOT_SOURCE_ID %in% as.numeric(rownames(spots_t1)[idx])]))

col_ages <- colorRampPalette(c("#303030","#909090","#bdbdbd","#fbb4b9","#c51b8a","#7a0177"))(201)
p <- scatterplot3d(spots_t1[idx,c("POSITION_X","POSITION_Y","POSITION_Z")], axis=F, tick.marks = F, label.tick.marks = F, grid=T, mar=c(0,0,0,0), pch=16, type = "p",
                   color = col_ages[spots_t1$FRAME[idx]+1], cex.symbols = 1)
for(i in 1:nrow(edges_t1)){
  lines(p$xyz.convert(spots_t1[as.character(edges_t1[i,]),c("POSITION_X","POSITION_Y","POSITION_Z")]), lwd = 0.75,
        col = col_ages[spots_t1[as.character(edges_t1[i,1]),"FRAME"]+1])
}
points(p$xyz.convert(spots_t1[spots_t1$FRAME == 0,c("POSITION_X","POSITION_Y","POSITION_Z")]),
       pch = 16, col = col_ages[1], cex = 3)
points(p$xyz.convert(spots_t1[spots_t1$FRAME == 200,c("POSITION_X","POSITION_Y","POSITION_Z")]),
       pch = 16, col = col_ages[201], cex = 1.5)



# make 3d scatter plot of all tracked nuclei spots in 0-65h
frames <- 0:130
col_trees <- c("#da2b84","#fa9fb5","#4b1d9f","#17713f")

idx <- sample(which(spots$FRAME %in% frames))

p <- scatterplot3d(spots[idx,c("POSITION_X","POSITION_Y","POSITION_Z")], axis=F, tick.marks = F, label.tick.marks = F, grid=T, mar=c(0,0,0,0), pch=16, type = "h",
                   color = paste0(col_trees[spots$Idx_tree[idx]],"30"), cex.symbols = 0.5)
points(p$xyz.convert(spots[spots$FRAME == min(frames),c("POSITION_X","POSITION_Y","POSITION_Z")]),
       pch = 21, col = "#000000", bg = col_trees[spots$Idx_tree[spots$FRAME == min(frames)]], cex = 2, lwd = 1)
points(p$xyz.convert(spots[spots$FRAME == max(frames),c("POSITION_X","POSITION_Y","POSITION_Z")]),
       pch = 21, col = "#303030", bg = col_trees[spots$Idx_tree[spots$FRAME == max(frames)]], cex = 1, lwd = 0.5)



# make movie of tracked nuclei along the 65h of tracking
frames <- 0:130
col_trees <- c("#da2b84","#fa9fb5","#4b1d9f","#17713f")

library(scatterplot3d)
temppng <- paste0(sapply(frames, function(x) tempfile()), ".png")
output_file <- paste0(getwd(), "/movie_trees_frame0-130.gif")
for(i in seq(length(frames))){
  coord <- spots[spots$FRAME == frames[i],c("POSITION_X","POSITION_Y","POSITION_Z")]
  lims <- apply(spots[spots$FRAME %in% frames,c("POSITION_X","POSITION_Y","POSITION_Z")], 2, function(x) c(min(x),max(x)))
  pt_colors <- col_trees[spots$Idx_tree[spots$FRAME == frames[i]]]
  png(temppng[i], height=1000, width=1000)
  par(cex=1)
  scatterplot3d(coord, box=F, axis=F, tick.marks = F, label.tick.marks = F, grid=T, pch=16, mar=c(0,0,0,0),
                xlim = lims[,1], ylim=lims[,2], zlim=lims[,3],
                color = pt_colors, cex.symbols = 4.5)
  dev.off()
}
system(paste0("convert -delay ", 12, " -loop 0 ", paste(temppng, collapse = " "), " ", output_file))
file.remove(temppng)



# inter-nuclear distances
dists_spots <- as.matrix(dist(spots[,c("POSITION_X","POSITION_Y","POSITION_Z")]))

same_trees <- matrix(rep(spots$Idx_tree, each = nrow(spots)), nrow = nrow(spots))
same_trees <- same_trees == t(same_trees)
ventricle <- c(1, 1, 1, 2)
same_ventricles <- matrix(rep(ventricle[spots$Idx_tree], each = nrow(spots)), nrow = nrow(spots))
same_ventricles <- same_ventricles == t(same_ventricles)

frames <- 0:130
dists_same_trees <- lapply(frames, function(f){
  idx <- which(spots$FRAME == f)
  this_frame <- matrix(FALSE, nrow = nrow(spots), ncol = nrow(spots))
  this_frame[idx,idx] <- TRUE
  return(as.numeric(dists_spots[same_trees & this_frame]))
})
dists_same_ventricle <- lapply(frames, function(f){
  idx <- which(spots$FRAME == f)
  this_frame <- matrix(FALSE, nrow = nrow(spots), ncol = nrow(spots))
  this_frame[idx,idx] <- TRUE
  return(as.numeric(dists_spots[(same_ventricles & !same_trees) & this_frame]))
})
dists_diff_ventricle <- lapply(frames, function(f){
  idx <- which(spots$FRAME == f)
  this_frame <- matrix(FALSE, nrow = nrow(spots), ncol = nrow(spots))
  this_frame[idx,idx] <- TRUE
  return(as.numeric(dists_spots[(!same_ventricles) & this_frame]))
})

plot_scatter_dists <- function(dists, ylim = NULL, ...){
  medians <- sapply(dists, median)
  uppers <- sapply(dists, quantile, 0.99)
  lowers <- sapply(dists, quantile, 0.01)
  if (is.null(ylim))
    ylim <- c(min(lowers), max(uppers))
  plot(seq(length(medians)), medians, type = "n", bty = "n", ylim = ylim, ...)
  polygon(c(seq(length(medians)), rev(seq(length(medians)))), c(uppers, rev(lowers)), border = NA, col = "#bdbdbd90")
  points(seq(length(medians)), medians, pch = 21, col = "#303030", bg = "#696969", cex = 1.5, lwd = 0.5)
  lines(predict(smooth.spline(seq(length(medians)), medians, df = 3), seq(length(medians))), col = "#000000")
}
layout(matrix(1:3, nrow=1)); par(mar=c(5,5,3,1))
plot_scatter_dists(dists_same_trees, ylim = c(0,450), xlab = "Hours", xaxt = "n", ylab = "Distances between nuclei", main = "Same trees", cex.main = 1.5, cex.lab = 1.2)
axis(side = 1, at = seq(0,120,20), labels = seq(0,60,10))
plot_scatter_dists(dists_same_ventricle, ylim = c(0,450), xlab = "Hours", xaxt = "n", ylab = "Distances between nuclei", main = "Different trees at the same lumen", cex.main = 1.5, cex.lab = 1.2)
axis(side = 1, at = seq(0,120,20), labels = seq(0,60,10))
plot_scatter_dists(dists_diff_ventricle, ylim = c(0,450), xlab = "Hours", xaxt = "n", ylab = "Distances between nuclei", main = "Different lumens", cex.main = 1.5, cex.lab = 1.2)
axis(side = 1, at = seq(0,120,20), labels = seq(0,60,10))
