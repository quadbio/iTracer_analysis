plotFeature <- function (coord,
                         values = NULL,
                         emphasize = NULL,
                         col = NULL,
                         colorPal = NULL,
                         col_alpha = NULL,
                         balance_col = TRUE,
                         mask_col = "#efefef50",
                         value_ceiling = NULL,
                         value_floor = NULL,
                         pt_border = FALSE,
                         col_border = "#303030",
                         lwd_border = 0.2,
                         main = NA,
                         axis = FALSE,
                         xlab = "Dim-1", 
                         ylab = "Dim-2",
                         cex = 1,
                         cex.main = 1,
                         cex.lab = 1,
                         cex.axis = 1,
                         random_order = TRUE,
                         sort_by_value = FALSE,
                         edges = NULL,
                         lwd.edges = 0.2,
                         col.edges = "#bdbdbd50",
                         do_label = FALSE,
                         cex.label = 1,
                         label_min_cell = 10,
                         label_round = FALSE,
                         label_round_cex = cex.label * 3,
                         label_round_col = "#303030",
                         label_round_lwd = 0.5,
                         label_round_bg = "#fefefe50",
                         do_legend = F,
                         legend_pos = "right",
                         legend_cex = 1,
                         legend_pt_cex = legend_cex * 2,
                         ...) 
{
  if (! is.null(col)){ do_legend = F; do_label = F }
  if (is.numeric(values)){ do_legend = F; do_label = F }
  
  if (is.null(col)) {
    if (is.numeric(values)) {
      values <- c(value_ceiling, value_floor, values)
      
      if (is.null(colorPal)) {
        colorPal <- grDevices::colorRampPalette(c("darkgreen", 
                                                  "yellow", "red"))
      } else if (is.character(colorPal)) {
        colorPal <- colorRampPalette(colorPal)
      }
      
      if (balance_col & min(values, na.rm=T) < 0 & max(values, na.rm=T) > 0){
        values <- c(-max(abs(values)), max(abs(values)), values)
        cellColor <- adjustcolor(colorPal(30), alpha = 0.8)[as.numeric(cut(values, 
                                                                           breaks = 30, right = F, include.lowest = T))]
        values <- values[-(1:2)]
        cellColor <- cellColor[-(1:2)]
      } else{
        cellColor <- adjustcolor(colorPal(30), alpha = 0.8)[as.numeric(cut(values, 
                                                                           breaks = 30, right = F, include.lowest = T))]
        if (min(values, na.rm = T) == 0) 
          cellColor[values == 0] <- "#bdbdbd30"
      }
      
      values <- values[(length(c(value_ceiling, value_floor))+1):length(values)]
      cellColor <- cellColor[(length(c(value_ceiling, value_floor))+1):length(cellColor)]
    }
    else {
      if (is.character(values)) values <- as.factor(values)
      
      cols <- NULL
      if (!is.null(names(colorPal)) & sum(values %in% names(colorPal))>0){
        cols <- colorPal
      } else{
        if (is.null(colorPal)) colorPal <- scales::hue_pal()
        if (is.character(colorPal)) colorPal <- colorRampPalette(colorPal)
        if (is.function(colorPal)) cols <- colorPal(length(levels(values)))
        if (is.null(names(cols))) cols <- setNames(cols, levels(values))
      }
      cellColor <- cols[as.character(values)]
    }
    
    if (!is.null(col_alpha))
      cellColor <- adjustcolor(cellColor, col_alpha)
  }
  else {
    if (length(col) == 1) 
      col <- rep(col, nrow(coord))
    cellColor <- col
  }
  col_border <- rep(col_border, length(cellColor))
  col_border[is.na(cellColor)] <- NA
  
  if (axis){
    plot(coord, type = "n", main = main, xlab = xlab, ylab = ylab, 
         cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, 
         ...)
  } else{
    plot(coord, type = "n", main = main, bty = "n", xlab = NA, ylab = NA, xaxt = "n", yaxt = "n",
         cex.main = cex.main, ...)
  }
  if (!is.null(edges) & (is.matrix(edges) | is.data.frame(edges))) {
    for (i in 1:nrow(edges)) lines(coord[as.numeric(edges[i, ]), 1], coord[as.numeric(edges[i, ]), 2],
                                   lwd = lwd.edges, 
                                   col = col.edges)
  }
  
  if (sum(is.na(cellColor))>0){
    idx <- which(!is.na(cellColor))
    if (is.null(emphasize))
      emphasize <- idx
    emphasize <- intersect(emphasize, idx)
  }
  
  if (is.null(emphasize)) {
    idx_order <- seq(nrow(coord))
    if (random_order){
      idx_order <- sample(idx_order)
    } else if (sort_by_value){
      idx_order <- order(values)
    }
    
    if (pt_border){
      points(coord[idx_order,], col = col_border[idx_order], bg = cellColor[idx_order], pch = 21, cex = cex, lwd = lwd_border)
    } else{
      points(coord[idx_order,], col = cellColor[idx_order], pch = 16, cex = cex)
    }
  }
  else {
    points(coord, col = mask_col, pch = 16, cex = cex)
    idx_order <- emphasize
    if (length(idx_order) > 0){
      if (random_order){
        idx_order <- sample(idx_order)
      } else if (sort_by_value){
        idx_order <- idx_order[order(values[idx_order])]
      }
      
      if (pt_border){
        points(coord[idx_order, ], col = col_border[idx_order], bg = cellColor[idx_order], pch = 21, cex = cex, lwd = lwd_border)
      } else{
        points(coord[idx_order, ], col = cellColor[idx_order], pch = 16, cex = cex)
      }
    }
  }
  
  if (is.factor(values) & do_label){
    labs <- table(values)
    labs <- names(which(labs >= label_min_cell))
    labs_coord <- t(sapply(labs, function(x) colMeans(as.data.frame(coord)[which(values == x),], na.rm = T)))
    if (label_round) points(labs_coord, pch = 21, cex = label_round_cex, lwd = label_round_lwd, col = label_round_col, bg = label_round_bg)
    text(labs_coord, labels = labs, cex = cex.label)
  }
  if (do_legend){
    legend(legend_pos, legend = names(cols), pch = 16, col = cols, pt.cex = legend_pt_cex, cex = legend_cex, bty="n")
  }
}

plotMultiFeatures <- function(coord, mat_values, main = NULL, ..., ncol = NULL, byrow = T, mar = c(1,1,3,1), par_cex = NULL, auto_layout = T){
  if (nrow(mat_values) != nrow(coord) & ncol(mat_values) == nrow(coord))
    mat_values <- t(mat_values)
  if (nrow(mat_values) != nrow(coord))
    stop("The number of cells in the value matrix has to match with the coordinates.")
  if (is.null(main))
    main <- colnames(mat_values)
  if (length(main) == 0)
    main <- ""
  if (length(main) == 1)
    main <- rep(main, ncol(mat_values))
  
  if (is.null(ncol))
    ncol <- ceiling(sqrt(ncol(mat_values)))
  nrow <- ceiling(ncol(mat_values) / ncol)
  if (auto_layout)
    layout(matrix(1:(nrow*ncol), nrow = nrow, byrow = byrow))
  par(mar = mar)
  if (! is.null(par_cex))
    par(cex = par_cex)
  
  for(i in 1:ncol(mat_values))
    plotFeature(coord, mat_values[,i], main = main[i], ...)
}
