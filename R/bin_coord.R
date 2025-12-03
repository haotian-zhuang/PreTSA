#' Optional spatial binning
#'
#' @param expr The normalized gene expression matrix. Rows represent genes and columns represent spots/cells
#' @param coord The matrix of spatial locations. Rows represent spots/cells. The column named "row" ("col") represents row (column) coordinates
#' @param num_bin Number of bins in each spatial dimension (100 by default)
#'
#' @returns A list of binned `expr` and `coord`
#' @export
#'
bin_coord <- function(expr, coord, num_bin = 100) {
  
  rbreak <- seq(from = min(coord[, "row"]), to = max(coord[, "row"]), length.out = num_bin + 1)
  cbreak <- seq(from = min(coord[, "col"]), to = max(coord[, "col"]), length.out = num_bin + 1)
  
  rbin <- cut(x = coord[, "row"], breaks = rbreak, labels = FALSE, include.lowest = TRUE, right = TRUE)
  cbin <- cut(x = coord[, "col"], breaks = cbreak, labels = FALSE, include.lowest = TRUE, right = TRUE)
  
  row_center <- (rbreak[-1] + rbreak[-length(rbreak)]) / 2
  col_center <- (cbreak[-1] + cbreak[-length(cbreak)]) / 2
  
  coord_agg <- expand.grid(row = row_center, col = col_center) # The first factors vary fastest
  rownames(coord_agg) <- paste0("Bin:", seq_len(nrow(coord_agg)))
  coord_agg <- as.matrix(coord_agg)
  
  bin_id <- (cbin - 1L) * num_bin + rbin
  bin_id <- paste0("Bin:", bin_id)
  names(bin_id) <- rownames(coord)
  bin_id <- factor(bin_id)
  
  bin_size <- as.numeric(table(bin_id))
  expr_agg <- rowsum(x = as.matrix(t(expr)), group = bin_id, reorder = FALSE)
  expr_agg <- expr_agg[levels(bin_id), ]
  expr_agg <- t(expr_agg / bin_size)
  
  coord_agg <- coord_agg[colnames(expr_agg), ]
  return(list(expr = expr_agg, coord = coord_agg)) 
}
