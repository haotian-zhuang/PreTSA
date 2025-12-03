#' Optional temporal binning
#'
#' @param expr The normalized gene expression matrix. Rows represent genes and columns represent cells
#' @param pseudotime The vector of user-provided pseudotime values
#' @param num_bin Number of bins (100 by default)
#'
#' @returns A list of binned `expr` and `pseudotime`
#' @export
#'
bin_pseudotime <- function(expr, pseudotime, num_bin = 100) {
  
  tbreak <- seq(from = min(pseudotime), to = max(pseudotime), length.out = num_bin + 1)
  
  tbin <- cut(x = pseudotime, breaks = tbreak, labels = FALSE, include.lowest = TRUE, right = TRUE)
  
  pseudotime_agg <- (tbreak[-1] + tbreak[-length(tbreak)]) / 2
  names(pseudotime_agg) <- paste0("Bin:", seq_len(length(pseudotime_agg)))
  
  bin_id <- paste0("Bin:", tbin)
  names(bin_id) <- names(pseudotime)
  bin_id <- factor(bin_id)
  
  bin_size <- as.numeric(table(bin_id))
  expr_agg <- rowsum(x = as.matrix(t(expr)), group = bin_id, reorder = FALSE)
  expr_agg <- expr_agg[levels(bin_id), ]
  expr_agg <- t(expr_agg / bin_size)
  
  pseudotime_agg <- pseudotime_agg[colnames(expr_agg)]
  return(list(expr = expr_agg, pseudotime = pseudotime_agg))
}
