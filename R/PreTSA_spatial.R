#' Test if the gene expression is associated with spatial locations
#'
#' @param expr The normalized gene expression matrix. Rows represent genes and columns represent spots/cells
#' @param coord The matrix of spatial locations. Rows represent spots/cells. The column named "row" ("col") represents row (column) coordinates
#' @param knot Number of knots (0 by default) or \code{"auto"} for automatic selection
#' @param maxknotallowed A user-defined maximum number of knots (5 by default)
#'
#' @return A data frame with the p-value (FDR) and test statistic for each gene (each row)
#' @export
#'
spatialTest <- function(expr, coord, knot = 0, maxknotallowed = 5) {
  
  expr <- expr[, rownames(coord), drop = F]
  if(knot != "auto") {
    
    knotnum <- rep(knot, nrow(expr))
    names(knotnum) <- rownames(expr)
    xrow <- cbind(1, splines::bs(coord[,'row'], intercept = F, df = knot+3))
    ycol <- cbind(1, splines::bs(coord[,'col'], intercept = F, df = knot+3))
    
    B <- xrow[,rep(1:ncol(xrow), ncol(xrow))] * ycol[,rep(1:ncol(ycol), each = ncol(ycol))]
    B <- B[, which(matrixStats::colSds(B)>0), drop = F]
    B <- cbind(1, B)
    rownames(B) <- colnames(expr)
    colnames(B) <- NULL
    tBB <- crossprod(B)
    rownames(tBB) <- colnames(tBB) <- NULL
    
    pred <- as.matrix(expr %*% B %*% tcrossprod(chol2inv(chol(tBB)), B))
    
    SSE <- Matrix::rowSums((expr - pred)^2)
    SST <- Matrix::rowSums(sweep(expr, 1, Matrix::rowMeans(expr), FUN = '-')^2)
    
    if (any(SST<SSE)) print(names(which(SST<SSE)))
    
    fstat <- ((SST - SSE)/(ncol(B) - 1))/(SSE/(nrow(B) - ncol(B)))
    fstat[which(Matrix::rowSums(expr) == 0)] <- 0
    
    pval <- stats::pf(q = fstat, df1 = ncol(B) - 1, df2 = nrow(B) - ncol(B), lower.tail = F)
    logpval <- stats::pf(q = fstat, df1 = ncol(B) - 1, df2 = nrow(B) - ncol(B), lower.tail = F, log.p = T)
    
    res <- data.frame(fstat = fstat, pval = pval, logpval = logpval)
    
  } else {
    
    knotnum0 <- 0:maxknotallowed
    names(knotnum0) <- knotnum0
    
    Blist <- lapply(knotnum0, function(numknot) {
      
      xrow <- cbind(1, splines::bs(coord[,'row'], intercept = F, df = numknot+3))
      ycol <- cbind(1, splines::bs(coord[,'col'], intercept = F, df = numknot+3))
      
      B <- xrow[,rep(1:ncol(xrow), ncol(xrow))] * ycol[,rep(1:ncol(ycol), each = ncol(ycol))]
      B <- B[, which(matrixStats::colSds(B)>0), drop = F]
      B <- cbind(1, B)
      rownames(B) <- colnames(expr)
      colnames(B) <- NULL
      tBB <- crossprod(B)
      rownames(tBB) <- colnames(tBB) <- NULL
      list(B = B, tBB = tBB)
    })
    names(Blist) <- as.character(knotnum0)
    
    testpos <- sapply(knotnum0, function(numknot) {
      
      tBB <- Blist[[as.character(numknot)]][['tBB']]
      !'try-error'%in%class(try(chol(tBB), silent = T))
      #matrixcalc::is.positive.definite(tBB)
    })
    
    if(mean(testpos) != 1) {
      maxknot <- which(testpos == F)[1] - 2
      knotnum0 <- 0:maxknot
      names(knotnum0) <- knotnum0 }
    
    expr <- Matrix::t(expr)
    
    bic <- sapply(knotnum0, Calbic, Blist = Blist, expr = expr)
    
    knotnum <- knotnum0[apply(bic, 1, which.min)]
    names(knotnum) <- rownames(bic)
    
    res <- lapply(unique(knotnum), function(k) {
      
      B <- Blist[[as.character(k)]][['B']]
      tBB <- Blist[[as.character(k)]][['tBB']]
      
      beta <- as.matrix(tcrossprod(chol2inv(chol(tBB)), B) %*% expr[, which(knotnum == k), drop = F])
      pred <- B %*% beta
      
      expr.sub <- expr[, colnames(pred), drop = F]
      
      SSE <- Matrix::colSums((expr.sub - pred)^2)
      SST <- Matrix::colSums(sweep(expr.sub, 2, Matrix::colMeans(expr.sub), FUN = '-')^2)
      
      if (any(SST<SSE)) print(names(which(SST<SSE)))
      
      fstat <- ((SST - SSE)/(ncol(B) - 1))/(SSE/(nrow(B) - ncol(B)))
      fstat[which(Matrix::colSums(expr.sub) == 0)] <- 0
      
      pval <- stats::pf(q = fstat, df1 = ncol(B) - 1, df2 = nrow(B) - ncol(B), lower.tail = F)
      logpval <- stats::pf(q = fstat, df1 = ncol(B) - 1, df2 = nrow(B) - ncol(B), lower.tail = F, log.p = T)
      
      data.frame(fstat = fstat, pval = pval, logpval = logpval)
    })
    
    res <- do.call(rbind, res)
    res <- res[colnames(expr), ]
  }
  
  res$fdr <- stats::p.adjust(res$pval, method = 'fdr')
  res$knotnum <- knotnum
  res <- res[, c("fdr", "logpval", "pval", "fstat", "knotnum")]
  res <- res[order(res$pval, -res$fstat), ]
  return(res)
}

#' Fit the gene expression along spatial locations
#'
#' @param expr The normalized gene expression matrix. Rows represent genes and columns represent spots/cells
#' @param coord The matrix of spatial locations. Rows represent spots/cells. The column named "row" ("col") represents row (column) coordinates
#' @param knot Number of knots (0 by default) or \code{"auto"} for automatic selection
#' @param maxknotallowed A user-defined maximum number of knots (5 by default)
#'
#' @return The fitted expression matrix. Rows represent genes and columns represent spots/cells
#' @export
#'
spatialFit <- function(expr, coord, knot = 0, maxknotallowed = 5) {
  
  expr <- expr[, rownames(coord), drop = F]
  if(knot != "auto") {
    
    knotnum <- rep(knot, nrow(expr))
    names(knotnum) <- rownames(expr)
    xrow <- cbind(1, splines::bs(coord[,'row'], intercept = F, df = knot+3))
    ycol <- cbind(1, splines::bs(coord[,'col'], intercept = F, df = knot+3))
    
    B <- xrow[,rep(1:ncol(xrow), ncol(xrow))] * ycol[,rep(1:ncol(ycol), each = ncol(ycol))]
    B <- B[, which(matrixStats::colSds(B)>0), drop = F]
    B <- cbind(1, B)
    rownames(B) <- colnames(expr)
    colnames(B) <- NULL
    tBB <- crossprod(B)
    rownames(tBB) <- colnames(tBB) <- NULL
    
    pred <- as.matrix(expr %*% B %*% tcrossprod(chol2inv(chol(tBB)), B))
    return(pred)
  } else {
  
  knotnum0 <- 0:maxknotallowed
  names(knotnum0) <- knotnum0
  
  Blist <- lapply(knotnum0, function(numknot) {
    
    xrow <- cbind(1, splines::bs(coord[,'row'], intercept = F, df = numknot+3))
    ycol <- cbind(1, splines::bs(coord[,'col'], intercept = F, df = numknot+3))
    
    B <- xrow[,rep(1:ncol(xrow), ncol(xrow))] * ycol[,rep(1:ncol(ycol), each = ncol(ycol))]
    B <- B[, which(matrixStats::colSds(B)>0), drop = F]
    B <- cbind(1, B)
    rownames(B) <- colnames(expr)
    colnames(B) <- NULL
    tBB <- crossprod(B)
    rownames(tBB) <- colnames(tBB) <- NULL
    list(B = B, tBB = tBB)
  })
  names(Blist) <- as.character(knotnum0)
  
  testpos <- sapply(knotnum0, function(numknot) {
    
    tBB <- Blist[[as.character(numknot)]][['tBB']]
    !'try-error'%in%class(try(chol(tBB), silent = T))
    #matrixcalc::is.positive.definite(tBB)
  })
  
  if(mean(testpos) != 1) {
    maxknot <- which(testpos == F)[1] - 2
    knotnum0 <- 0:maxknot
    names(knotnum0) <- knotnum0 }
  
  expr <- Matrix::t(expr)
  
  bic <- sapply(knotnum0, Calbic, Blist = Blist, expr = expr)
  
  knotnum <- knotnum0[apply(bic, 1, which.min)]
  names(knotnum) <- rownames(bic)
  print(table(knotnum))
  
  pred <- lapply(unique(knotnum), function(k) {
    
    B <- Blist[[as.character(k)]][['B']]
    tBB <- Blist[[as.character(k)]][['tBB']]
    
    beta <- as.matrix(tcrossprod(chol2inv(chol(tBB)), B) %*% expr[, which(knotnum == k), drop = F])
    pred <- B %*% beta
    pred
  })
  
  pred <- do.call(cbind, pred)
  pred <- pred[, colnames(expr)]
  pred <- Matrix::t(pred)
  return(pred)
  }
}

Calbic <- function(numknot, Blist, expr) {
  
  B <- Blist[[as.character(numknot)]][['B']]
  tBB <- Blist[[as.character(numknot)]][['tBB']]
  
  beta <- as.matrix(tcrossprod(chol2inv(chol(tBB)), B) %*% expr)
  pred <- B %*% beta
  mse <- Matrix::colMeans((expr - pred)^2)
  bic <- nrow(B)*(1+log(2*pi)+log(mse)) + log(nrow(B))*(ncol(B)+1) # stats::AIC(lm(), k = log(nrow(B)))
  return(bic)
}
