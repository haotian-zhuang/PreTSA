#' Test if the gene expression is associated with pseudotime values
#'
#' @param expr The normalized gene expression matrix. Rows represent genes and columns represent cells
#' @param pseudotime The vector of user-provided pseudotime values
#' @param pseudotime_permute A list of permuted pseudotime values from subsampled cells. Each element in the list has the same format of the argument `pseudotime`
#' @param knot Whether to select the optimal number of knots automatically
#' @param maxknotallowed A user-defined maximum number of knots (10 by default)
#'
#' @return A data frame with the p-value and FDR for each gene (each row)
#' @export
#'
temporalTest<- function(expr, pseudotime, pseudotime_permute = NULL, knot = F, maxknotallowed = 10) {
  
  if(is.null(pseudotime_permute)) {
    return(temporalTest.fix(expr = expr, pseudotime = pseudotime,
                            knot = knot, maxknotallowed = maxknotallowed))
  } else {
    fstat.ori <- Calfstat(expr = expr, pseudotime = pseudotime,
                          knot = knot, maxknotallowed = maxknotallowed)
    
    fstat.perm <- sapply(pseudotime_permute, function(i) {
      Calfstat(expr = expr, pseudotime = i,
               knot = knot, maxknotallowed = maxknotallowed)
    })
    
    pval.empirical <- (Matrix::rowSums(fstat.perm >= fstat.ori) + 1)/(ncol(fstat.perm) + 1)
    
    pval.parametric <- sapply(1:nrow(fstat.perm), function(i) {
      if(!'try-error'%in%class(try(suppressWarnings(fitdistrplus::fitdist(fstat.perm[i, ], 'gamma')), silent = T))) {
        fit.gamma <- suppressWarnings(fitdistrplus::fitdist(fstat.perm[i, ], 'gamma'))
        return(stats::pgamma(fstat.ori[i], shape = fit.gamma$estimate[1], rate = fit.gamma$estimate[2], lower.tail = F))
      } else {
        return(NA)
      }
    })
    
    fdr.empirical <- stats::p.adjust(pval.empirical, method = 'fdr')
    fdr.parametric <- stats::p.adjust(pval.parametric, method = 'fdr')
    res <- data.frame(fdr.parametric = fdr.parametric, pval.parametric = pval.parametric,
                      fdr.empirical = fdr.empirical, pval.empirical = pval.empirical)
    return(res)
  }
}

Calfstat <- function(expr, pseudotime, knot = F, maxknotallowed = 10){
  
  expr <- expr[, names(pseudotime), drop = F]
  if(knot == F) {
    
    knotnum <- rep(0, nrow(expr))
    names(knotnum) <- rownames(expr)
    B <- splines::bs(pseudotime, intercept = F, df = 3)
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
    return(fstat)
    
  } else {
    
    knotnum0 <- 0:maxknotallowed
    names(knotnum0) <- knotnum0
    
    Blist <- lapply(knotnum0, function(numknot) {
      
      B <- splines::bs(pseudotime, intercept = F, df = numknot+3)
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
    
    fstat <- lapply(unique(knotnum), function(k) {
      
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
      fstat
    })
    
    fstat <- unlist(fstat)
    fstat <- fstat[colnames(expr)]
    return(fstat)
  }
}

temporalTest.fix <- function(expr, pseudotime, knot = F, maxknotallowed = 10) {
  
  expr <- expr[, names(pseudotime), drop = F]
  if(knot == F) {
    
    knotnum <- rep(0, nrow(expr))
    names(knotnum) <- rownames(expr)
    B <- splines::bs(pseudotime, intercept = F, df = 3)
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
    pval <- stats::pf(q = fstat, df1 = ncol(B) - 1, df2 = nrow(B) - ncol(B), lower.tail = F)
    
    pval[which(Matrix::rowSums(expr) == 0)] <- 1
    
  } else {
    
    knotnum0 <- 0:maxknotallowed
    names(knotnum0) <- knotnum0
    
    Blist <- lapply(knotnum0, function(numknot) {
      
      B <- splines::bs(pseudotime, intercept = F, df = numknot+3)
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
    
    pval <- lapply(unique(knotnum), function(k) {
      
      B <- Blist[[as.character(k)]][['B']]
      tBB <- Blist[[as.character(k)]][['tBB']]
      
      beta <- as.matrix(tcrossprod(chol2inv(chol(tBB)), B) %*% expr[, which(knotnum == k), drop = F])
      pred <- B %*% beta
      
      expr.sub <- expr[, colnames(pred), drop = F]
      
      SSE <- Matrix::colSums((expr.sub - pred)^2)
      SST <- Matrix::colSums(sweep(expr.sub, 2, Matrix::colMeans(expr.sub), FUN = '-')^2)
      
      if (any(SST<SSE)) print(names(which(SST<SSE)))
      
      fstat <- ((SST - SSE)/(ncol(B) - 1))/(SSE/(nrow(B) - ncol(B)))
      pval <- stats::pf(q = fstat, df1 = ncol(B) - 1, df2 = nrow(B) - ncol(B), lower.tail = F)
      
      pval[which(Matrix::colSums(expr.sub) == 0)] <- 1
      pval
    })
    
    pval <- unlist(pval)
    pval <- pval[colnames(expr)]
  }
  
  fdr <- stats::p.adjust(pval, method = 'fdr')
  res <- data.frame(fdr = fdr, pval = pval, knotnum = knotnum)
  return(res)
}

#' Fit the gene expression along pseudotime values
#'
#' @param expr The normalized gene expression matrix. Rows represent genes and columns represent cells
#' @param pseudotime The vector of user-provided pseudotime values
#' @param knot Whether to select the optimal number of knots automatically
#' @param maxknotallowed A user-defined maximum number of knots (10 by default)
#'
#' @return The fitted expression matrix. Rows represent genes and columns represent cells
#' @export
#'
temporalFit <- function(expr, pseudotime, knot = F, maxknotallowed = 10) {
  
  expr <- expr[, names(pseudotime), drop = F]
  if(knot == F) {
    
    knotnum <- rep(0, nrow(expr))
    names(knotnum) <- rownames(expr)
    B <- splines::bs(pseudotime, intercept = F, df = 3)
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
    
    B <- splines::bs(pseudotime, intercept = F, df = numknot+3)
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

Calbic<- function(numknot, Blist, expr) {
  
  B <- Blist[[as.character(numknot)]][['B']]
  tBB <- Blist[[as.character(numknot)]][['tBB']]
  
  beta <- as.matrix(tcrossprod(chol2inv(chol(tBB)), B) %*% expr)
  pred <- B %*% beta
  mse <- Matrix::colMeans((expr - pred)^2)
  bic <- nrow(B)*(1+log(2*pi)+log(mse)) + log(nrow(B))*(ncol(B)+1) # stats::AIC(lm(), k = log(nrow(B)))
  return(bic)
}
