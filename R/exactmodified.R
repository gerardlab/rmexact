tetexactmodified <- function(y, log_p = FALSE, frac = FALSE, midpval = FALSE) {
  TOL <- sqrt(.Machine$double.eps)
  stopifnot(length(y) == 5)
  stopifnot(is.logical(log_p), length(log_p) == 1)
  stopifnot(is.logical(frac), length(frac) == 1)
  
  n <- sum(y)
  
  if (!frac) {
    umax <- min(n - y[[3]] - y[[1]], floor(y[[2]] / 2), floor((n - y[[3]] - y[[4]])/2), y[[5]])
    umin <- max(-y[[1]], ceiling((y[[3]] + y[[2]] - n)/2), ceiling(-y[[4]]/2), y[[3]] + y[[5]] - n)
    m <- umax - umin + 1
  } else {
    umax <- min(n - y[[3]] - y[[1]], y[[2]] / 2, (n - y[[3]] - y[[4]])/2, y[[5]])
    umin <- max(-y[[1]], (y[[3]] + y[[2]] - n)/2, -y[[4]]/2, y[[3]] + y[[5]] - n)
    m <- floor(umax - umin) + 1
  }
  
  stopifnot(m > 0)
  if (abs(m - 1) < TOL) {
    if (!log_p) {
      return(1)
    } else {
      return(0)
    }
  }
  
  ## current y
  cy <- c(y[[1]] + umax, y[[2]] - 2 * umax, y[[3]], y[[4]] + 2 * umax, y[[5]] - umax)
  hvec <- rep(NA_real_, length.out = m)
  hvec[[1]] <- condprob(cy)
  for (i in 2:m) {
    cy <- tetdown(cy)
    hvec[[i]] <- condprob(cy)
  }
  hsum <- log_sum_exp(hvec)
  hvec <- hvec - hsum
  
  hobs <- condprob(y)
  hobs <- hobs - hsum
  
  pval <- log_sum_exp(hvec[hvec <= hobs])
  
  if (!log_p) {
    pval <- exp(pval)
  }
  if (!midpval) {
    pval <- exp((1/2) * hobs) + pval 
  }
  
  return(pval)
}