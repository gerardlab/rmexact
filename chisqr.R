#'Calculate Chi-squared for tetraploids
#'
#' @param nvec vector of length 5
#' @param log_pvalue if the user desires to get the log of the p value
#' @return the p-value or the log of the p-value depending on log_pvalue
#' @example
#'nvec <- c(0,6,30,37,27)
#'rmchisq(nvec, log_pvalue = TRUE)
rmchisq <- function(nvec, log_pvalue = FALSE) {
  stopifnot(nvec == as.vector(nvec))
  counts <- nvec
  pvec <- hwep::rmlike(nvec = nvec)$p
  qvec <- gfreq(pvec = pvec)
  ecounts <- sum(nvec) * qvec
  q <- sum((counts - ecounts)^2/ecounts)
  result <- pchisq(q = q, df = 2, lower.tail = FALSE)
  if (log_pvalue == TRUE) {
    result <- log(result)
  }
  return(result)
}
