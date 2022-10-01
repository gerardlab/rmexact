#'Calculate Chi-squared for tetraploid
#'
#' @param nvec vector of length 5
#' @param log_pvalue if the user desires to get the log of the p value
#' @return the p-value or the log of the p-value depending on log_pvalue
#' @examples
#'nvec <- c(2, 1, 2, 3, 0)
#'rmchisq(nvec, log_pvalue = FALSE)
#'hwep::rmlike(nvec = nvec, thresh = 0)$p_rm
#'@export
rmchisq <- function(nvec, log_pvalue = FALSE) {
  stopifnot(nvec == as.vector(nvec))
  pvec <- hwep::rmlike(nvec = nvec)$p
  qvec <- gfreq(pvec = pvec)
  ecounts <- sum(nvec) * qvec
  q <- sum((nvec - ecounts)^2/ecounts)
  result <- pchisq(q = q, df = 2, lower.tail = FALSE, log.p = TRUE)
  if (!log_pvalue) {
    result <- exp(result)
  }
  return(result)
}
