#'Calculate Chi-squared for tetraploid
#'
#' @param nvec vector of length 3
#' @return the p-value
#' @example 
#'
rmchisq <- function(nvec, log_p = TRUE) {
  if (log_p = TRUE) {
  n <- sum(nvec)
  mle <- hwep::rmlike(nvec = nvec)$p
  q <- gfreq(pvec = pvec)
  ecounts <- n * q #Expected
  p <- sum(((nvec[1] - ecounts[1])^2 / ecounts[1]) + ((nvec[2] - ecounts[2])^2 / ecounts[2]) + ((nvec[3] - ecounts[3])^2 / ecounts[3]))
  p_value <- pchisq(q = p, df = 1)
  result <- log(p_value)
  } else {
  n <- sum(nvec)
  mle <- hwep::rmlike(nvec = nvec)$p
  q <- gfreq(pvec = pvec)
  ecounts <- n * q #Expected
  p <- sum(((nvec[1] - ecounts[1])^2 / ecounts[1]) + ((nvec[2] - ecounts[2])^2 / ecounts[2]) + ((nvec[3] - ecounts[3])^2 / ecounts[3]))
  result <- pchisq(q = p, df = 1) 
  }
  return(result)
}