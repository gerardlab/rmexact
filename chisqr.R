#'Calculate Chi-squared for tetraploid
#'
#' @param nvec vector of length 3
#' @return the p-value
#' @example
#'
#'ifelse(log_p = TRUE, log_p = log_p, p <- log_p)
#'#we will have to run rmultinom first
#'result <- rmultinom(n = 1, size = 10, prob = c(0.3,0.1,0.6))
library(rmexact)
#if log_p is true, we would like to get the log of the p_value
rmchisq <- function(pvec, log_p = FALSE, N) {
  freq <- gfreq(pvec = pvec)
  result <- rmultinom(n = 1, size = N, prob = freq)
  nvec <- as.vector(nvec)
  mle <- hwep::rmlike(nvec = nvec)$p
  freq1 <- gfreq(pvec = mle)
  ecounts <- N * freq1 #Expected
  p <- sum(((nvec[1] - ecounts[1])^2 / ecounts[1]),
            ((nvec[2] - ecounts[2])^2 / ecounts[2]),
             ((nvec[3] - ecounts[3])^2 / ecounts[3]),
             ((nvec[4] - ecounts[4])^2 / ecounts[4]),
             ((nvec[5] - ecounts[5])^2 / ecounts[5]))
  result <- pchisq(q = p, df = 1)
  if (log_p == TRUE) {
    freq <- gfreq(pvec = pvec)
    result <- rmultinom(n = 1, size = N, prob = freq)
    nvec <- as.vector(nvec)
    mle <- hwep::rmlike(nvec = nvec)$p
    freq1 <- gfreq(pvec = mle)
    ecounts <- N * freq1 #Expected
    p <- sum(((nvec[1] - ecounts[1])^2 / ecounts[1]),
               ((nvec[2] - ecounts[2])^2 / ecounts[2]),
               ((nvec[3] - ecounts[3])^2 / ecounts[3]),
               ((nvec[4] - ecounts[4])^2 / ecounts[4]),
               ((nvec[5] - ecounts[5])^2 / ecounts[5]))
    p_value <- pchisq(q = p, df = 1)
    result <- log(p_value)
  }
  return(result)
}
rmchisq(pvec = loc, log_p = FALSE, N = 10)
