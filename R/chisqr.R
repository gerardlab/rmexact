#'Calculate Chi-squared for tetraploid
#'
#' @param pvec vector of length 3
#' @param log_p if the user desires to get the log of the p value
#' @param N the sample size
#' @return the p-value
#' @example
#'rmchisq(pvec = loc, log_p = FALSE, N = 10)
#'rmchisq(pvec = loc, log_p = TRUE, N = 10)

rmchisq <- function(pvec, log_p = FALSE, N) {
  freq <- gfreq(pvec = pvec)
  result <- rmultinom(n = 1, size = N, prob = freq)
  nvec <- as.vector(nvec) #nvec = genotype counts
  mle <- hwep::rmlike(nvec = nvec)$p  #Extract mle
  freq1 <- gfreq(pvec = mle)
  ecounts <- N * freq1 #Expected counts
  p <- sum(((nvec[1] - ecounts[1])^2 / ecounts[1]),
            ((nvec[2] - ecounts[2])^2 / ecounts[2]),
             ((nvec[3] - ecounts[3])^2 / ecounts[3]),
             ((nvec[4] - ecounts[4])^2 / ecounts[4]),
             ((nvec[5] - ecounts[5])^2 / ecounts[5]))
  result <- pchisq(q = p, df = 1)
  if (log_p == TRUE) {
    freq <- gfreq(pvec = pvec)
    result <- rmultinom(n = 1, size = N, prob = freq)
    nvec <- as.vector(nvec) #nvec = genotype counts
    mle <- hwep::rmlike(nvec = nvec)$p
    freq1 <- gfreq(pvec = mle)
    ecounts <- N * freq1 #Expected counts
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
