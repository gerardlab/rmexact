#' Universal conference
#' @example
#' set.seed(1)
#' pvec <- c(1/3, 1/3, 1/3)
#' n <- 100
#' freq <- rmexact::gfreq(pvec = pvec)
#' group1 <- c(rmultinom(n = 1, size = n/2, prob = freq))
#' group2 <- c(rmultinom(n = 1, size = n/2, prob = freq))
#' universal(group1 = as.vector(group1), group2 = as.vector(group2))
##'@export
universal <- function(group1, group2, log_pvalue = FALSE) {
  stopifnot(group1 == as.vector(group1))
  stopifnot(group2 == as.vector(group2))
  mle1 <- hwep::rmlike(nvec = group1)$p
  qvec1 <- rmexact::gfreq(pvec = mle1)
  mle2 <- hwep::rmlike(nvec = group2)$p
  qvec2 <- rmexact::gfreq(pvec = mle2)
  Tn <- prod(group1 * qvec2)/prod(group1 * qvec1)
  Tnswap <- prod(group2 * qvec2)/prod(group2 * qvec1)
  sn <- (Tn + Tnswap)/2
  return(sn)
}


universal(group1 = as.vector(group1), group2 = as.vector(group2))
