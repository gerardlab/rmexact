#' Universal conference
#' @examples 
#' set.seed(1)
#' pvec <- c(1/3, 1/3, 1/3)
#' n <- 100
#' freq <- rmexact::gfreq(pvec = pvec)
#' nvec <- c(rmultinom(n = 1, size = n, prob = freq))
#' nvec <- c(1000, 100, 0, 1000, 0)
#' universal(nvec = nvec)
##'@export
universal <- function(nvec, log_pvalue = FALSE) {
  ploidy <- length(nvec) - 1
  mapply(0:ploidy, nvec, FUN = rep) |>
    unlist() |>
    sample(size = round(sum(nvec) / 2), replace = FALSE) |>
    factor(levels = 0:ploidy) |>
    table() |>
    as.vector() ->
    group1
  group2 <- nvec - group1
  
  mle1 <- hwep::rmlike(nvec = group1)$p
  qvec1 <- rmexact::gfreq(pvec = mle1)
  mle2 <- hwep::rmlike(nvec = group2)$p
  qvec2 <- rmexact::gfreq(pvec = mle2)
  lrs1 <- dmultinom(x = group1, size = sum(group1), prob = qvec1, log = TRUE) - 
    dmultinom(x = group1, size = sum(group1), prob = qvec2, log = TRUE)
    
  lrs2 <- dmultinom(x = group2, size = sum(group2), prob = qvec2, log = TRUE) -
    dmultinom(x = group2, size = sum(group2), prob = qvec1, log = TRUE)
  
  lpval <- min(log(2) - log_sum_exp(c(lrs1, lrs2)), 0)
  
  if (!log_pvalue) {
    lpval <- exp(lpval)
  }
  return(lpval)
}


