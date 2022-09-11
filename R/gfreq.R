#'Calculate genotype frequency from gametes frequencies
#'
#' @param pvec A vector of length three (gamete frequencies) that sum to 1
#'
#' @return Genotype frequency.
#' @example 
#' pvec <- c(0.1,0.4,0.5)
#' gfreq(pvec)
#' gfreq <- function(pvec) {
qvec <- stats::convolve(pvec, rev(pvec), type = "open")
qvec[qvec < 0] <- 0
qvec <- qvec / sum(qvec)
return(qvec)
}