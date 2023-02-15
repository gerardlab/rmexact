#' Calculate genotype frequencies from gamete frequencies
#'
#' These genotype frequencies are calculated  under the assumption of
#' random mating where each frequency is generated through a convolution
#' of the gamete frequencies.
#'
#' @param pvec A vector of length three (gamete frequencies) that sum to 1
#'
#' @return Genotype frequency.
#' @examples
#' pvec <- c(0.1,0.4,0.5)
#' gfreq(pvec)
#' @export
gfreq <- function(pvec) {
  qvec <- stats::convolve(pvec, rev(pvec), type = "open")
  qvec[qvec < 0] <- 0
  qvec <- qvec / sum(qvec)
  return(qvec)
}

