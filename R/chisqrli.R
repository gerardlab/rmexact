#' Chi-Squared based on Li Paper
#'
#' Using Conditions on sufficient statistics to calculate the chi-squared
#' p-value.
#'
#
#' If you have genotype likelihoods, I would first convert those to genotype
#' posteriors (e.g. using \code{ldsep::gl_to_gp()}) to obtain \code{gp}, then
#' run \code{y <- colSums(gp, na.rm = TRUE)}. Then use the \code{frac = TRUE}?
#'
#' @param y The vector of genotype counts (of of length 5). \code{y[k]} is
#'     the number of individuals with genotype \code{k-1}.
#' @param log_p A logical. Should we return the log p-value (\code{TRUE}) or not
#'     (\code{FALSE})?
#' @param frac A logical. Does \code{y} contain fractional genotypes (e.g.
#'     because \code{y} is calculated by summing up posterior probabilities)?
#'
#' @author David Gerard and Karene Matoka Nana
#'
#' @examples
#' set.seed(1)
#'
#' # Generate data according to random mating
#' p <- c(0.5, 0.3, 0.2)
#' q <- stats::convolve(p, rev(p), type = "open")
#' y <- c(stats::rmultinom(n = 1, size = 200, prob = q))
#' tetexact(y = y)
#'
#' # Generate data different from random mating
#' q <- c(0.4, 0.3, 0.2, 0.1, 0.0)
#' y <- c(stats::rmultinom(n = 1, size = 200, prob = q))
#' tetexact(y = y)
#'
#' @export
chisqrli <- function(y, log_p = FALSE, frac = FALSE) {
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

  ## Get conditional distribution of counts
  cy <- c(y[[1]] + umax, y[[2]] - 2 * umax, y[[3]], y[[4]] + 2 * umax, y[[5]] - umax)
  hvec <- rep(NA_real_, length.out = m)
  hvec[[1]] <- condprob(cy)
  for (i in 2:m) {
    cy <- tetdown(cy)
    hvec[[i]] <- condprob(cy)
  }
  hsum <- log_sum_exp(hvec)
  hvec <- hvec - hsum

  ## Go through loop again to calculate expectation
  cy <- c(y[[1]] + umax, y[[2]] - 2 * umax, y[[3]], y[[4]] + 2 * umax, y[[5]] - umax)
  ek <- cy * exp(hvec[[1]])
  for (i in 2:m) {
    cy <- tetdown(cy)
    ek <- ek + cy * exp(hvec[[i]])
  }

  ## Chi-squared test and p-value
  X2 <- sum((y - ek)^2/ek)
  pval <- stats::pchisq(q = X2, df = 1, lower.tail = FALSE, log.p = TRUE)
  if (!log_p) {
    pval <- exp(pval)
  }
  return(pval)
}
