#' Move genotype counts up one.
#'
#' @param y A vector of length 5. \code{y[k]} is the number of individuals with
#'     genotype \code{k}.
#'
#' @author David Gerard
#'
#' @noRd
tetup <- function(y) {
  stopifnot(length(y) == 5)
  y + c(1, -2, 0, 2, -1)
}

#' Move genotype counts down one
#'
#' @param y A vector of length 5. \code{y[k]} is the number of individuals with
#'     genotype \code{k}.
#'
#' @author David Gerard
#'
#' @noRd
tetdown <- function(y) {
  stopifnot(length(y) == 5)
  y + c(-1, 2, 0, -2, 1)
}

#' Proportional Conditional probability
#'
#' Proportional to actual conditional probability (given sufficient statistics
#' for random mating in tetraploids).
#'
#' @param y A vector of length 5. \code{y[k]} is the number of individuals with
#'     genotype \code{k}.
#' @param log_p A logical. Should we return the log (\code{TRUE}) or not
#'     (\code{FALSE})?
#'
#'  @author David Gerard
#'
#' @noRd
condprob <- function(y, log_p = TRUE) {
  stopifnot(length(y) == 5)
  stopifnot(is.logical(log_p), length(log_p) == 1)
  cp <- lgamma(sum(y) + 1) + y[[2]] * log(2) + y[[4]] * log(2) - sum(lgamma(y + 1))
  if (!log_p) {
    cp <- exp(cp)
  }
  return(cp)
}

#' Exact test for tetraploid random mating
#'
#' Conditions on sufficient statistics to calculate an exact p-value
#' against the null of random mating.
#'
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
#' ## small change when using fractional counts
#' y <- y + runif(5)
#' tetexact(y = y)
#'
#' @export
tetexact <- function(y, log_p = FALSE, frac = FALSE) {
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

  ## current y
  cy <- c(y[[1]] + umax, y[[2]] - 2 * umax, y[[3]], y[[4]] + 2 * umax, y[[5]] - umax)
  hvec <- rep(NA_real_, length.out = m)
  hvec[[1]] <- condprob(cy)
  for (i in 2:m) {
    cy <- tetdown(cy)
    hvec[[i]] <- condprob(cy)
  }
  hsum <- log_sum_exp(hvec)
  hvec <- hvec - hsum

  hobs <- condprob(y)
  hobs <- hobs - hsum

  pval <- log_sum_exp(hvec[hvec <= hobs])

  if (!log_p) {
    pval <- exp(pval)
  }

  return(pval)
}


#' Fuzzy p-values for random mating using genotype likelihoods.
#'
#' Returns samples from the distribution of p-values. This distribution
#' is collectively called the "fuzzy p-value" (Geyer and Meeden, 2005).
#'
#' @param gp The matrix of genotype posteriors. \code{gp[i,k]} is the
#'     posterior probability that individual \code{i} has genotype
#'     \code{k - 1}.
#' @param nrep The number of type to run \code{\link{tetexact}()}.
#' @param log_p A logical. Should we return the log p-value (\code{TRUE}) or not
#'     (\code{FALSE})?
#'
#' @return A sample from the distribution of the abstract randomized p-value.
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{C. J. Geyer and G. D. Meeden. Fuzzy and Randomized Confidence Intervals and P-Values. \emph{Statistical Science}, 20(4):358 – 366, 2005. \doi{10.1214/088342305000000340}.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(99)
#'
#' # Generate data according to random mating
#' p <- c(0.5, 0.3, 0.2)
#' q <- stats::convolve(p, rev(p), type = "open")
#' y <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' gp <- hwep::simgl(nvec = y, rdepth = 1000, ret = "gp")
#' pval1 <- tetgp(gp)
#' hist(pval1, breaks = 20, xlim = c(0, 1)) ## no evidence against null
#' vpv(pval1) ## vacuumed p-value
#'
#' # Generate data different from random mating
#' q <- c(0.4, 0.3, 0.2, 0.1, 0.0)
#' y <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' gp <- hwep::simgl(nvec = y, rdepth = 1000, ret = "gp")
#' pval2 <- tetgp(gp)
#' hist(pval2, breaks = 20, xlim = c(0, 1)) ## evidence against null
#' vpv(pval2) ## vacuumed p-value
#' }
#'
#' @export
tetgp <- function(gp, nrep = 1000, log_p = FALSE) {
  TOL <- sqrt(.Machine$double.eps)
  stopifnot(ncol(gp) == 5)
  stopifnot(abs(rowSums(gp) - 1) < TOL)
  pvec <- rep(NA_real_, length.out = nrep)

  for (i in seq_len(nrep)) {
    y <- table(factor(apply(gp, 1, sample, x = 0:4, size = 1, replace = FALSE), 0:4))
    pvec[[i]] <- tetexact(y = y, log_p = log_p, frac = FALSE)
  }

  return(pvec)
}


#' Vacuum the fuzzy p-values
#'
#' @param pvals A vector of sampled abstract randomized p-values
#'      (Geyer and Meeden, 2005).
#'
#' @return The vacuumed p-value (VP-value)
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{C. J. Geyer and G. D. Meeden. Fuzzy and Randomized Confidence Intervals and P-Values. \emph{Statistical Science}, 20(4):358 – 366, 2005. \doi{10.1214/088342305000000340}.}
#' }
#'
#' @examples
#' set.seed(1)
#' pvals <- runif(1000)
#' vpv(pvals)
#'
#' @export
vpv <- function(pvals) {
  pvals <- sort(pvals)
  n <- length(pvals)
  which_q <- match(TRUE, pvals >= 1 - seq_len(n) / n)
  if (is.na(which_q)) {
    return(1)
  } else {
    return(pvals[[which_q]])
  }
}
