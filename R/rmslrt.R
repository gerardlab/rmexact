#' Split likelihood ratio test for random mating
#'
#' Implements the split likelihood ratio test of Wasserman et al. (2020
#' when testing for random mating in autopolyploids.
#'
#' @param nvec The vector of genotype counts
#' @param sprop The proportional split.
#' @param log_p Should we return the log p-value?
#' @param nrep Number of iterations
#'
#' @return the p-value or the log of the p-value
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' p <- c(1/3, 1/3, 1/3)
#' q <- stats::convolve(p, rev(p), type = "open")
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' rmslrt(nvec, sprop = 0.5, nrep = 100)
#'
#' q <- c(0.25, 0.2, 0.1, 0.2, 0.25)
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' rmslrt(nvec, sprop = 0.5, nrep = 100)
#' }
#'
#' @references
#' \itemize{
#'   \item{Wasserman, L., Ramdas, A., & Balakrishnan, S. (2020). Universal inference. \emph{Proceedings of the National Academy of Sciences}, 117(29), 16880-16890. \doi{10.1073/pnas.1922664117}}
#' }
#'
#' @export
rmslrt <- function(nvec, sprop = 0.5, log_p = FALSE, nrep = 10) {
  ## Set params ----
  TOL <- sqrt(.Machine$double.eps)
  ploidy <- length(nvec) - 1
  nind <- sum(nvec)

  ## check args ----
  stopifnot(ploidy %% 2 < TOL)
  stopifnot(sprop > 1 / nind,
            sprop < (nind - 1) / nind,
            length(sprop) == 1)
  stopifnot(is.logical(log_p), length(log_p) == 1)

  ## Split ----
  lrtvec <- rep(NA_real_, length.out = nrep)
  for (i in seq_len(nrep)) {
    mapply(FUN = rep, x = 0:ploidy, times = nvec) |>
      unlist() |>
      sample(size = round(nind * sprop), replace = FALSE) |>
      factor(levels = 0:ploidy) |>
      table() |>
      as.vector() ->
      nvec1
    nvec2 <- nvec - nvec1

    ## LRT ----
    q1 <- nvec1 / sum(nvec1)
    p2 <- hwep::rmlike(nvec = nvec2, thresh = 0)$p
    q2 <- stats::convolve(p2, rev(p2), type = "open")
    q2[q2 < 0] <- 0
    q2 <- q2 / sum(q2)

    lrtvec[[i]] <- stats::dmultinom(x = nvec2, prob = q1, log = TRUE) -
      stats::dmultinom(x = nvec2, prob = q2, log = TRUE)
  }
  lpval <- -log_sum_exp(lrtvec) + log(nrep)
  lpval <- min(lpval, 0)

  if(!log_p) {
    lpval <- exp(lpval)
  }
  return(list(q1 = q1, q2 = q2, nvec1 = nvec1, nvec2 = nvec2, pval = lpval))
}
