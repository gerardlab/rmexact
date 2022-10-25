#' Split likelihood ratio test for random mating
#'
#' @param nvec The vector of genotype counts
#' @param sprop The proportional split.
#' @param log_p Should we return the log p-value?
#'
#' @author David Gerard and Karene Matoka Nana
#'
#' @examples
#' set.seed(1)
#' p <- c(1/3, 1/3, 1/3)
#' q <- stats::convolve(p, rev(p), type = "open")
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' rmslrt(nvec, sprop = 0.5)
#'
#' q <- c(0.25, 0.2, 0.1, 0.2, 0.25)
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' rmslrt(nvec, sprop = 0.5)
#'
#' @export
rmslrt <- function(nvec, sprop = 0.5, log_p = FALSE) {
  ## Set params ----
  TOL <- sqrt(.Machine$double.eps)
  ploidy <- length(nvec) - 1
  nind <- sum(nvec)
  
  ## check args ----
  stopifnot(-./ploidy %% 2 < TOL)
  stopifnot(sprop > 1 / nind,
            sprop < (nind - 1) / nind,
            length(sprop) == 1)
  stopifnot(is.logical(log_p), length(log_p) == 1)
  
  ## Split ----
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
  
  lpval <- -stats::dmultinom(x = nvec2, prob = q1, log = TRUE) +
    stats::dmultinom(x = nvec2, prob = q2, log = TRUE)
  lpval <- min(lpval, 0)
  
  if(!log_p) {
    lpval <- exp(lpval)
  }
  return(c(q1, q2, nvec1, nvec2, lpval))
}
