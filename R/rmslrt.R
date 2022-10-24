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
  return(lpval)
}

#Individual Example
q <- c(0.25, 0.2, 0.1, 0.2, 0.25)
nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
sprop <- 0.5
#First step
TOL <- sqrt(.Machine$double.eps)
ploidy <- length(nvec) - 1
nind <- sum(nvec)


##Second step
mapply(FUN = rep, x = 0:ploidy, times = nvec) |>
  unlist() |>
  sample(size = round(nind * sprop), replace = FALSE) |>
  factor(levels = 0:ploidy) |>
  table() |>
  as.vector() ->
  nvec1
nvec2 <- nvec - nvec1

##Third Step
q1 <- nvec1 / sum(nvec1)
p2 <- hwep::rmlike(nvec = nvec2, thresh = 0)$p
q2 <- stats::convolve(p2, rev(p2), type = "open")
q2[q2 < 0] <- 0
q2 <- q2 / sum(q2)

lpval <- -stats::dmultinom(x = nvec2, prob = q1, log = TRUE) +
  stats::dmultinom(x = nvec2, prob = q2, log = TRUE)
lpval <- min(lpval, 0)

#' nvec2 values are 11 10 4 10 15
#' q1 0.24 0.24 0.10 0.26 0.16
#' q2 0.1260489 0.1490833 0.3529680 0.1826663 0.1892335
#' The test statistics have to be q1 and q2
#' Maybe run the test in a loop, extracting the test statistics each iteration, 
#' average the log of the test statistics, then calculate the p-value:
plisttrial <- list(c(0.333, 0.333, 0.333),
              c(0, 1, 0),
              c(1/2, 0, 1/2),
              c(0.1, 0.3, 0.6))

nvectrial <- c(10, 100)
seedvec <- seq_len(10)
trialrmslrt <- expand.grid(p = plisttrial,
                       n = nvectrial,
                       seed = seedvec)
for(i in seq_len(nrow(trialrmslrt))) {
  set.seed(paramdf$seed[[i]])
  p <- paramdf$p[[i]]
  n <- paramdf$n[[i]]
  freq <- rmexact::gfreq(pvec = p)
  result <- c(rmultinom(n = 1, size = n, prob = freq))
  teststatistics1 <- rmslrt(nvec, sprop = 0.5)$q1 #Extract the first test statistics
  teststatistics2 <- rmslrt(nvec, sprop = 0.5)$q2 #Extract the second test statistics
  avg <- log(q1 + q2) - log(2) #Average the log of the test statistics
  lpval <- -stats::dmultinom(x = nvec2, prob = q1, log = TRUE) +
    stats::dmultinom(x = nvec2, prob = q2, log = TRUE)
}

#Simple Test
first <- c(2,5,7)
second <- c(1,5,8)
res1 <- log(first)
res2 <- log(second)
log(first + second) - log(2)

#'Due to the failures of running universal, we may need to first run `rmslrt()`,
#'inside a for loop, then extract the test statistics. I am thinking that I could 
#'directly use the calculations used inside `universal()` to calculate the pvals.
#'Or I would have to do it manually and then used that calculations. 
