#' Exact tests for random mating in autopolyploids.
#'
#' Right now, we only have an implementation for tetraploids. But we can
#' account for genotype uncertainty through genotype posteriors.
#'
#' The main functions are
#' \describe{
#'   \item{\code{\link{tetexact}()}}{Calculate the exact p-value when genotypes are known.}
#'   \item{\code{\link{tetgp}()}}{Calculate the fuzzy p-value when using genotype posteriors.}
#'   \item{\code{\link{vpv}()}}{Calculate the vacuumed p-value.}
#' }
#'
#' @keywords internal
#' @aliases rmexact
#'
#' @author David Gerard and Karene Matoka Nana
#'
#' @importFrom hwep simgl
#'
"_PACKAGE"
