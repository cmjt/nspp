# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Calculating distances between points subject to periodic boundary constraints.
#'
#' Calculates pairwise distances between points subject to periodic boundary constraints.
#'
#' @inheritParams fit.ns
pbc_distances <- function(points, lims) {
    .Call('_nspp_pbc_distances', PACKAGE = 'nspp', points, lims)
}

