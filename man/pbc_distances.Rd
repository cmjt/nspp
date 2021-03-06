% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pbc_distances}
\alias{pbc_distances}
\title{Calculating distances between points subject to periodic boundary constraints.}
\usage{
pbc_distances(points, lims)
}
\arguments{
\item{points}{A matrix containing locations of observed points,
where each row corresponds to a point and each column corresponds
to a dimension.}

\item{lims}{A matrix with two columns, corresponding to the upper
and lower limits of each dimension, respectively.}
}
\description{
Calculates pairwise distances between points subject to periodic boundary constraints.
}

