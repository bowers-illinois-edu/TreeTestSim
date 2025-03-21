% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_table_approach_alpha_varies.R
\name{local_min_p}
\alias{local_min_p}
\title{Unadjusted local step}
\usage{
local_min_p(pvals_children, alpha = 0.05)
}
\arguments{
\item{pvals_children}{Numeric vector of child p-values.}

\item{alpha}{Numeric scalar of alpha}
}
\value{
A single numeric value: the highest p<=alpha among the children of the parent node or the minimum of the p-values of the children if they are all greater than alpha
}
\description{
Given \eqn{k} child p-values, return the highest p-value below alpha. This is the same as respecting the global gating criteria of not continuing to test along any path where p > alpha, but doesn't do a local, parent-level adjustment. This is mostly useful for simulation and study purposes although with k=2 and low l  (for example, less than 5 or so) this will also tend to control the FWER.
}
\details{
The idea is to not adjust so that we can compare the consequences of not adjusting to adjusting
}
\examples{
local_min_p(c(0.01, 0.04, 0.10, 0.20)) ## should return .04
local_min_p(c(0.10, 0.20)) ## should return .04

}
