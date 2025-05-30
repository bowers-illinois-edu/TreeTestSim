% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_blocks_simulations.R
\name{tau_norm}
\alias{tau_norm}
\alias{tau_norm_outliers}
\alias{tau_norm_covariate}
\alias{tau_norm_covariate_outliers}
\alias{tau_norm_covariate_cont}
\alias{tau_norm_covariate_levels}
\alias{tau_null}
\title{Causal Effect Functions}
\usage{
tau_norm(ybase, tau_sds, covariate)

tau_norm_outliers(ybase, tau_sds, covariate)

tau_norm_covariate(ybase, tau_sds, covariate)

tau_norm_covariate_outliers(ybase, tau_sds, covariate)

tau_norm_covariate_cont(ybase, tau_sds, covariate)

tau_norm_covariate_levels(ybase, tau_sds, covariate)

tau_null(ybase, tau_sds, covariate)
}
\arguments{
\item{ybase}{Is a vector of potential outcome to control}

\item{tau_sds}{is the number of sds of ybase to shift the distribution}

\item{covariate}{Contains information about covariates currently a character
name of a column in idat. It is NULL if not used. Mostly it is a vector the
same length as ybase.}
}
\value{
A vector of individual level causal effects (taus) that we will add
to ybase (potential outcome to control) to get y1var or potential outcome to
treatment.
}
\description{
These functions create individual level causal effects, \eqn{\tau_i}, that can be
combined with a potential outcome to control to create a potential outcome
to treatment.
}
\section{Functions}{
\itemize{
\item \code{tau_norm()}: A basic function with no outliers

\item \code{tau_norm_outliers()}: Draws from a Normal but also adds a few outliers.

\item \code{tau_norm_covariate()}: A basic function that specifies a tau_sds*2 size effect if cov>median(cov) and otherwise is a tau_sds/2 size effect. The idea is to keep the average individual effect size the same as other functions --- i.e. about tau_sds --- but to make a strong but simple relationship with a covariate.

\item \code{tau_norm_covariate_outliers()}: A basic function that specifies a tau_sds size effect if cov>median(cov) and otherwise is a tau_sds/4 size effect

\item \code{tau_norm_covariate_cont()}: A basic function that specifies a tau_sds size effect if cov>median(cov) and otherwise is a tau_sds/4 size effect

\item \code{tau_norm_covariate_levels()}: A basic function that specifies a tau_sds size effect that varies, randomly, by level of covariate. We envision this to be used with covariate that have relatively few levels into which the outcome values can be grouped.

\item \code{tau_null()}: Represent the no effects situation

}}
