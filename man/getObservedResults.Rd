% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getObservedResults.R
\name{getObservedResults}
\alias{getObservedResults}
\title{Obtains the observed treatment effects for differential binding candidate regions.}
\usage{
getObservedResults(
  experiment,
  includeBetas = F,
  numCores = 1,
  clusterExists = F,
  verbose = T
)
}
\arguments{
\item{experiment}{A RangedSummarizedExperiment object describing the experiment.}

\item{includeBetas}{A boolean indicating whether or not to record the raw treatment effect (beta coefficient) and its corresponding p-value and FDR. By default, only considers the t-value of the treatment effect.}

\item{numCores}{The number of cores to use for computation. If \code{==1}, then computation is done in serial. Otherwise, it is run in parallel.}

\item{clusterExists}{Boolean indicating if a cluster has already been set up for parallelization. Passed as TRUE when called from \code{\link{permutationTest}} to avoid re-constructing the cluster.}

\item{verbose}{Boolean indicating whether or not to print out messages indicating progress.}
}
\description{
Obtains the observed treatment effects for differential binding candidate regions.
}
