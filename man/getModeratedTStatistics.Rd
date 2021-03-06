% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getModeratedTStatistics.R
\name{getModeratedTStatistics}
\alias{getModeratedTStatistics}
\title{Obtains the moderated t statistics.}
\usage{
getModeratedTStatistics(experiment)
}
\arguments{
\item{experiment}{A RangedSummarizedExperiment object describing the experiment.}
}
\description{
Obtains the moderated t statistics for an experiment following the same procedure as in limma. Called by \link{\code{getObservedResults()}}.
For details refer to "Linear Models and Empirical Bayes Methods for Assessing Differential Expression in Microarray Experiments" by Smyth (2004).
}
