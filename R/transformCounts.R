#' Transforms the read counts associated with the experiment. Stores the result in a metadata field of the input \code{experiment}.
#'
#' @param experiment A RangedSummarizedExperiment object describing the experiment.
#' @param FUN A function implementing the transformation. Should take in a matrix (i.e. the count assay) and output a matrix of the same dimensions, preserving ordering of rows/columns. Default is to log-transform counts.
#'
#' @export
#'
transformCounts <- function(experiment, FUN = function(x) log(x + 0.5)){
  #Error check
  if((class(experiment) != "RangedSummarizedExperiment")){
    stop("experiment is not a RangedSummarizedExperiment object.")
  }

  metadata(experiment)$transformedCounts <- FUN(assay(experiment))
  return(experiment)
}
