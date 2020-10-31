#' Formats the results from getting read counts into an experiment for the permutation test. Also performs error checking.
#'
#' @param experiment A RangedSummarizedExperiment object from \code{windowCounts()} in \code{csaw}, or a comparable function. Must have a count assay in \code{assay(experiment)} and the corresponding GRanges in \code{rowRanges(experiment)}.
#' @param sampleSheet A data frame describing the samples in the experiment. At minimum, one column is needed to identify the experimental condition of each sample. A blocking variable and other covariates may also be included.
#' @param blackList A GRanges object containing genomic regions to exclude from further analysis.
#' @param expCondition A string indicating the name of the experimental condition variable in \code{sampleSheet}. Differential binding will be assessed across experimental condition. Permutations will be made across experimental condition. Only two groups (e.g. treatment and control) are currently permitted.
#' @param block A string indicating the name of the blocking variable in \code{sampleSheet}. A blocking variable can be used to identify matched pairs, if present, or a similar setup. In the permutation test, samples are only permuted within their respective blocks.
#' @param covariates A vector of strings indicating the names of any other covariates in \code{sampleSheet} to be controlled for in the GLS regressions.
#'
#' @export
#'
formatExperiment <- function(assay, rowRanges, sampleSheet, blacklist = NULL, expCondition, block = NULL, covariates = NULL){
    #Error checking
    if((class(assay) != "matrix")){
        stop("countAssay must be a matrix.")
    }

    if((class(rowRanges) != "GRanges")){
        stop("countAssay must be a GRanges object.")
    }

    if(!is.null(blacklist)){
        if((class(blacklist) != "GRanges")){
            stop("blacklist must be a GRanges object.")
        }
    }

    experiment <- SummarizedExperiment(assays=list(counts=assay),
                                       rowRanges=rowRanges,
                                       colData=sampleSheet)

    metadata(experiment)$blacklist <- blacklist
    metadata(experiment)$expCondition <- expCondition
    metadata(experiment)$block <- block
    metadata(experiment)$covariates <- covariates

    #Check that the experimental condition is in the sample sheet.
    if(!(metadata(experiment)$expCondition %in% colnames(colData(experiment)))){
        stop(paste("expCondition variable listed in metadata(experiment)$expCondition:", metadata(experiment)$expCondition, " is not in the sample sheet, colData(experiment)."))
    }

    #Check that the blocking variable is in the sample sheet.
    if(!is.null(metadata(experiment)$block)){
        if(!(metadata(experiment)$block %in% colnames(colData(experiment)))){
            stop(paste("Blocking variable listed in metadata(experiment)$block:", metadata(experiment)$block, " is not in the sample sheet, colData(experiment)."))
        }
    }

    #Check that the covariates are in the sample sheet.
    if(!is.null(metadata(experiment)$covariates)){
        if(!all(metadata(experiment)$covariates %in% colnames(colData(experiment)))){
            stop(paste("Covariates listed in metadata(experiment)$covariates are not all in the sample sheet, colData(experiment)."))
        }

        #Check if any covariates in the sample sheet say "expCondition".
        if(length(grep("expCondition", metadata(experiment)$covariates)) > 0){
            stop("No covariates used in metadata(experiment)$covariates can contain 'expCondition' in their name.")
        }
    }

    return(experiment)
}
