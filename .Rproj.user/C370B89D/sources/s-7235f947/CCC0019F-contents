#' Checks if a permutation is a valid one to include in the permutation test.
#
#' @param experiment A RangedSummarizedExperiment object describing the experiment.
#' @param permutation A vector reflecting a permutation of the samples in the experiment (the original ordering is that which corresponds with the rows of the \code{sampleSheet} data frame in \code{experiment}).
#' @param permutationType A string indicating whether to obtain all ("all") permutations or only "balanced" permutations. For instance, in a two treatment vs. two control setup, "balanced" permutations would only include those permutations where exactly one treatment sample is swapped with one control sample. "balanced" is generally recommended.
#'
#' @export
#'
checkPermutation <- function(experiment, permutation, permutationType = "balanced"){
  if(length(permutation) != nrow(colData(experiment))){
    stop("Number of elements in the suggested permutation does not match the number of samples implied by the sampleSheet data frame.")
  }

  #Record the number of samples that should be swapped from either experimental group.
  #If permutationType is "all", any number of swaps up to the size of the smaller
  #experimental group is allowed. If "balanced", then the number of swaps is allow to be either the floor
  #or the ceiling of half the size of the smaller experimental group (excluding 0 if an experimental group
  #happens to only have 1 sample).
  sizeSmallerExpCond <- min(table(colData(experiment)[,metadata(experiment)$expCondition]))

  numSwap <- 0:sizeSmallerExpCond
  if(permutationType == "balanced"){
    numSwap <- unique(c(floor(sizeSmallerExpCond/2), ceiling(sizeSmallerExpCond/2)))
    if(sizeSmallerExpCond == 1){
      numSwap <- 1
    }
  }

  #Check how many swaps have occured across experimental condition.
  numCondDisagreement <- sum(colData(experiment)[,metadata(experiment)$expCondition] != colData(experiment)[,metadata(experiment)$expCondition][permutation])
  if(!((numCondDisagreement/2) %in% numSwap)){
    return(FALSE)
  }


  #Check that no swaps have been made outside of the blocking factor.
  if(!is.null(metadata(experiment)$block)){
    numBlockDisagreement <- sum(colData(experiment)[,metadata(experiment)$block] != colData(experiment)[,metadata(experiment)$block][permutation])
    if(numBlockDisagreement != 0){
      return(FALSE)
    }
  }

  #Check that no swaps within experimental condition.
  #e.g. We wouldn't want a permutation that swaps two treatment samples. We only want to
  #consider permutations that only destroy the association between experimental condition and binding affinity.
  mismatches <- which(permutation != 1:nrow(colData(experiment)))

  for(i in mismatches){
    if(colData(experiment)[,metadata(experiment)$expCondition][i] == colData(experiment)[,metadata(experiment)$expCondition][permutation][i]){
      return(FALSE)
    }
  }

  #If no rule has been violated, the permutation is valid.
  return(TRUE)
}
