#' Obtains permutations of samples that both cross experimental condition and stay within blocks. This function is not particularly efficient for experiments with many samples. It may stall if the number of permutations requested is too large.
#
#' @param experiment A RangedSummarizedExperiment object describing the experiment.
#' @param permutationType A string indicating whether to obtain all ("all") permutations, or only "balanced" permutations. For instance, in a two treatment vs. two control setup, "balanced" permutations would only include those permutations where exactly one treatment sample is swapped with one control sample. "balanced" is generally recommended.
#' @param numRandomPerm An integer indicating the number of permutations to randomly sample if not "all" permutations or every "balanced" permutation should be returned. If \code{NA} (the default), every requested permutation is returned.
#' @param numCores The number of cores to use for computation. If \code{==1}, then computation is done in serial. Otherwise, it is run in parallel. Running in parallel is only useful for experiments with more than a few samples per experimental group.
#' @param verbose Boolean indicating whether or not to print out messages indicating progress.
#'
#' @export
#'
getPermutations <- function(experiment, permutationType = "balanced", numRandomPerm = NA, numCores = 1, verbose=T){
  #If there is a blocking factor, could use permute package to get all possible permutations respecting the blocking factor
  #if(!is.null(metadata(experiment)$block)){
  #  blocks <- colData(experiment)[,metadata(experiment)$block]
  #  permutations <- as.matrix(shuffleSet(length(blocks),
  #                                       control=how(blocks=blocks, complete=T, maxperm = factorial(length(blocks))),
  #                                       quietly = T))
  #
  #  #Add back in the observed order of samples.
  #  permutations <- rbind(1:length(blocks), permutations)
  #}

  #Generally, it's much more computationally feasible to generate permutations
  #manually by considering those permutations across experimental condition. The permute package doesn't
  #allow for this - it can restrict permutations to within groups but not across/between groups.

  #Matrix to store permutations
  permutations <- matrix(1:nrow(colData(experiment)), nrow=1, ncol=nrow(colData(experiment)))

  #Get the sample indices for each experimental condition
  groups <- unique(colData(experiment)[,metadata(experiment)$expCondition])
  if(length(groups) != 2){
    stop("There are more than or fewer than two groups associated with the identified experimental condition in the sampleSheet data frame. There must be exactly two groups.")
  }
  group1Members <- unname(which(colData(experiment)[,metadata(experiment)$expCondition] == groups[1]))
  group2Members <- unname(which(colData(experiment)[,metadata(experiment)$expCondition] == groups[2]))

  #Get all of the permutations that involve swapping i samples in one experimental condition with
  #i in the other. For instance, with 6 samples in each group, we can choose one from each to swap,
  #giving (6 choose 1)*(6 choose 1) = 36 permutations. We can choose 2 from each group, giving
  #(6 choose 2)*(6 choose 2) = 225 permutations, and so on for i = 1 to i=6. i = 0 refers to the observed data
  #which is included by default.
  for(i in 1:min(length(group1Members), length(group2Members))){
    group1Samples <- combn(x=group1Members, m=i, simplify = F)
    group2Samples <- combn(x=group2Members, m=i, simplify = F)
    swaps <- expand.grid(group1Samples, group2Samples)

    #Append the corresponding permutations to the permutation matrix
    for(j in 1:nrow(swaps)){
      permutationTemp <- 1:nrow(colData(experiment))
      permutationTemp[c(unlist(swaps[j, 1]),  unlist(swaps[j,2]))] <- c(unlist(swaps[j, 2]),  unlist(swaps[j,1]))
      permutations <- rbind(permutations, unname(permutationTemp))
    }
  }

  if(numCores < 1){
    message("Inputted numCores is less than 1, so it is being set equal to 1.")
    numCores <- 1
  }


  #Run in parallel if requested.
  if(numCores > 1){
    #Set up cluster.
    if(numCores > detectCores()){
      numCores <- max(1, detectCores()-1)
    }

    if(verbose){
      message("Using ", numCores, " cores. Checking which permutations are valid.")
    }

    cl <- makeCluster(numCores)
    registerDoParallel(cl)

    keep <- foreach(i=1:nrow(permutations), .combine=c) %dopar% {
      checkPermutation(experiment=experiment, permutation=permutations[i,], permutationType = permutationType)
    }

    stopCluster(cl)

  }

  #Else, run in serial
  else{
    if(verbose){
      message("Running in serial. Checking which permutations are valid.")
    }

    #Record which permutations are valid (based on the given permutationType) and keep them.
    keep <- apply(permutations, FUN=checkPermutation, MARGIN = 1, experiment=experiment, permutationType=permutationType)
  }

  permutations <- permutations[keep,]


  #If only a random subset was requested, get a random subset.
  if(!is.na(numRandomPerm)){
    if(numRandomPerm >= nrow(permutations)){
      message("Desired number of random permutations is at least as large as the total number of valid permutations. Returning all valid permutations.")
      return(permutations)
    }
    permutations <- permutations[sort(sample(1:nrow(permutations), size=numRandomPerm)),]
  }



  return(permutations)
}
