#' Performs a permutation test to assess for differential binding.
#'
#' @param experiment A RangedSummarizedExperiment object describing the experiment.
#' @param permutations A matrix where each row corresponds to a permutation of the samples (i.e. permutation of the columns in the \code{counts} matrix).
#' @param effectType The type of effect of experimental condition on differential binding to look for, or how the p-value will be computed. Options include \code{"negative"} (p-value is proportion of null effects <= observed effect), \code{"positive"} (p-value is proportion of null effects >= observed effect), or \code{"two-sided"} (p-value is proportion of null effects in absolute value >= observed effect in absolute value). The option chosen here should correspond with how the lower and upper thresholds are chosen in \code{\link{getRegions}}.
#' @param includeBetas A boolean indicating whether or not to record the raw treatment effect (beta coefficient) and its corresponding p-value and FDR. By default, only considers the t-value of the treatment effect. Will be forced to be between 1 and \code{detectCores()}.
#' @param numCores The number of cores to use for computation. If \code{==1}, then computation is done in serial. Otherwise, it is run in parallel.
#' @param verbose Boolean indicating whether or not to print out messages indicating progress.
#'
#' @export
#'
permutationTest <- function(experiment, permutations = NA, effectType, includeBetas = F, numCores = 1, verbose = T){
  #Error check
  if((class(experiment) != "RangedSummarizedExperiment")){
    stop("experiment is not a RangedSummarizedExperiment object.")
  }

  if(class(permutations) != "matrix"){
    stop("permutations argument must be a matrix where each row is a permutation of the samples (i.e. permutation of columns in the count matrix).")
  }

  if(class(permutations) == "matrix"){
    if(ncol(permutations) != ncol(assay(experiment))){
      stop("permutations matrix must have the same number of columns as the matrix of read counts.")
    }
  }

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


  #Check effect type and thresholds chosen.
  if(metadata(experiment)$thresholdType == "percentile"){
    if((metadata(experiment)$lower != 0) & (metadata(experiment)$upper != 1) & ((0.5 - metadata(experiment)$lower) != (metadata(experiment)$upper - 0.5))){
      warning("For percentile threshold values, it is recommended that the lower and upper thresholds either be symmetrical (i.e. 0.5 - lower == upper - 0.5) or asymmetrical such that only one direction is considered (e.g. lower = 0.025, upper = 1, or lower = 0, upper = 0.95, etc.).")
    }
    if(((0.5 - metadata(experiment)$lower) == (metadata(experiment)$upper - 0.5)) & (effectType != "two-sided")){
      warning("For symmetric lower and upper thresholds, using a 'two-sided' effectType is suggested.")
    }
  }

  if(metadata(experiment)$thresholdType == "raw"){
    if((metadata(experiment)$lower > -Inf) & (metadata(experiment)$upper < Inf) & ((-metadata(experiment)$lower) != metadata(experiment)$upper)){
      warning("For raw threshold values, it is recommended that the lower and upper thresholds either be symmetrical (i.e. -lower == upper) or asymmetrical such that only one direction is considered (e.g. lower = 0, upper = Inf, or lower = -Inf, upper = 2, etc.).")
    }
    if(((-metadata(experiment)$lower) == (metadata(experiment)$upper)) & (effectType != "two-sided")){
      warning("For symmetric lower and upper thresholds, using a 'two-sided' effectType is suggested.")
    }
  }

  if(numCores < 1){
    message("Inputted numCores is less than 1, so it is being set equal to 1.")
    numCores <- 1
  }

  #Get the observed effects of experimental condition.
  if(verbose){
    message("Getting observed effects for the candidate regions.")
  }

  #Use boolean to record if cluster is created
  clusterExists <- F

  if(numCores > 1){
    if(verbose){
      message("Setting up cluster to run in parallel.")
    }

    #Set up cluster if it doesn't already exists
    if(numCores > detectCores()){
      numCores <- max(1, detectCores()-1)
    }

    if(verbose){
      message("Using ", numCores, " cores.")
    }

    cl <- makeCluster(numCores)
    registerDoParallel(cl)

    clusterExists <- T
  }

  experiment <- getObservedResults(experiment, includeBetas = includeBetas, numCores = numCores, clusterExists = clusterExists, verbose = verbose)

  #Check how many regressions failed to fit.
  missing <- is.na(metadata(experiment)$regions$tValue)
  if(sum(missing)/length(metadata(experiment)$regions) >= 0.01){
    warning(paste0(round(sum(missing)/length(metadata(experiment)$regions)*100, 2), "% of the models for the observed data failed to fit.
                   \nIf covariates are being used, they may cause issues with fitting regressions. You may want to use fewer or no covariates.
                   \nOtherwise, you may want to manually inspect your data."))
  }

  #Perform the permutation test, iterating over permutations.
  nullBetas <- nulltValues <- nullModtValues <- nullModtValuesFromBioVar <- list()
  for(i in 1:nrow(permutations)){
    if(verbose){
      message("Analyzing permutation ", i, " out of ", nrow(permutations))
    }

    #Copy over the experiment, permute the counts.
    experimentPerm <- experiment
    assay(experimentPerm) <- assay(experiment)[,permutations[i,]]
    if(!is.null(metadata(experiment)$transformedCounts)){
      metadata(experimentPerm)$transformedCounts <- metadata(experimentPerm)$transformedCounts[,permutations[i,]]
    }

    #Get new candidate regions using the same parameters as the observed results, get the corresponding null results.
    experimentPerm <- getRegions(experimentPerm, lower=metadata(experiment)$lower, upper = metadata(experiment)$upper,
                                 thresholdType = metadata(experiment)$thresholdType, minWidth = metadata(experiment)$minWidth,
                                 maxGap = metadata(experiment)$maxGap)
    experimentPerm <- getObservedResults(experimentPerm, includeBetas = includeBetas, numCores = numCores, clusterExists = clusterExists, verbose = verbose)

    #Append the results to those already obtained.
    nulltValues[[i]] <- metadata(experimentPerm)$regions$tValue
    nullModtValues[[i]] <- metadata(experimentPerm)$regions$modtValue
    nullModtValuesFromBioVar[[i]] <- metadata(experimentPerm)$regions$modtValueFromBioVar
    if(includeBetas){
      nullBetas[[i]] <- metadata(experimentPerm)$regions$beta
    }
  }

  #Stop the cluster if present
  if(clusterExists){
    stopCluster(cl)
  }

  #Reformat into vectors
  nulltValues <- unlist(nulltValues)
  nullModtValues <- unlist(nullModtValues)
  nullModtValuesFromBioVar <- unlist(nullModtValuesFromBioVar)
  if(includeBetas){
    nullBetas <- unlist(nullBetas)
  }

  #Report NA values (arise due to failure in generalized least squares fit).
  missing <- is.na(nulltValues)

  if(verbose){
    message(sum(missing), " out of ", length(nulltValues) ," null models failed to fit and are being excluded from analysis.")
  }

  if(sum(missing)/length(nulltValues) >= 0.01){
    warning(paste0(round(sum(missing)/length(nulltValues)*100, 2), "% of the null models failed to fit.
                   \nIf covariates are being used, they may cause issues with fitting regressions. You may want to use fewer or no covariates.
                   \nOtherwise, you may want to manually inspect your data."))
  }

  #Remove missing values.
  nulltValues <- nulltValues[!missing]
  nullModtValues <- nullModtValues[!missing]
  nullModtValuesFromBioVar <- nullModtValuesFromBioVar[!missing]
  if(includeBetas){
    nullBetas <- nullBetas[!missing]
  }

  #Obtain the p-values and FDRs for beta coefficients, for t values, and for moderated t values
  if(includeBetas){
    elementMetadata(metadata(experiment)$regions)[["betaPValue"]] <- sapply(metadata(experiment)$regions$beta, FUN=findPermPValue, nullValues=nullBetas, effectType = effectType)
    elementMetadata(metadata(experiment)$regions)[["betaFDR"]] <- p.adjust(metadata(experiment)$regions$betaPValue, method="BH")
  }

  elementMetadata(metadata(experiment)$regions)[["tPValue"]] <- sapply(metadata(experiment)$regions$tValue, FUN=findPermPValue, nullValues=nulltValues, effectType = effectType)
  elementMetadata(metadata(experiment)$regions)[["tFDR"]] <- p.adjust(metadata(experiment)$regions$tPValue, method="BH")

  elementMetadata(metadata(experiment)$regions)[["modtPValue"]] <- sapply(metadata(experiment)$regions$modtValue, FUN=findPermPValue, nullValues=nullModtValues, effectType = effectType)
  elementMetadata(metadata(experiment)$regions)[["modtFDR"]] <- p.adjust(metadata(experiment)$regions$modtPValue, method="BH")

  #Revert the following code when done with lm() rather than gls()
  message("Calculating permutation p-values and FDRs for moderated t statistics based on shrinking biological variance. NAs will be ignored in the calculations.")
  missing <- which(is.na(nullModtValuesFromBioVar) | is.nan(nullModtValuesFromBioVar))
  message("Ignoring ", length(missing), " NA or NaNs null values.")
  elementMetadata(metadata(experiment)$regions)[["modtFromBioVarPValue"]] <- sapply(metadata(experiment)$regions$modtValueFromBioVar, FUN=findPermPValue, nullValues=nullModtValuesFromBioVar, effectType = effectType, ignoreNAs = T)
  elementMetadata(metadata(experiment)$regions)[["modtFromBioVarFDR"]] <- p.adjust(metadata(experiment)$regions$modtFromBioVarPValue, method="BH")

  return(experiment)
}
