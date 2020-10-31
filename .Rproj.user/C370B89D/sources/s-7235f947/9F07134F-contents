#' Obtains the observed treatment effects for differential binding candidate regions.
#'
#' @param experiment A RangedSummarizedExperiment object describing the experiment.
#' @param includeBetas A boolean indicating whether or not to record the raw treatment effect (beta coefficient) and its corresponding p-value and FDR. By default, only considers the t-value of the treatment effect.
#' @param numCores The number of cores to use for computation. If \code{==1}, then computation is done in serial. Otherwise, it is run in parallel.
#' @param clusterExists Boolean indicating if a cluster has already been set up for parallelization. Passed as TRUE when called from \code{\link{permutationTest}} to avoid re-constructing the cluster.
#' @param verbose Boolean indicating whether or not to print out messages indicating progress.
#'
#' @export
#'
getObservedResults <- function(experiment, includeBetas = F, numCores = 1, clusterExists = F, verbose = T){
  #Error check
  if((class(experiment) != "RangedSummarizedExperiment")){
    stop("experiment is not a RangedSummarizedExperiment object.")
  }

  if(numCores < 1){
    message("Inputted numCores is less than 1, so it is being set equal to 1.")
    numCores <- 1
  }

  #Get the observed treatment effects. Can record both the beta and the t-value.
  results <- matrix(NA, nrow = length(metadata(experiment)$regions), ncol= 9)

  #Procedure in serial
  if(numCores == 1){
    for(i in 1:length(metadata(experiment)$regions)){
      if(i %% 1000 == 0){
        if(verbose){
          message(paste("Getting effect for region", i, "out of", length(metadata(experiment)$regions)))
        }
      }
      #If running in serial, verbose can be TRUE. If running in parallel, verbose must be false to avoid issues with tryCatch().
      results[i,] <- getSingleRegionEffect(experiment, metadata(experiment)$regions[i], verbose=((numCores == 1) & verbose))
    }
  }

  #Procedure in parallel
  else{
    if(verbose & !(clusterExists)){
      message("Setting up cluster and running in parallel.")
    }

    #Set up cluster if it doesn't already exists
    if((numCores > detectCores()) & !(clusterExists)){
      numCores <- max(1, detectCores()-1)
    }

    if(!clusterExists){
      if(verbose){
        message("Using ", numCores, " cores.")
      }

      cl <- makeCluster(numCores)
      registerDoParallel(cl)
    }

    #Get treatment effects in parallel
    results <- foreach(i=1:length(metadata(experiment)$regions), .combine=rbind) %dopar% {
      #If running in serial, verbose can be TRUE. If running in parallel, verbose must be false to avoid issues with tryCatch().
      getSingleRegionEffect(experiment, metadata(experiment)$regions[i], verbose = F)
    }

    if(!clusterExists){
      stopCluster(cl)
    }
  }

  #Name the results columns in accordance with what getSingleRegionEffect() feeds back.
  #colnames(results) <- c("tValue", "residDF", "residSE", "betaCoef", "betaCoefSE")
  colnames(results) <- c("tValue", "residDF", "residSE", "betaCoef", "betaCoefSE", "withinSS", "withinDF", "betweenSS", "betweenDF")

  #Record the results in the dbRanges attribute of the experiment.
  if(includeBetas){
    elementMetadata(metadata(experiment)$regions)[["beta"]] <- results[,"betaCoef"]
    elementMetadata(metadata(experiment)$regions)[["tValue"]] <- results[,"tValue"]
  }
  else{
    elementMetadata(metadata(experiment)$regions)[["tValue"]] <- results[,"tValue"]
  }

  #Also record information relevant for getting the moderated t statistic.
  metadata(experiment)$modtInfo <- results[,c("residDF", "residSE", "betaCoef", "betaCoefSE", "withinSS", "withinDF", "betweenSS", "betweenDF")]

  #Get and store the moderated t values based on shrinking var(beta.hat)
  #message("Getting moderated t stats")
  modtValue <- getModeratedTStatistics(experiment)
  metadata(experiment)$modtInfo <-  cbind(metadata(experiment)$modtInfo, modtValue)
  elementMetadata(metadata(experiment)$regions)[["modtValue"]]<- modtValue

  #Revert the following code when done with lm() rather than gls(). Or adjust it - this function assumes
  #an lm is used specifically when calculating the standard error.
  #Get and store the moderated t values based on shrinking biological variance only, and recalculating var(beta.hat)
  modtValueFromBioVar <- getModeratedTStatisticsFromBioVar(experiment, verbose)
  metadata(experiment)$modtInfo <- cbind(metadata(experiment)$modtInfo, modtValueFromBioVar)
  elementMetadata(metadata(experiment)$regions)[["modtValueFromBioVar"]] <- modtValueFromBioVar


  return(experiment)
}
