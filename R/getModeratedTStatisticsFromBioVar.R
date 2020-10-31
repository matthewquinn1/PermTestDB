#' Obtains the moderated t statistics by shrinking only the biological variance.
#'
#' Assumes that lm is used rather than gls when calculating standard errors.
#'
#' @param experiment A RangedSummarizedExperiment object describing the experiment.
#' @param verbose Boolean indicating whether or not to print out messages indicating progress.
#'
#' @export
#'
getModeratedTStatisticsFromBioVar <- function(experiment, verbose = T){
  #Shrink the biological variance
  modBioVar <- squeezeVar(var=metadata(experiment)$modtInfo[,"betweenSS"]/metadata(experiment)$modtInfo[,"betweenDF"],
                          df = metadata(experiment)$modtInfo[,"betweenDF"])
  modtInfo <- cbind(metadata(experiment)$modtInfo,
                    rep(modBioVar$df.prior, nrow(metadata(experiment)$modtInfo)),
                    modBioVar$var.post)
  colnames(modtInfo) <-  c(colnames(metadata(experiment)$modtInfo), "modBioVarPriorDF", "modBioVar")
  metadata(experiment)$modtInfo <- modtInfo #Temporary

  #Go through each region, recalculate the standard error (assumes lm), get the new t statistic.
  #The size of the regression model or number of coefficients (p) has already been subtracted off
  #from withinDF by getSingleRegionEffect(). Thus, the p doesn't need to be accounted for below.
  #Addition of 1 is simply to match up with formulas in my notes. 1 was already subtracted from the DFs
  #in getSingleRegionEffect() so it's added back here only to get canceled when getting newResidVar.
  nStar <- metadata(experiment)$modtInfo[,"withinDF"] + metadata(experiment)$modtInfo[,"betweenDF"] + metadata(experiment)$modtInfo[,"modBioVarPriorDF"] + 1
  kStar <- metadata(experiment)$modtInfo[,"betweenDF"] + metadata(experiment)$modtInfo[,"modBioVarPriorDF"] + 1
  tauSquare <- metadata(experiment)$modtInfo[,"withinSS"]/metadata(experiment)$modtInfo[,"withinDF"]
  #sigmaSquare <- metadata(experiment)$modtInfo[,"betweenSS"]/metadata(experiment)$modtInfo[,"betweenDF"]
  sigmaSquare <- modBioVar$var.post
  newResidVar <- ((nStar - kStar)/(nStar - 1))*tauSquare + ((kStar - 1)/(nStar - 1))*sigmaSquare #Number of regression coefficients, p, already accounted for by "withinDF" above.

  #For reference, the original MSE (the mean squared residual) before shrinkage of biological variance.
  #This is not the same as
  #(metadata(experiment)$modtInfo[,"withinSS"] +  metadata(experiment)$modtInfo[,"betweenSS"])/( metadata(experiment)$modtInfo[,"withinDF"] + metadata(experiment)$modtInfo[,"betweenDF"])
  #(i.e. withinSS + betweenSS/(withinDF + betweenSS))
  #because, compared to the original regression, one more DF is lost when performing anova on the residuals themselves
  #with respect to sample (which is not the original regression).
  #Instead, it's equal to withinSS + betweenSS/(withinDF + betweenSS + 1).
  #originalResidVar <- (metadata(experiment)$modtInfo[,"withinSS"] +  metadata(experiment)$modtInfo[,"betweenSS"])/( metadata(experiment)$modtInfo[,"withinDF"] + metadata(experiment)$modtInfo[,"betweenDF"]+1)
  originalResidVar <- metadata(experiment)$modtInfo[,"residSE"]^2
  #summary(newResidVar/originalResidVar)

  #First method for getting moderated t statistics: Iterate over regions, re-estimate se(beta.hat)
  #Much slower than second approach below. Both give the same results.
  #results <- rep(NA, length(metadata(experiment)$regions))

  #Iterate over regions, getting the standard error and t statistics
  #for(i in 1:length(metadata(experiment)$regions)){
  #  if(i %% 1000 == 0){
  #    if(verbose){
  #      message(paste("Getting moderated statistics for region", i, "out of", length(metadata(experiment)$regions)))
  #    }
  #  }

    #Grab a region
  #  region <- metadata(experiment)$regions[i]

    #Get the data for this region
  #  regionBinIndices <- metadata(experiment)$regionsBins[region$index,"startBin"]:metadata(experiment)$regionsBins[region$index,"endBin"]
  #  if(!is.null(metadata(experiment)$transformedCounts)){
  #    regionDat <- metadata(experiment)$transformedCounts[regionBinIndices,,drop=F] #drop=F maintains matrix format when only grabbing one bin
  #  }
  #  else{
  #    regionDat <- assay(experiment)[regionBinIndices,,drop=F] #drop=F maintains matrix format when only grabbing one bin
  #  }
  #  rownames(regionDat) <- regionBinIndices
  #  colnames(regionDat) <- paste("Sample", 1:ncol(regionDat))

    #Need to reshape data into long format with columns: count, bin, condition, block, other covariatess
  #  regionDatLong <- melt(regionDat) #Lists all counts from first sample, then second, etc.
  #  colnames(regionDatLong) <- c("bin", "sample", "count")

    #Add the other variables.
  #  regionDatLong$expCondition <- as.factor(rep(colData(experiment)[,metadata(experiment)$expCondition], each=nrow(regionDat)))
  #  if(!is.null(metadata(experiment)$block)){
  #    regionDatLong$block <- as.factor(rep(colData(experiment)[,metadata(experiment)$block], each=nrow(regionDat)))
  #  }
  #  if(!is.null(metadata(experiment)$covariates)){
  #    regionDatLong[,metadata(experiment)$covariates] <- as.matrix(rep(colData(experiment)[,metadata(experiment)$covariates], each=nrow(regionDat)))
  #  }

    #Relevel experimental condition to ensure consistency with region-finding treatment effect estimate.
  #  regionDatLong$expCondition <- relevel(regionDatLong$expCondition, ref=toString(unique(colData(experiment)[,metadata(experiment)$expCondition])[1]))

    #Convert sample to factor.
  #  regionDatLong$sample <- as.factor(regionDatLong$sample)

    #Revert or adjust the following code when done with lm() rather than gls(). Assumes blocking variable is present.
    #Should look more like code in getSingleRegionEffect()
  #  if(!is.null(metadata(experiment)$block)){
  #    formula <- paste("count ~ factor(bin) + block + expCondition + ", paste(metadata(experiment)$covariates, collapse=" + "), " - 1",sep = "")
  #    formula <- as.formula(gsub("+  -", "-", formula, fixed=T)) #In case there are no additional covariates to include
  #  }

    #Get the design Matrix
  #  X <- model.matrix(formula, regionDatLong)
  #  expConditionDesign <- grep("expCondition", colnames(X), value=T)

    #Get the standard error for the treatment coefficient
  #  newSEBeta <- sqrt(newResidVar[i]*diag(solve(t(X)%*%X)))[expConditionDesign]
  #  beta <- metadata(experiment)$modtInfo[,"betaCoef"][i] #for some reason, extracting the beta from "region" object directly doesn't work.
  #  #cat("Region", i, "with new se", sqrt(newResidVar[i]*diag(solve(t(X)%*%X)))[expConditionDesign], "and beta", beta)
  #  if(!is.null(beta/newSEBeta)){
  #    results[i] <- beta/newSEBeta
  #  }
  #  else{
  #    results[i] <- NA
  #  }
  #}

  #Revised approach to getting moderated t statistics (avoid the matrix calculations)

  #Record the sqrt(var.hat(y))/sqrt(var.hat*(y)) where var.hat(y) is the original estimate of the variance
  #of y (or the total error in the regression model) and var.hat*(y) is the posterior estimate of this
  #variance after shrinking the biological variance.
  changeInSE <- sqrt(newResidVar/originalResidVar)
  newSE <- changeInSE*metadata(experiment)$modtInfo[,"betaCoefSE"]
  resultsRevised <- metadata(experiment)$modtInfo[,"betaCoef"]/newSE
  #summary(changeInSE)
  #summary(metadata(experiment)$regions$tValue/resultsRevised)
  #summary(metadata(experiment)$regions$tValue/metadata(experiment)$regions$modtValueFromBioVar)

  #results from above are stored in metadata(experiment)$regions$modtValueFromBioVar
  #resultsRevised - metadata(experiment)$regions$modtValueFromBioVar

  return(resultsRevised)
  }
