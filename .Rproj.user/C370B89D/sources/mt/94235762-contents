#' Gets the treatment effect on read counts for a single region.
#'
#' @param experiment A RangedSummarizedExperiment object describing the experiment.
#' @param region A single GRanges region.
#' @param verbose Boolean indicating whether or not to print out messages indicating progress. If this function is being used in parallel, \code{verbose} is set to \code{FALSE} to avoid issues when using \code{tryCatch()}.
#'
#' @export
#'
getSingleRegionEffect <- function(experiment, region, verbose = T){
  #Grab the candidate region's corresponding counts.
  regionBinIndices <- metadata(experiment)$regionsBins[region$index,"startBin"]:metadata(experiment)$regionsBins[region$index,"endBin"]
  if(!is.null(metadata(experiment)$transformedCounts)){
    regionDat <- metadata(experiment)$transformedCounts[regionBinIndices,,drop=F] #drop=F maintains matrix format when only grabbing one bin
  }
  else{
    regionDat <- assay(experiment)[regionBinIndices,,drop=F] #drop=F maintains matrix format when only grabbing one bin
  }
  rownames(regionDat) <- regionBinIndices
  colnames(regionDat) <- paste("Sample", 1:ncol(regionDat))

  #Need to reshape data into long format with columns: count, bin, condition, block, other covariatess
  regionDatLong <- melt(regionDat) #Lists all counts from first sample, then second, etc.
  colnames(regionDatLong) <- c("bin", "sample", "count")

  #Add the other variables.
  regionDatLong$expCondition <- as.factor(rep(colData(experiment)[,metadata(experiment)$expCondition], each=nrow(regionDat)))
  if(!is.null(metadata(experiment)$block)){
    regionDatLong$block <- as.factor(rep(colData(experiment)[,metadata(experiment)$block], each=nrow(regionDat)))
  }
  if(!is.null(metadata(experiment)$covariates)){
    regionDatLong[,metadata(experiment)$covariates] <- as.matrix(rep(colData(experiment)[,metadata(experiment)$covariates], each=nrow(regionDat)))
  }

  #Relevel experimental condition to ensure consistency with region-finding treatment effect estimate.
  regionDatLong$expCondition <- relevel(regionDatLong$expCondition, ref=toString(unique(colData(experiment)[,metadata(experiment)$expCondition])[1]))

  #Convert sample to factor.
  regionDatLong$sample <- as.factor(regionDatLong$sample)

  #Fit GLS.
  fit.gls <- tryCatch({
    if(length(unique(regionDatLong$bin)) > 1){
      #Bin index must be treated as a factor when used a covariate, but as numerical for the correlation.
      if(!is.null(metadata(experiment)$block)){
        formula <- paste("count ~ factor(bin) + block + expCondition + ", paste(metadata(experiment)$covariates, collapse=" + "), " - 1",sep = "")
        formula <- as.formula(gsub("+  -", "-", formula, fixed=T)) #In case there are no additional covariates to include
        #fit <- gls(formula, data=regionDatLong, correlation = corAR1(form = ~bin|sample), control = list(returnObject=T))
        fit <- lm(formula, data=regionDatLong)
      }
      else{
        formula <- paste("count ~ factor(bin) + expCondition + ", paste(metadata(experiment)$covariates, collapse=" + "), " - 1",sep = "")
        formula <- as.formula(gsub("+  -", "-", formula, fixed=T)) #In case there are no additional covariates to include
        fit <- gls(formula, data=regionDatLong, correlation = corAR1(form = ~bin|sample), control = list(returnObject=T))
      }
    }
    else{ #If the data correspond to only one bin
      if(!is.null(metadata(experiment)$block)){
        formula <- paste("count ~ block + expCondition + ", paste(metadata(experiment)$covariates, collapse=" + "), " - 1",sep = "")
        formula <- as.formula(gsub("+  -", "-", formula, fixed=T)) #In case there are no additional covariates to include
        fit <- gls(formula, data=regionDatLong, correlation = corAR1(form = ~bin|sample), control = list(returnObject=T))
      }
      else{
        formula <- paste("count ~ expCondition + ", paste(metadata(experiment)$covariates, collapse=" + "), " - 1",sep = "")
        formula <- as.formula(gsub("+  -", "-", formula, fixed=T)) #In case there are no additional covariates to include
        fit <- gls(formula, data=regionDatLong, correlation = corAR1(form = ~bin|sample), control = list(returnObject=T))
      }
    }
  }, error = function(cond) {
    #Can't include message() calls when running in parallel (tryCatch() doesn't work if they're included).
    #if(verbose){
    #  message(paste0("Original error message for region:"))
    #  message(cond)
    #  message("\nReturning NA for this region.")
    #}


    #message("Returning a treatment effect without including a correlation structure instead for this region.")
    ## Choose a return value in case of error
    #if(length(unique(regionDatLong$bin)) > 1){
    #  fit <- gls(count ~ factor(bin) + block + expCondition -1, data=regionDatLong, control = list(returnObject=T))
    #}
    #else{
    #  fit <- gls(count ~ block + expCondition -1, data=regionDatLong, control = list(returnObject=T))
    #}
    #return(fit)
    return(NA)
  })

  if((class(fit.gls) != "gls") & (class(fit.gls) != "lm")){
    return(rep(NA, 9))
  }

  #Record information for obtaining moderated t statistics as done in limma
  #https://www.degruyter.com/view/journals/sagmb/6/1/article-sagmb.2007.6.1.1252.xml.xml
  expConditionCoef <- grep("expCondition", names(fit.gls$coefficients), value=T)
  residDF <- nrow(regionDatLong) - length(coefficients(fit.gls))
  residSE <- summary(fit.gls)$sigma #For some reason, not equal to sqrt(sum(fit.gls$residuals^2)/residDF) but close (for lm it is).

  #Revert the following code when done with lm() rather than gls()
  #betaCoef <- summary(fit.gls)$tTable[expConditionCoef,c("Value")]
  betaCoef <- summary(fit.gls)$coefficients[expConditionCoef,c("Estimate")]
  #betaCoefSE <- summary(fit.gls)$tTable[expConditionCoef,c("Std.Error")]
  betaCoefSE <- summary(fit.gls)$coefficients[expConditionCoef,c("Std. Error")]

  #Revert the following code when done with lm() rather than gls()
  #Extract the t-value for the effect associated with the experimental condition.
  #expConditionEffect <- summary(fit.gls)$tTable[expConditionCoef,c("t-value")]
  expConditionEffect <- summary(fit.gls)$coefficients[expConditionCoef,c("t value")]
  #anova(fit.gls) would return sequential sum of squares, at least when using lm. Don't want that.
  #Want sum of squares within and across samples, based on residuals (e.g. after fitting model).
  residDataFrame <- data.frame(regResiduals=resid(fit.gls), sample=regionDatLong$sample, block=regionDatLong$block)
  anovaResults <- aov(regResiduals ~ sample, data=residDataFrame) #The "Residuals" column in this output refers to residuals from ANOVA, not from the regression.
  withinSS <- summary(anovaResults)[[1]]["Residuals","Sum Sq"]
  withinDF <- summary(anovaResults)[[1]]["Residuals","Df"] - length(coef(fit)) #Subtract off number of coefficients in model
  betweenSS <- summary(anovaResults)[[1]]["sample","Sum Sq"]
  betweenDF <- summary(anovaResults)[[1]]["sample","Df"]


  #Return all of the results
  regressionResults <- c(expConditionEffect, residDF, residSE, betaCoef, betaCoefSE, withinSS, withinDF, betweenSS, betweenDF)
  names(regressionResults) <- c("tValue", "residDF", "residSE", "betaCoef", "betaCoefSE", "withinSS", "withinDF", "betweenSS", "betweenDF")

  return(regressionResults)
}

