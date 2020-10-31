#' Obtains the moderated t statistics.
#'
#' Obtains the moderated t statistics for an experiment following the same procedure as in limma. Called by \link{\code{getObservedResults()}}.
#' For details refer to "Linear Models and Empirical Bayes Methods for Assessing Differential Expression in Microarray Experiments" by Smyth (2004).
#'
#' @param experiment A RangedSummarizedExperiment object describing the experiment.
#'
#' @export
#'
getModeratedTStatistics <- function(experiment){
  #Follows the procedures in https://www.degruyter.com/view/journals/sagmb/6/1/article-sagmb.2007.6.1.1252.xml.xml

  #Record results with missing values, may not actually be necessary
  #missing <- unique(which(is.na(metadata(experiment)$modtInfo), arr.ind=T)[,1])

  #Obtain z_g = log s_g^2
  z <- log((metadata(experiment)$modtInfo[,"residSE"])^2)

  #Obtain e_g = z_g - digamma(d_g/2) + log(d_g/2)
  d <- metadata(experiment)$modtInfo[,"residDF"]
  e <- z - digamma(d/2) + log(d/2)
  eBar <- mean(e, na.rm=T)

  #Get number of regions with a t statistic
  n <- sum(!is.na(metadata(experiment)$regions$tValue))

  #Need to get estimate of d0 using Newton's method
  #Right hand side (RHS) is:
  RHS <- mean( ((e - eBar)^2)*(n/(n-1)) - trigamma(d/2), na.rm=T)

  if(RHS < 0){
    #If this occurs, d0 is set to Inf, s0 to sqrt(exp(eBar)) and sTildeSquare is just equal to s0^2 later
    d0 <- Inf
    s0 <- sqrt(exp(eBar))
  }
  else{
    #Run Newton's method (see Appendix of the paper referenced above)
    #Need to solve trigamma(d0/2) = RHS
    y <- 0.5 + 1/RHS
    delta <- trigamma(y)*(1 - trigamma(y)/RHS)/psigamma(y, deriv = 2)

    while((-delta/y) > (10^(-6))){
      y <- y + delta
      delta <- trigamma(y)*(1 - trigamma(y)/RHS)/psigamma(y, deriv = 2)
    }

    d0 <- y*2

    #Get s0 estimate
    s0 <- sqrt( exp( eBar + digamma(d0/2) - log(d0/2) ) )
  }

  #Calculate the moderated t values
  sgSquare <- metadata(experiment)$modtInfo[,"residSE"]^2
  if(RHS < 0){
    sgTildeSquare <- s0^2
  }
  else{
    sgTildeSquare <- (d0*(s0^2) + d*sgSquare)/(d0 + d)
  }
  vg <- (metadata(experiment)$modtInfo[,"betaCoefSE"]^2)/sgSquare
  modtValue <- metadata(experiment)$modtInfo[,"betaCoef"]/(sqrt(sgTildeSquare)*sqrt(vg))

  #cor(modtValue, metadata(experiment)$regions$tValue) #Around 0.86 in Simulated 2x2

  #Return the moderated t values
  return(modtValue)
}
