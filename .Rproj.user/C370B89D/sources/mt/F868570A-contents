#' Finds candidate regions to assess for differential binding.
#'
#' @param experiment A RangedSummarizedExperiment object describing the experiment.
#' @param lower A lower threshold. Regions with treatment effects below this threshold will analyzed for possible differential binding. In general, either \code{lower} and \code{upper} should provide symmetry (equally large effects in either direction are considered) or they should be asymmetrical such that only negative or only positive effects are considered.
#' @param upper An upper threshold. Regions with treatment effects above this threshold will analyzed for possible differential binding. In general, either \code{lower} and \code{upper} should provide symmetry (equally large effects in either direction are considered) or they should be asymmetrical such that only negative or only positive effects are considered.
#' @param thresholdType Possible values include \code{"percentile"} and \code{"raw"}. If \code{"percentile"}, the \code{lower} and \code{upper} values will be considered to be percentile cutoffs for treatment effect. If \code{"raw"}, they will instead be considered to be raw value cutoffs for treatment effects.
#' @param minWidth Candidate regions narrow than this width will be excluded from analysis. It is strongly recommended to make this larger than the bin/width width used for generating read counts. If left \code{NA} (the default), it will be set to to the median width of the \code{ranges} attribute of \code{experiment} plus one base pair.
#' @param maxGap Candidate regions for differential binding analysis containing gaps larger than this value (in base pairs) will be split into multiple regions, which will be analyzed separately.
#' @param verbose Boolean indicating whether or not to print out messages indicating progress.
#'
#' @export
#'
getRegions <- function(experiment, lower = 0, upper = 1, thresholdType = "percentile", minWidth = NA, maxGap = 1500, verbose = T){
  #Error check
  if((class(experiment) != "RangedSummarizedExperiment")){
    stop("experiment is not a RangedSummarizedExperiment object.")
  }

  if(is.null(metadata(experiment)$expCondition)){
    stop("An experimental condition (one of the columns of colData(experiment)) must be entered in metadata(experiment)$expCondition")
  }

  if(!(metadata(experiment)$expCondition %in% colnames(colData(experiment)))){
    stop("An experimental condition (one of the columns of colData(experiment)) must be entered in metadata(experiment)$expCondition")
  }

  if(is.na(minWidth)){
    minWidth <- median(width(rowRanges(experiment))) + 1
  }
  else{
    if(minWidth <= median(width(rowRanges(experiment)))){
      warning("The minimum region width is no more than the median bin/window width used for generating read counts. This may lead to errors later.")
    }
  }

  #Check lower and upper thresholds.
  if(thresholdType == "percentile"){
    if((lower != 0) & (upper != 1) & ((0.5 - lower) != (upper - 0.5))){
      warning("For percentile threshold values, it is recommended that the lower and upper thresholds either be symmetrical (i.e. 0.5 - lower == upper - 0.5) or asymmetrical such that only one direction is considered (e.g. lower = 0.025, upper = 1, or lower = 0, upper = 0.95, etc.).")
    }
  }
  if(thresholdType == "raw"){
    if((lower > -Inf) & (upper < Inf) & ((-lower) != upper)){
      warning("For raw threshold values, it is recommended that the lower and upper thresholds either be symmetrical (i.e. -lower == upper) or asymmetrical such that only one direction is considered (e.g. lower = 0, upper = Inf, or lower = -Inf, upper = 2, etc.).")
    }
  }

  #Take difference in sample means across experimental condition by bin.
  groups <- unique(colData(experiment)[,metadata(experiment)$expCondition])
  if(length(groups) != 2){
    stop("There are more than or fewer than two groups associated with the identified experimental condition in the sampleSheet data frame. There must be exactly two groups.")
  }
  group1Members <- unname(which(colData(experiment)[,metadata(experiment)$expCondition] == groups[1]))
  group2Members <- unname(which(colData(experiment)[,metadata(experiment)$expCondition] == groups[2]))

  message(paste0("Treatment effect taken as: ", metadata(experiment)$expCondition, "==", groups[2]," - ",  metadata(experiment)$expCondition, "==", groups[1]))
  diffMean <- rowMeans(assay(experiment)[,group2Members]) - rowMeans(assay(experiment)[,group1Members])

  #Adjust the threshold values if they are given as percentiles.
  if(thresholdType == "percentile"){
    lowerAdjusted <- quantile(diffMean, lower)
    upperAdjusted <- quantile(diffMean, upper)
  }

  #Find pertinent bins/windows based on the difference in sample means.
  outsideThresholds <- which((diffMean < lowerAdjusted) | (diffMean > upperAdjusted))

  #Find region starting positions
  starts <- outsideThresholds[which(!((outsideThresholds - 1) %in% outsideThresholds))]

  #Find region ending positions
  ends <- outsideThresholds[which(!((outsideThresholds + 1) %in% outsideThresholds))]

  #Check where any gaps larger than maxGap are.
  #This vector records the bin index where the gap BEGINS, not finishes.
  gaps <- which(start(rowRanges(experiment))[-1] > end(rowRanges(experiment))[-length(rowRanges(experiment))] + maxGap)

  #Also include chromosomal gaps that happen to not be caught.
  chrGaps <- which(as.vector(seqnames(rowRanges(experiment))[-1] != seqnames(rowRanges(experiment))[-length(rowRanges(experiment))]))
  gaps <- sort(c(gaps, chrGaps))

  #Iterate over the different gaps found, break up regions accordingly.
  i=1
  lastRegionEndBin <- ends[length(ends)]
  lastRegionStartBin <- starts[length(starts)]


  while((i <= length(gaps)) & (gaps[i] < lastRegionEndBin)){
    gap <- gaps[i]

    #Find the region that would contain the gap.
    indexInStartsOrEnds <- match(TRUE, ends - gap > 0)

    #If the region contains the gap, split it up.
    if((starts[indexInStartsOrEnds] <= gap) & (ends[indexInStartsOrEnds] > gap)){
      if(gap != lastRegionStartBin){
        starts <- c(starts[1:indexInStartsOrEnds], gap+1, starts[(indexInStartsOrEnds+1):length(starts)])
      }
      else{
        starts <- c(starts[1:indexInStartsOrEnds], gap+1)
      }
      ends <- c(ends[1:indexInStartsOrEnds-1], gap, ends[(indexInStartsOrEnds):length(ends)])
    }

    i <- i + 1
  }

  #Construct a genomic ranges object corresponding to the candidate regions.
  chromosomes <- seqnames(rowRanges(experiment))[starts]
  chromosomesCheck <- seqnames(rowRanges(experiment))[ends]
  regionStart <- start(rowRanges(experiment))[starts]
  regionEnd <- end(rowRanges(experiment))[ends]
  strands <- strand(rowRanges(experiment))[starts]


  df <- data.frame(chr=chromosomes, start=regionStart, end=regionEnd, strand=strands)

  #Check if any regions cross a chromosome. If so, remove them and print a message.
  errors <- which(as.vector(chromosomes != chromosomesCheck))
  if(length(errors) > 0){
    message(paste("Candidate Region(s)", errors, "cross chromosomes and are being removed."))
    df <- df[-errors,]
  }

  candidateRegions <- makeGRangesFromDataFrame(df)  # strand value "." is replaced with "*"
  candidateRegionsBins <- cbind(starts, ends) #Bin/window indices associated with each region.

  #Remove candidate regions that are not at least minWidth wide
  keep <- which(width(candidateRegions) >= minWidth)
  candidateRegions <- candidateRegions[keep]
  candidateRegionsBins <- candidateRegionsBins[keep,]

  #Remove regions that overlap with blacklist
  blacklist <- metadata(experiment)$blacklist

  if(!is.null(blacklist)){
    dontKeep <- subjectHits(findOverlaps(blacklist, candidateRegions))

    if(length(dontKeep) > 0){
      candidateRegions <- candidateRegions[-dontKeep]
      candidateRegionsBins <- candidateRegionsBins[-dontKeep,]

      if(verbose){
        message("Removing ", length(dontKeep), " candidate regions due to overlap with the blacklist.")
      }
    }
  }

  colnames(candidateRegionsBins) <- c("startBin", "endBin")


  #Record an index column
  elementMetadata(candidateRegions)[["index"]] <- 1:length(candidateRegions) #Preserve index number of region.

  #Return a dbListResults object with the candidate regions included.
  #if(class(experiment) == "dbList"){
  #  return(dbListResults(experiment, dbRanges = candidateRegions, dbRangesBins = candidateRegionsBins,
  #                       lower = lower, upper = upper, thresholdType = thresholdType, minWidth = minWidth,
  #                       maxGap = maxGap))
  #}

  #If not class dbList, then it's already of class dbListResults.
  #experiment@dbRanges <- candidateRegions
  #experiment@dbRangesBins <- candidateRegionsBins
  #experiment@lower <- lower
  #experiment@upper <- upper
  #experiment@thresholdType <- thresholdType
  #experiment@minWidth <- minWidth
  #experiment@maxGap <- maxGap

  #Add candidate regions and other info to experiment metadata
  metadata(experiment)$regions <- candidateRegions
  metadata(experiment)$regionsBins <- candidateRegionsBins
  metadata(experiment)$lower <- lower
  metadata(experiment)$upper <- upper
  metadata(experiment)$thresholdType <- thresholdType
  metadata(experiment)$minWidth <- minWidth
  metadata(experiment)$maxGap <- maxGap

  return(experiment)

}
