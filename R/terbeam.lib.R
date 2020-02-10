################################################################
#  _____       _                   ____            _           #
# |_   _|     | |                 / __ \          (_)          #
#   | |  _ __ | |_ ___  __ _ _ __| |  | |_ __ ___  _  ___ ___  #
#   | | | '_ \| __/ _ \/ _` | '__| |  | | '_ ` _ \| |/ __/ __| #
#   | |_| | | | ||  __| (_| | |  | |__| | | | | | | | (__\__ \ #
# |_____|_| |_|\__\___|\__, |_|   \____/|_| |_| |_|_|\___|___/ #
#                       __/ |                                  #
#                      |___/                                   #
################################################################

################################################################
# @script: terbeam.R                                          
# @author: Lucas Robin
# @author: Jean-Daniel Zucker
# @date: August 2016                                                    
################################################################


generateAllSingleFeatureModel <- function(X, y, clf, ind.sub = NULL, eval.all = FALSE)
{
  
  # if we wish to focus on a subset
  if(!is.null(ind.sub))
  {
    # sanity check
    if(any(is.na(match(ind.sub, seq_len(nrow(X))))))
    {
      stop("generateAllSingleFeatureModel: index does not match with X")
    }
    
  }else
  {
    # otherwise take all
    ind.sub <- seq_len(nrow(X))
  }
  
  # if we want to generate all the models
  if(clf$params$testAllSigns)
  {
    listOfVecPos <- as.list(ind.sub)
    clfPos <- clf
    clf$coeffs <- abs(clf$coeffs)
    popPos <- listOfSparseVecToListOfModels(X, y, clfPos, v = listOfVecPos)
    
    listOfVecNeg <- as.list(ind.sub)
    clfNeg <- clf
    clf$coeffs <- -abs(clf$coeffs)
    popNeg <- listOfSparseVecToListOfModels(X, y, clgNeg, v = listOfVecNeg)
    
    pop <- c(popPos, popNeg)
    
  } else
  {
    listOfVec <- as.list(ind.sub)
    pop <- listOfSparseVecToListOfModels(X, y, clf, v = listOfVec)
  }
  
  # return the population
  return(pop)
}



#' Generate every possible combination of a list of features and evaluate them
#' 
#' @title generateAllCombinations
#' @importFrom  gtools combinations
generateAllCombinations <- function(X, y, clf, ind.features.to.keep, sparsity, allFeatures)
{
  # Generate a matrix containing every combination of size sparsity of the 
  # features contained in ind.features.to.keep
  pop_vec <- combinations(n = length(ind.features.to.keep), r = sparsity, v = ind.features.to.keep)
  
  if(clf$params$language == "ratio" & clf$params$objective == "auc")
  {
    ind <- (apply(pop_vec,1,min) < length(ind.features.to.keep)/2) & (apply(pop_vec, 1, max) > length(ind.features.to.keep)/2)
    pop_vec <- pop_vec[ind,]
  }
  
  # if no more models are found
  if(class(pop_vec) != "matrix")
  {
    return(NULL)
  }
  
  # transform pop_vec from matrix to list of sparseVec objects
  pop_vec <- lapply(seq_len(nrow(pop_vec)), function(i, pop_vec) {pop_vec[i,]}, pop_vec)
  # transform the sparseVec objects onto predomics objects
  pop <- listOfSparseVecToListOfModels(X, y, clf, pop_vec)
  return(pop)
}


countEachFeatureApperance <- function(clf, allFeatures, pop, best, veryBest)
{
  if(!isPopulation(pop))
  {
    return(NULL)
  }
  
  featuresInPop      <- unique(populationGet_X(element2get = "names_", toVec = TRUE, na.rm = TRUE)(pop = pop))
  #featuresInPop      <- featuresInPop[order(featuresInPop)]
  featuresInBest     <- unique(populationGet_X(element2get = "names_", toVec = TRUE, na.rm = TRUE)(pop = best))
  featuresInVeryBest <- unique(populationGet_X(element2get = "names_", toVec = TRUE, na.rm = TRUE)(pop = veryBest))
  nbInBest           <- sapply(featuresInPop, 
                               function(x, featuresInBest) 
                               {sum(x == featuresInBest)},
                               featuresInBest)
  nbInVeryBest       <- sapply(featuresInPop, 
                               function(x, featuresInBest)
                               {sum(x == featuresInVeryBest)}, featuresInBest)
  # normalize to 1
  freqInBest <- nbInBest / max(nbInBest, na.rm = TRUE)
  
  feat.coeff <- clf$coeffs_[featuresInPop]
  
  res <- data.frame(nbInBest, nbInVeryBest, freqInBest, feat.coeff)
  res <- res[order(res[, "nbInVeryBest"], res[, "nbInBest"], decreasing = TRUE),]
  
  return(res)
}


getFeatures2Keep <- function(clf, featuresApperance, threshold = 0.01)
{
  # added threshold as the minimum percentage of models where a feature apprers in the best models
  if(is.null(featuresApperance))
  {
    return(NULL)
  }
  
  if(clf$params$language == "ratio" & clf$params$objective == "auc")
  {
    fa.pos <- featuresApperance[featuresApperance$feat.coeff == "1",]
    fa.neg <- featuresApperance[featuresApperance$feat.coeff == "-1",]
    
    tookeep.pos <- fa.pos$nbInVeryBest !=0 | fa.pos$freqInBest > threshold #| 1:clf$params$nbBest/2
    mask <- rep(FALSE,nrow(fa.pos)); mask[1:clf$params$nbBest/2] <- TRUE
    tookeep.pos <- tookeep.pos | mask
    
    tookeep.neg <- fa.neg$nbInVeryBest !=0 | fa.neg$freqInBest > threshold
    mask <- rep(FALSE,nrow(fa.neg)); mask[1:clf$params$nbBest/2] <- TRUE
    tookeep.neg <- tookeep.neg | mask
    
    tokeep <- c(rownames(fa.pos)[tookeep.pos], rownames(fa.neg)[tookeep.neg])
  }else
  {
    tokeep <- (featuresApperance[,"nbInVeryBest"] !=0  | featuresApperance[,"freqInBest"] > threshold)  
    tokeep <- rownames(featuresApperance[tokeep,])
  }
  
  return(tokeep)
}

