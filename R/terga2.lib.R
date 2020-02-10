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
# @script: terga2.lib.R                                          
# @author: Edi Prifti
# @author: Lucas Robin
# @date: August 2016                                                    
################################################################


######################################################
###               CLASS POPULATION               #####
######################################################
#' @import foreach
#' @import doRNG
population2 <- function(X, y, clf, featEval = NULL)
{
  # initialize featuEval
  if(is.null(featEval))
  {
    featEval <-rep(NA, length(rownames(X)))
    names(featEval) <- rownames(X)
  }
  
  pop <- list()
  if(!clf$params$randomSigns) # if we need to specify signs
  {
    if(clf$params$language=="bin" | clf$params$language=="bininter") 
    {
      signs  <- rep(1, nrow(X))
    }else # for the other languages get the signs
    {
      if(!is.null(clf$coeffs_))
      {
        signs <- clf$coeffs_
      }else # we compute here otherwise
      {
        signs  <- getSign(X = X, y = y, clf = clf)  # +1, -1 for the ternary learning    
      }
    }
  } else 
  {
    signs = NULL
  }
  
  if(clf$params$size_pop_random!=0){
    # create individuals using the signs
    if(clf$params$parallel.local)
    {
      #print("OUI!!!!!!!!!!!!! paralle!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      pop <- foreach(i = 1:clf$params$size_pop_random) %dorng%
      {
        individual(X, y, clf, signs = signs, eval.all = FALSE)
      }
    } else
    {
      for(i in 1:clf$params$size_pop_random)
      {
        pop[[i]] <- individual(X, y, clf, signs = signs, eval.all = FALSE)
      }
    }
  }

  
  # If there are individuals to inject in the population 
  # passed by the clf
  if(!(length(clf$params$in_pop)==1 & all(clf$params$in_pop=="NULL")))
  {
    pop_add <- clf$params$in_pop
    pop[(length(pop) +1):(length(pop) + length(pop_add))] <- pop_add
    
    if(clf$params$verbose)
    {
      cat("\t\tseeding population with ",length(pop_add)," injected individuals","\n")
    }
  }
  
  # or from a file
  if(clf$params$popSourceFile!="NULL")
  {
    if(length(clf$params$popSourceFile)==1)
    {
      popFromFile <- loadPopulation(clf$params$popSourceFile[[1]])
      pop[(length(pop) +1):(length(pop) + length(popFromFile))] <- popFromFile
      
      if(clf$params$verbose)
      {
        cat("\t\tseeding population with ",length(popFromFile)," individuals loaded from file","\n")
      }
    }else
    {
      # if more than one file
      for(i in 1:length(clf$params$popSourceFile))
      {
        popFromFile <- loadPopulation(clf$params$popSourceFile[[i]])
        pop[(length(pop) +1):(length(pop) + length(popFromFile))] <- popFromFile
        
        if(clf$params$verbose)
        {
          cat("\t\tseeding population with ",length(popFromFile)," individuals loaded from file ",i,"\n")
        }
      } # end for
    } # end more files
  } # end load from files
  
  # generate randomly created individuals
  if(clf$params$size_pop_random!=0)
  {
    # make other individuals
    if(!clf$params$evolve_k1)
    {
      # generate all the single feature models
      singleFeatMod <- generateAllSingleFeatureModel(X, y, clf, eval.all = FALSE)
      # add them to the main initial population
      pop[(length(pop) +1):(length(pop) + length(singleFeatMod))] <- singleFeatMod
      
      if(clf$params$verbose)
      {
        cat("\t\tseeding population with ",length(singleFeatMod)," randomly created individuals ","\n")
      }
    }
  } # end random generation
 
  if(!isPopulation(pop))
  {
    return(NULL)
  }
  
  # evaluate only the population but only the evalToFit and what it needs to be computed
  pop <- evaluatePopulation(X, y, clf, pop, eval.all = FALSE)
  pop <- unique(pop) # keep only unique samples
  
  # sort the population by fitness
  pop <- sortPopulation(pop, evalToOrder = "fit_")

  if(clf$params$verbose)
  {
    cat("\t\tThe final population is composed of ",length(pop)," individuals","\n")
  }
  return(pop)
}


# One working version of the individual function which is used in population2
individual_vec_v1 <- function(clf, signs = NULL) # changer à dense
{
  prob_1 <- clf$params$sparsity.mean/clf$params$size_world
  prob <- c(prob_1/2, 1 - prob_1, prob_1/2)
  
  # Check binary language
  if(clf$params$language=="bin" | clf$params$language=="bininter") # for binary language
  {
    res <- sample(1, clf$params$size_world, replace = TRUE, prob = prob)
  }else # for ternary language
  {
    res <- sample(-1:1, clf$params$size_world, replace = TRUE, prob = prob)  
  }
  
  effectiveSpar <- length(which(res != 0))
  if(effectiveSpar < clf$params$sparsity.min)
  {
    if(clf$params$sparsity.min < clf$params$sparsity.max)
    {
      ind2Add.nb <- sample((clf$params$sparsity.min - effectiveSpar):(clf$params$sparsity.max - effectiveSpar), 1)
    } else
    {
      ind2Add.nb <- clf$params$sparsity.min - effectiveSpar
    }
    notNull <- which(res != 0)
    world <- (1:clf$params$size_world)
    if(length(notNull) > 0)
      world <- world[-notNull]
    
    ind2Add <- sample(world, ind2Add.nb)
    
    # Check binary language
    if(clf$params$language=="bin" | clf$params$language=="bininter") # for binary language
    {
      res[ind2Add] <- sample(c(1), length(ind2Add), replace = TRUE)
    }else # for ternary language
    {
      res[ind2Add] <- sample(c(-1,1), length(ind2Add), replace = TRUE)
    }

  } else if(effectiveSpar > clf$params$sparsity.max) 
  {
    if(clf$params$sparsity.min < clf$params$sparsity.max)
    {
      ind2rm.nb <- sample((effectiveSpar - clf$params$sparsity.max):(effectiveSpar - clf$params$sparsity.min), 1)
    } else 
    {
      ind2rm.nb <- effectiveSpar - clf$params$sparsity.max
    }
    notNull <- which(res != 0)
    ind2Rm <- sample(notNull, ind2rm.nb)
    res[ind2Rm] <- 0
  }
  return(res)
}


# Another one, they are here mainly to be examples
individual_vec_v2 <- function(clf, signs = NULL)
{
  if(clf$params$size_world == "NULL")
  {
    if(clf$params$warnings) warning("individual_vec_v2: the size_world is not set while expected")
    return(NULL)
  }
  
  if(clf$params$sparsity.min != clf$params$sparsity.max)
  {
    notNull.number    <- sample(clf$params$sparsity.min:clf$params$sparsity.max, 1)
  } else
  {
    notNull.number    <- clf$params$sparsity.mean
  }
  notNull.ind         <- sample(1:clf$params$size_world, notNull.number)
  res                 <- rep(0, clf$params$size_world)
  
  if(is.null(signs)) 
  {
    # Check binary language
    if(clf$params$language=="bin" | clf$params$language=="bininter") # for binary language
    {
      res[notNull.ind]  <- sample(c(1), notNull.number, replace = TRUE)
    }else # for ternary language
    {
      res[notNull.ind]  <- sample(c(-1, 1), notNull.number, replace = TRUE)
    }
    
  } else 
  {
    res[notNull.ind]  <- signs[notNull.ind]
  }
  return(res)
}


####### Watch out, this function returns a sparse vector. NOT WORKING FOR THE MOMENT
individual_vec_v3 <- function(clf, signs = NULL)
{
  notNull.number   <- sample(clf$params$sparsity.min:clf$params$sparsity.max, 1)
  notNull.ind      <- sample(1:clf$params$size_world, notNull.number)
  
  # Check binary language
  if(clf$params$language=="bin" | clf$params$language=="bininter") # for binary language
  {
    res <- sample(c(1), notNull.number, replace = TRUE)
  }else # for ternary language
  {
    res <- sample(c(-1, 1), notNull.number, replace = TRUE)
  }
  return(list(indices = notNull.ind, coeffs = res))
}


####### Tests effectué sur les deux fonctions précédentes : #################################################################
# microbenchmark::microbenchmark(                                                                                           #
#   test <- foreach(k = 1:1000) %do% {individual_vec_v1(clf.terga)},                                                        #
#   test <- foreach(k = 1:1000) %do% {individual_vec_v2(clf.terga)},                                                        #
#   test <- foreach(k = 1:1000) %do% {individual_vec_v3(clf.terga)},                                                        #
#   times = 100                                                                                                             #
# )                                                                                                                         #
# Unit: milliseconds                                                                                                        #
#                                                            expr      min       lq     mean   median       uq max    neval #
# test <- foreach(k = 1:1000) %do% {individual_vec_v1(clf.terga)} 243.6722 249.3285 254.2703 252.2820 255.5914 325.4101 100 #
# test <- foreach(k = 1:1000) %do% {individual_vec_v2(clf.terga)} 246.8161 253.1758 260.6690 255.8954 260.4250 370.5941 100 #
# test <- foreach(k = 1:1000) %do% {individual_vec_v3(clf.terga)} 244.4738 250.4082 257.0116 252.7833 256.8256 357.6657 100 #
#############################################################################################################################


#' Tag individuals for parenting
#' @description Function to add the tag "selected" to the best individuals of the population
#' @param clf: the classifier object
#' @param pop: the population on which we want to add the tag
#' @param nbToSelect: the number of individuals we are going to select in the population
#' @return the population given as an input with `nbToSelect` bests individuals with `$selected = TRUE`
tag_SelectElite <- function(clf, pop, nbToSelect)
{
  pop <- sortPopulation(pop, evalToOrder = "fit_")
  
  if(is.null(pop))
  {
    warning("... ... resetTags: the population is empty\n")
    return(NULL)
  }
  
  nbSelected <- 0
  for(i in 1:length(pop))
  {
    if(nbSelected < nbToSelect)
    {
      if(!pop[[i]]$selected)
      {
        pop[[i]]$selected <- TRUE
        nbSelected <- nbSelected +1
      }
    }
  }
  return(pop)
}


#' Randomly tag selected individuals parenting
#' @description This function turns the selected switch on when an individual is 
#' selected to survive the generation and be among the pool of parents for the 
#' next generation.
#' @param clf: The classifier object
#' @param pop: The population on which the selection process will be performed.
#' @param nbToSelect: the number of individuals we are going to select in the population
#' @return the population given as an input with `nbToSelect` individuals with `selected = TRUE`
tag_SelectRandom <- function(clf, pop, nbToSelect)
{
  if(is.null(pop))
  {
    if(clf$params$warnings) warning("resetTags: the population is empty")
    return(NULL)
  }
  
  indexNonSelected <- list()
  for(i in 1:length(pop))
  {
    if(!pop[[i]]$selected)
    {
      indexNonSelected[[length(indexNonSelected) +1]] <- i
    }
  }
  
  if(length(indexNonSelected) < nbToSelect)
  {
    selection <- sample(indexNonSelected, length(indexNonSelected))  
  }else
  {
    selection <- sample(indexNonSelected, nbToSelect)
  }
  
  if(clf$params$parallel.local)
  {
    pop <- foreach(i = 1:length(pop)) %dorng%
    {
      if(i %in% selection)
      {
        pop[[i]]$selected <- TRUE
      }
      return(pop[[i]])
    }
  } else
  {
    for(i in 1:length(pop))
    {
      if(i %in% selection)
      {
        pop[[i]]$selected <- TRUE
      }
    }
  }
  return(pop)
}

#' Tag the couples
#' @description This function selects constitutes the couples that will give the
#' next generation individuals by adding the couple id on the mate attribute.
#' @param pop: The population on which we the couples are being constituted;
#' @param parents: The parent candidate individuals from which the couples will 
#' be selected.
#' @return The parent population with the couple tags set.
tag_Couples <- function(pop, parents)
{
  if(is.null(pop))
  {
    if(clf$params$warnings) warning("resetTags: the population is empty")
    return(NULL)
  }
  
  for(i in 1:(length(parents)/2))
  {
    if(length(parents) > 1)
    {
      couple <- sample(parents, 2)
      couple <- couple[order(couple)]
      parents <- parents[-which(parents == couple[[1]])]
      parents <- parents[-which(parents == couple[[2]])]
      
      pop[[couple[1]]]$mate <- couple[2]
      pop[[couple[2]]]$mate <- -couple[1]
    }
  }
  return(pop)
}


#' Tag individuals for mutation
#' @description Function to add the tag "toBeMutated" to a randomly sampled part of 
#' the population. Some individuals (the best ones) will be protected from the mutation
#' so that genetic decline does not happen.
#' @param pop: The population on which the individuals to be mutated will be selected
#' @param mutate_size: The number of individuals to mutate
#' @param protected: The index of individuals which should not be mutated.
#' @return The population given as an input with `mutate_size` individuals with `toBeMutated = TRUE`
tag_ToBeMutated <- function(pop, mutate_size, protected = NULL)
{
  if(is.null(pop))
  {
    if(clf$params$warnings) warning("resetTags: the population is empty")
    return(NULL)
  }
  
  all.indexes <- seq_len(length(pop))
  
  if(!is.null(protected))
  {
    # sanity check
    if(!all(protected %in% all.indexes))
    {
      stop("tag_ToBeMutated: protected indices don't match the population")
    }
    possible.to.mutate <- (length(pop)-length(protected))
    
# TODO fix negative delta ?????
    delta = round(mutate_size) - possible.to.mutate
    
    if(round(mutate_size) > possible.to.mutate)
    {
      if(clf$params$warnings) warning(paste("tag_ToBeMutated: the protected individuals don't leave enough to mutate: nb protected",
                    length(protected), " and nb available", possible.to.mutate,"mutating less ..."))
      protected <- protected[1:(length(protected) - delta)] # stop protecting some to assure the mutation quota.
    }
    tosample.indexes <- which(!all.indexes %in% protected)
  }else # if none is protected
  {
    tosample.indexes <- all.indexes
  }

  # select mutate_size indivuduals out of the not protected
  indexToBeMutated <- sample(tosample.indexes, mutate_size)
  # set the indexes to be mutated
  for(i in indexToBeMutated)
  {
    pop[[i]]$toBeMutated <- TRUE
  }
  
  return(pop)
}



#' Add `selected` tag using elite and random selection
#' @description This function combines \link[predomics]{tag_SelectElite} and \link[predomics]{tag_SelectRandom}
#' to tag the desired individuals in a population following the proportion given in the clf
#' @param X: Unused but still here for compatibility
#' @param y: same as `X`
#' @param clf: the classifier object where parameters are defined
#' @param pop: the population on which we want to apply the selection
#' @return the population with the tag `selected` on some of the individuals
tag_select <- function(X, y, clf, pop)
{
  pop <- sortPopulation(pop, evalToOrder = "fit_")
  
  nb2BeMutated <- ceiling(clf$params$select_perc * clf$params$size_pop/100)
  nbByMethod   <- sapply(clf$params$select_percByMethod,        # On supose qu'il n'y a que 2 methodes de selection,
                         function(x) {x * nb2BeMutated/100})    # et qu'elles sont : elite et random
  selectElite   <- 1:nbByMethod[1]
  
  selectRandom <- sample(length(pop), nbByMethod[2])
  for(x in selectElite)
  {
    # At the current version an individual can be crossed with itself
    pop[[x]]$mate[length(pop[[x]]$mate)+1] <- sample(unique(selectElite, selectRandom), 1) 
  }
  for(x in selectRandom)
  {
    # At the current version an individual can be crossed with itself
    pop[[x]]$mate[length(pop[[x]]$mate)+1] <- sample(unique(selectElite, selectRandom), 1)
  }
  
  # tag to be mutated
  pop <- tag_ToBeMutated(pop = pop, mutate_size = (clf$params$size_pop * clf$params$mutate_size / 100), protected = NULL)
  
  return(pop)
}



#' Return list of parents
#' @description Wraper function that returns the list of parents (normaly taged with 
#' \link[predomics]{tag_SelectRandom})
#' @param pop: population list
#' @return a list of individuals (= a population) containing all the selected parents
get_Parents <- function(pop)
{
  return(which(populationGet_X(element2get = "selected", toVec = TRUE, na.rm = TRUE)(pop)))
}


#' Return list of individuals to mutate
#' @description Wraper function that returns the list of individuals that are going to go through 
#' the mutate function (normaly taged with \link[predomics]{tag_SelectRandom})
#' @param pop: population list
#' @return a list of individuals (= a population) containing all the individuals selected for mutation
get_IndividualToBeMutated <- function(pop)
{
  return(which(populationGet_X(element2get = "toBeMutated", toVec = TRUE, na.rm = TRUE)(pop)))
}


#' Resets selection, mutation and mate tags to inactive
#' @description Resets selection, mutation and mate tags to inactive
#' @param pop: The population to be evolved
#' @param selected: set to (default:FALSE) the selected attribute, if not null.
#' @param toBeMutated: set to (default:FALSE) the selected attribute, if not null.
#' @param mate: set to (default:-1) the selected attribute, if not null.
#' @return A modified population
resetTags <- function(pop, selected = FALSE, toBeMutated = FALSE, mate = -1)
{
  
  if(is.null(pop))
  {
    if(clf$params$warnings) warning("resetTags: the population is empty")
    return(NULL)
  }
  
  if(!is.list(pop))
  {
    stop("resetTags: the population is not a list")
  }
  
  for(i in 1:length(pop))
  {
    if(!is.null(selected)) pop[[i]]$selected <- selected
    if(!is.null(toBeMutated)) pop[[i]]$toBeMutated <- toBeMutated
    if(!is.null(mate)) pop[[i]]$mate <- mate
  }
  return(pop)
}


evolve1 <- function(X, y, clf, pop, featEval) 
{
  if(anyNA(lapply(pop, function(x) x$fit_)))
    pop <- evaluatePopulation(X, y, clf, pop, eval.all = FALSE)
  
  # 1. Evaluate the current population
  # evaluation <- getFitPopulation(pop)
  evaluation <- populationGet_X(element2get = "fit_", toVec = TRUE, na.rm = TRUE)(pop)
  best_individual <- pop[[which.max(abs(evaluation))]]
  evolved_pop <- pop
  
  #   if (clf$params$verbose) cat("\tBest individual: ", best_individual$indices_,
  #                               "\t\tAUC:\t", best_individual$fit_, "\n")
  
  # 2 Select the parents of the next generation (indexes in the population)
  if(clf$params$cross)
  {
    parents <- select(clf, pop)
    evolved_pop <- crossing2(X, y, clf, evolved_pop, parents)
    evolved_pop <- evaluatePopulation(X, y, clf, evolved_pop, eval.all = FALSE)
    evaluation <- populationGet_X(element2get = "fit_", toVec = TRUE, na.rm = TRUE)(evolved_pop)
  }
  
  # 4. Mutate the new population
  if(clf$params$mutate)
  {
    selection     <- !(1:length(evaluation)) %in% which.max(abs(evaluation))
    selection.ind <- sample(which(selection), ceiling(sum(selection)*clf$params$mutate_size/100))
    evolved_pop   <- mutate2(X, y, clf, evolved_pop, selection.ind, featEval)
    evolved_pop   <- evaluatePopulation(X, y, clf, evolved_pop, force.re.evaluation = TRUE, eval.all = FALSE)
    evaluation <- populationGet_X(element2get = "fit_", toVec = TRUE, na.rm = TRUE)(evolved_pop)
  }
  
  best_individual <- evolved_pop[[which.max(abs(evaluation))]]
  if (clf$params$verbose) 
  {
    if(isModel(best_individual))
    {
      try(print(printModel(mod = best_individual, method = clf$params$print_ind_method, score = "fit_")), silent = TRUE)
    }
  }
  
  return(evolved_pop)
}



#' @title Second version of the evolve method
#' @description This evolve method realize the selection of the parents with 
#' two methods (for the moment) : elite and random. The it tags every selected 
#' individual with the index of it's mate in the population. Then the individual 
#' which sould be mutated are tagged. for each individual of the population we 
#' check if they need to be crossed and/or mutated and the we apply the operations 
#' and create a new individual for each operation applied.
#' @param X: Dataset to classify
#' @param y: Variable to predict
#' @param clf: The claffifier  object containing the different parameters
#' @param pop: The population to be evolved
#' @param featEval: A dataframe with the evaluation of each variable of the dataset, 
#' used for some mutator.
#' @return A list with the size_pop bests of the combination of the old population 
#' with the new population
evolve2m <- function(X, y, clf, pop, featEval, generation)
{
  pop <- sortPopulation(pop, evalToOrder = "fit_")
  if(!isPopulation(obj = pop)) return(NULL)
  
  if(clf$params$debug) {print(paste("evolve2m after sort:", length(pop)))}
  
  nb2BeMutated <- ceiling(clf$params$select_perc * clf$params$size_pop/100)
  nbByMethod   <- lapply(clf$params$select_percByMethod,        # On supose qu'il n'y a que 2 methodes de selection,
                         function(x) {x * nb2BeMutated/100})    # et qu'elles sont : elite et random
  
  pop <- resetTags(pop)
  
  # for(i in 1:length(nbByMethod))
  # {
  pop <- tag_SelectElite(clf, pop, nbByMethod[[1]])
  pop <- tag_SelectRandom(clf, pop, nbByMethod[[2]])
  # }
  indexListOfAvailableParents <- get_Parents(pop)
  pop <- tag_Couples(pop, indexListOfAvailableParents)
  
  mutate_size <- (length(pop) * clf$params$mutate_size / 100)
  pop <- tag_ToBeMutated(pop = pop, mutate_size = mutate_size, protected = NULL)
  
  # indexListOfIndividualToBeMutated <- get_IndividualToBeMutated(pop)
  
  if(clf$params$debug) {print(paste("evolve2m after tagging:", length(pop)))}
  
  if(clf$params$parallel.local)
  {
    newPop <- foreach(i = 1:length(pop)) %dorng%
    {
      child <- list()
      mutatedIndiv <- list()
      
      if(myAssertNotNullNorNa(pop[[i]]$selected) & myAssertNotNullNorNa(pop[[i]]$mate) & myAssertNotNullNorNa(pop[[i]]$toBeMutated))
      {
        # Appliquer les operations sur pop[[i]]
        if(pop[[i]]$selected && (pop[[i]]$mate > 0))
        {
          child_1 <- clf$functions$crosser(X, y, clf, pop[[i]], pop[[pop[[i]]$mate]])
          
          # if parents have the same language
          if(pop[[i]]$language == pop[[pop[[i]]$mate]]$language)
          {
            # make only one child
            child[[length(child) +1]]<-evaluateModel(mod = child_1, X = X, y = y, clf = clf, eval.all = FALSE, force.re.evaluation = TRUE)
          }else
          {
            # otherwise make two identical with each of the parent's languages
            child_2 <-child_1
            child_1$language <- pop[[i]]$language
            child_2$language <- pop[[pop[[i]]$mate]]$language
            child[[length(child) +1]]<-evaluateModel(mod = child1, X = X, y = y, clf = clf, eval.all = FALSE, force.re.evaluation = TRUE)
            child[[length(child) +1]]<-evaluateModel(mod = child2, X = X, y = y, clf = clf, eval.all = FALSE, force.re.evaluation = TRUE)
          }
        }
        
        if(pop[[i]]$toBeMutated)
        {
          mutatedIndiv <- clf$functions$mutator(X, y, clf, pop, pop[[i]], 1:clf$params$size_world, featEval)
          mutatedIndiv <- evaluateModel(mod = mutatedIndiv, X = X, y = y, clf = clf, eval.all = FALSE, force.re.evaluation = TRUE)
        }
      } # end iff test
      
      return(list(child = child, mutatedIndiv = mutatedIndiv))
    }
    
    mutatedPop <- lapply(newPop, function(x)
    {
      if(length(x$mutatedIndiv) > 0)
      {
        return(x$mutatedIndiv)
      } else
      {
        return(NULL)
      }
    })
    mutatedPop <- mutatedPop[!sapply(mutatedPop, is.null)]
    
    children <- list()
    for(l in 1:length(newPop))
    {
      if(length(newPop[[l]]$child) > 0)
      {
        for(m in 1:length(newPop[[l]]$child))
        {
          children[[length(children)+1]] <- newPop[[l]]$child[[m]]
        }
      }
    }
    
    newPop <- c(children, mutatedPop)
    
  }else # if not parallel
  {
    newPop <- list()
    # for each of the size_pop best individuals (smaller then the real population size)
    for(i in 1:length(pop))
    {
      if(myAssertNotNullNorNa(pop[[i]]$selected) & myAssertNotNullNorNa(pop[[i]]$mate) & myAssertNotNullNorNa(pop[[i]]$toBeMutated))
      {
        if(pop[[i]]$selected && (pop[[i]]$mate > 0))
        {
          # create the first child
          child_1 <- clf$functions$crosser(X, y, clf, pop[[i]], pop[[pop[[i]]$mate]])
          
          # if parents languages are the same we keep only one child
          if(pop[[i]]$language == pop[[pop[[i]]$mate]]$language)
          {
            child <- evaluateModel(mod = child_1, X = X, y = y, clf = clf, eval.all = FALSE, force.re.evaluation = TRUE)
            newPop[[length(newPop) +1]]  <- child
            
          }else # if they have different languages (we create two identical children, each with one language)
          {
            child_2 <- child_1
            child_1$language <- pop[[i]]$language
            child_2$language <- pop[[pop[[i]]$mate]]$language
            newPop[[length(newPop) +1]] <- evaluateModel(mod = child_1, X = X, y = y, clf = clf, eval.all = FALSE, force.re.evaluation = TRUE)  
            newPop[[length(newPop) +1]] <- evaluateModel(mod = child_2, X = X, y = y, clf = clf, eval.all = FALSE, force.re.evaluation = TRUE)
          }
        }
        
        if(pop[[i]]$toBeMutated)
        {
          mutatedIndiv <- clf$functions$mutator(X, y, clf, pop, individual_to_be_mutated = pop[[i]], all_genes = c(1:clf$params$size_world), featEval)
          newPop[[length(newPop) +1]] <- evaluateModel(mod = mutatedIndiv, X = X, y = y, clf = clf, eval.all = FALSE, force.re.evaluation = TRUE)  
        }
      }
    } # end for each model in the population size
  } # end if no parallel
  
  if(clf$params$debug) {print(paste("evolve2m after crossing/mutation:", length(pop)))}

  pop <- sortPopulation(pop, evalToOrder = "fit_")
  pop <- pop[1:min(length(pop),round(clf$params$size_pop*2/3))] # get the 2/3 population or the size of the pop if smaller
  
  if(clf$params$debug) {print(paste("evolve2m after selection:", length(pop)))}
  
  # set population size bigger at the first generation
  if(generation==0) { clf$params$size_pop <- round(clf$params$size_pop*3/2) }
  
  n_rest <- clf$params$size_pop - length(pop)
  # add in the rest of the population 1/3, the best models
  pop <- c(pop, sortPopulation(newPop, evalToOrder = "fit_")[1:n_rest])
  #pop[(length(pop)+1):(length(pop)+ n_rest)] <- sortPopulation(newPop, evalToOrder = "fit_")[1:n_rest]
  # new2keep <-  selectIndividualToKeep(clf, newPop)
  # pop[(length(pop)+1):(length(pop)+ length(new2keep))] <- new2keep
  
  # add the whole new population
  #pop[(length(pop)+1):(length(pop) + length(newPop))] <- newPop
  pop <- c(pop, newPop)
  pop <- cleanPopulation(pop = pop, clf = clf)
  pop <- unique(pop) # clean
  pop <- sortPopulation(pop, evalToOrder = "fit_") # sort
  pop <- resetTags(pop) # reset
  
  if(clf$params$debug) {print(paste("evolve2m after merging the new pop:", length(pop)))}
  
  return(pop)
}



#' @title Another version of the evolve method that distingishes the different lanugages
#' @description  In a nutshell this evolve method will select the parents using 2 methods (elite and random).
#' Then it tags every selected individual with the index of it's mate in the population. 
#' Then the individuals which are randomly selected and tagged as to be mutated. For each individual of the 
#' population we check if they need to be crossed and/or mutated and the we apply the operations and 
#' create a new individual for each operation applied. If two parents are of different languages the algorithm 
#' will produce two children that are the same with each having one of the languages.
#' This as been as a pipeline that could be paralellized with the most efficiency. We could summarize it the following way : 
#' IN -> Evaluation -> initialisation of the tags -> Tagging -> Crossing and mutation -> OUT (newpop)
#' @param X: dataset to classify
#' @param y: variable to predict
#' @param clf: an object containing the different parameters of the classifier
#' @param pop: the population to be evolved
#' @param featEval: a dataset with the evaluation of each variable of the dataset, used for some mutator.
#' @return A list with every new individuals
evolve3m <- function(X, y, clf, pop, featEval) 
{
  
  if(clf$params$debug) # DEBUG
  {
    pop.lengths <- list()
    pop.lengths["pop"] <- length(pop)
  }
  # evaluate the population
  pop <- evaluatePopulation(X, y, clf, pop, eval.all = FALSE)
  
  # set two empty parameters to each individual (mate and toBeMutated)
  for(i in 1:length(pop))
  {
    pop[[i]]$mate <- c()
    pop[[i]]$toBeMutated <- FALSE
  }
  
  # Tag the population
  pop <- tag_select(X, y, clf, pop)
  
  ###################################################################################################################                         
  # function(x) {x * nb2BeMutated/100})    # et qu'elles sont : elite et random
  
  if(clf$params$parallel.local)
  {
    newPop <- foreach(i = 1:length(pop)) %dorng%
    {
      children <- list()
      mutatedIndiv <- list()
      # Appliquer les operations sur pop[[i]]
      if(length(pop[[i]]$mate) > 0)
      {
        for(j in 1:length(pop[[i]]$mate))
        {

          
          child<-clf$functions$crosser(X, y, clf, pop[[i]],  pop[[pop[[i]]$mate[j]]])
          if(pop[[i]]$language==pop[[pop[[i]]$mate[j]]]$language){
            children[[length(children) +1]]<-evaluateModel(mod = child, X = X, y = y, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE)
          }else
          {
            child_2 <-child
            child$language <- pop[[i]]$language
            child_2$language <- pop[[pop[[i]]$mate[j]]]$language
            children[[length(children) +1]]<-evaluateModel(mod = child, X = X, y = y, clf = clf, child, eval.all = TRUE, force.re.evaluation = TRUE)
            children[[length(children) +1]]<-evaluateModel(mod = child_2, X = X, y = y, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE)
          }
          
        }
      }
      
      if(pop[[i]]$toBeMutated)
      {
        mutatedIndiv <- clf$functions$mutator(X, y, clf, pop, pop[[i]], 1:clf$params$size_world, featEval)
        mutatedIndiv <- evaluateModel(mod = mutatedIndiv, X = X, y = y, clf = clf, eval.all = TRUE)
      }
      
      return(list(children = children, mutatedIndiv = mutatedIndiv))
    }
    
    mutatedPop <- lapply(newPop, function(x) 
    {
      if(length(x$mutatedIndiv) > 0)
      {
        return(x$mutatedIndiv)
      } else
      {
        return(NULL)
      }
    })
    mutatedPop <- mutatedPop[!sapply(mutatedPop, is.null)]
    
    if(clf$params$debug) # DEBUG
    {
      pop.lengths["mutated_pop"] <- length(mutatedPop)
    }
    
    children <- list()
    for(l in 1:length(newPop))
    {
      if(length(newPop[[l]]$children) > 0)
      {
        for(m in 1:length(newPop[[l]]$children))
        {
          children[[length(children)+1]] <- newPop[[l]]$children[[m]]
        }
      }
    }
    cat("Taille de la population : ", length(children) + length(mutatedPop), "\n")
    newPop <- c(children, mutatedPop)
    
    if(clf$params$debug) # DEBUG
    {
      pop.lengths["newPop"] <- length(newPop)
    }
    
  } else # else if not in parallel
  {
    newPop <- list()
    for(i in 1:length(pop))
    {
      if(length(pop[[i]]$mate) > 0) # if there at least one mate this is to become a parent
      {
        for(j in 1:length(pop[[i]]$mate)) # for all mates
        {
          if(pop[[i]]$language == pop[[pop[[i]]$mate[j]]]$language){
            child <- clf$functions$crosser(X=X, y=y, clf=clf, parent1=pop[[i]],  parent2=pop[[pop[[i]]$mate[j]]])
            newPop[[length(newPop) +1]] <- evaluateModel(mod = child, X = X, y = y, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE)
          }else # if parents are different (in terms of languages)
          { # compute two metisse childs
            parent1 <- pop[[i]]
            parent2 <- pop[[pop[[i]]$mate[j]]]
            # create a crossed individual
            child_1 <- clf$functions$crosser(X=X, y=y, clf=clf, parent1=parent1,  parent2=parent2)
            # duplicate it
            child_2 <- child_1
            child_1$languge <- parent1$language
            child_2$languge <- parent2$language
            newPop[[length(newPop) +1]] <- evaluateModel(mod = child_1, X = X, y = y, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE)
            newPop[[length(newPop) +1]] <- evaluateModel(mod = child_2, X = X, y = y, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE)
          }
          
        } # end mate loop
      }
      
      if(pop[[i]]$toBeMutated)
      {
        mutatedIndiv <- clf$functions$mutator(X, y, clf, pop, pop[[i]], 1:clf$params$size_world, featEval)
        newPop[[length(newPop) +1]] <- evaluateModel(mod = mutatedIndiv, X = X, y = y, clf = clf, eval.all = TRUE)
      } # end mutation
    }
    
    if(clf$params$debug) # DEBUG
    {
      pop.lengths["newPop"] <- length(newPop)
    }
  }
  
  
  #t<-length(newPop)
  #newPop[(length(newPop)+1):(length(newPop)+length(pop))]<-pop
  #newPop <- unique(newPop)
  #newPop <- sortPopulation(newPop, evalToOrder = "fit_")[1:t]
  
  if(clf$params$debug) # DEBUG
  {
    print(unlist(pop.lengths))
  }
  
  return(newPop)
}



mutator_v1 <- function(X, y, clf, pop, individual_to_be_mutated, all_genes, featEval)
{
  # degree of mutation (number of genes to be mutated)
  perc <- ceiling(individual_to_be_mutated$eval.sparsity * clf$params$mutate_rate/100)
  
  # the size of the reservoir
  size_reservoir <- clf$params$size_world - individual_to_be_mutated$eval.sparsity
  
  if(perc > size_reservoir)  # if the gene reservoir is smaller than the number of genes to mutate, 
  {                          # lower the mutation rate to the size of the reservoir
    perc <- size_reservoir
  }
  # identify the idexes of the genes to be mutated
  index_genes_to_mutate <- sample(individual_to_be_mutated$eval.sparsity, perc)
  # the unique genes remaining in the reservoir that are not in the untouched part of the genome
  # It is to be noted that the genes to mutate are part of this reservoir and some of them can be 
  # picked again, which lowers the percentage of mutation
  remaining_unique_genes <- all_genes[!all_genes %in% individual_to_be_mutated$indices_[-index_genes_to_mutate]]
  # draw new genes (some of them might be old)
  new_genes <- sample(x = remaining_unique_genes, size = perc)
  # the mutated individual
  mutated_individual <- individual_to_be_mutated$indices_
  mutated_individual[index_genes_to_mutate] <- new_genes
  mutated_individual <- sort(mutated_individual)

  return(individual(X, y, clf, ind = mutated_individual))
}



mutator_v2 <- function(X, y, clf, pop, individual_to_be_mutated, all_genes, featEval) #This one does not use the mutate rate
{
  ## start idea for feature mutation based on their evaluation 
  featEval[is.na(featEval)] <- rep(max(featEval, na.rm = TRUE), length(featEval[is.na(featEval)])) # for the NA fill them up
  featProba <- featEval/sum(featEval)
  
  if(any(featEval < 0))
  {
    stop("mutator_v2: Negative probabiliy in featProba")
  }
  
  # which mutation event to apply
  set.seed(clf$params$current_seed + sum(individual_to_be_mutated$indices_)) # fix the seed this individual
  
  if(clf$params$language == "bin" | clf$params$language == "bininter")
  {
    mutation2apply <- sample(c("add", "rm"), 1)
  }else
  {
    mutation2apply <- sample(c("modif", "add", "rm"), 1)
  }
  
  if(max(clf$params$sparsity) == min(clf$params$sparsity))
  {
    mutation2apply <- "modif"
  } else if((mutation2apply == "add") &&
            (length(individual_to_be_mutated$indices_) >= max(clf$params$sparsity)))
  {
    mutation2apply <- "rm"
  } else if((mutation2apply == "rm") && 
            (length(individual_to_be_mutated$indices_) <= min(clf$params$sparsity)))
  {
    mutation2apply <- "add"
  }
  switch(mutation2apply, 
         modif= # mutate
         {
           # select the feature to modify
           gene2Modify <- sample.int(length(individual_to_be_mutated$indices_), 1)
           # create a copy to manipulate
           mutated_individual <- individual_to_be_mutated
           # change the coefficient of that index
           mutated_individual$coeffs_[[gene2Modify]] <- individual_to_be_mutated$coeffs_[[gene2Modify]] * -1
           return(mutated_individual)
         },
         add= # insert
         {
           # which genes are outside the individuals from which we can draw
           remaining_unique_genes <- all_genes[!all_genes %in% individual_to_be_mutated$indices_]
           
           # If enough genes are left to draw from
           if(length(remaining_unique_genes) > 1) 
           {
             if(length(unique(remaining_unique_genes)) < length(remaining_unique_genes))
             {
               stop("mutator_v2: list of genes from which to choose contains identical values")
             }
             
             # sample a gene do to add
             set.seed(clf$params$current_seed + sum(individual_to_be_mutated$indices_)) # fix the seed this individual  
             gene2Add <- sample(all_genes, 1, prob = featProba)
             
             if(!gene2Add %in% individual_to_be_mutated$indices_)
             {
               indices <- individual_to_be_mutated$indices_
               indices[[length(indices) +1 ]] <- gene2Add
               indices <- indices[order(indices)]
               return(individual(X, y, clf, ind = indices))
             } else
             {
               return(individual_to_be_mutated)
             }
           } else
           {
             if(clf$params$warnings) warning("mutator_v2: Not enough features remaining for the insertion mutation. Nothing happened")
             return(individual_to_be_mutated)
           }
         },
         rm= # delete
         {
           # sample a gene do to remove
           set.seed(clf$params$current_seed + sum(individual_to_be_mutated$indices_)) # fix the seed this individual 
           gene2Rm <- sample.int(length(individual_to_be_mutated$indices_), 1)
           mutated_individual <- individual_to_be_mutated
           mutated_individual$indices_ <- individual_to_be_mutated$indices_[-gene2Rm]
           mutated_individual$coeffs_  <- individual_to_be_mutated$coeffs_[-gene2Rm]
           mutated_individual$names_   <- individual_to_be_mutated$names_[-gene2Rm]
           return(mutated_individual)
         }
  ) # end switch
}


rankFeatures <- function(X, y, clf, pop, featEval)
{
  featuresInModels <- populationGet_X(element2get = "indices_", toVec = FALSE, na.rm = FALSE)(pop)
  for(modIndice in 1:length(featuresInModels))
  {
    listFeat <- featuresInModels[[modIndice]]
    for(i in 1:length(listFeat))
    {
      #featEval[listFeat[i]] <- mean(c(featEval[listFeat[i]], pop[[modIndice]][[clf$params$evalToFit]]), na.rm = TRUE)
      featEval[listFeat[i]] <- mean(c(featEval[listFeat[i]], pop[[modIndice]]$fit_), na.rm = TRUE) # to take into account the penalized fit
      featEval <- featEval
      i <- i + 1
    }
  }
  
  if(any(featEval < 0 & !is.na(featEval)))
  {
    stop("rankFeatures, negative evaluation")
  }
  
  return(featEval)
}



mutate2 <- function(X, y, clf, pop, selection, featEval) 
{
  mutated_pop <- pop                                          ###
  all_genes <- (1:clf$params$size_world) # all genes       ######
                                                        #####################################
  featEval <- rankFeatures(X, y, clf, pop, featEval) ##################### Work to do here ##
                                                        #####################################
  #print(summary(featEval))                                ######
  if(clf$params$parallel.local)                               ###
  { # if // computing
    res <- foreach(i = selection)  %dorng% 
      {
        individual_to_be_mutated <- pop[[i]]
        mutated_indiv <- clf$functions$mutator(X, y, clf, pop, individual_to_be_mutated, all_genes, featEval)
        return(mutated_indiv)
      }
    mutated_pop[selection] <- res 
  } else
  {
    for (i in selection)
    {
      individual_to_be_mutated <- pop[[i]]
      mutated_pop[[i]] <- clf$functions$mutator(X, y, clf, pop, individual_to_be_mutated, all_genes, featEval)
    }
  }
  
  return(mutated_pop)
}



crossingIndividual_v1 <- function(X, y, clf, parent1, parent2) #not working for the moment
{
  parents_gene_reservoir <- sort(unique(c(parent1$indices_, parent2$indices_)))
  # unique vars
  
  if(parent1$eval.sparsity != parent2$eval.sparsity)
  {
    child.sparsity <- sample(parent1$eval.sparsity:parent2$eval.sparsity, 1)
  } else
  {
    child.sparsity <- parent1$eval.sparsity - sample(0:parent1$eval.sparsity, 1)
  }
  
  indices <- sort(parents_gene_reservoir[sample(x = length(parents_gene_reservoir), 
                                     size = child.sparsity, replace = FALSE)])
  return(individual(X, y, clf, ind = indices))
}

crossingIndividual_v2 <- function(X, y, clf, parent1, parent2)
{
  size <- nrow(X)
  parent1_dense <- individualGetDenseVec(parent1, size)
  parent2_dense <- individualGetDenseVec(parent2, size)
  
  child <- c(parent1_dense[1:(size/2)], parent2_dense[(size/2):size])
  
  ###### Controle de la sparsité ######
  sparsity <- length(which(child != 0))
  if(sparsity < clf$params$sparsity.min)
  {
    reservoir <- parent1_dense + parent2_dense
    reservoir <- which(reservoir != 0)
    genes2add <- sample(reservoir, clf$params$sparsity.min - sparsity) ### A changer si ça marche
                                                                       ### afin de garder une sparsité homogène
    child[genes2add] <- sample(c(-1, 1), clf$params$sparsity.min - sparsity, replace = TRUE)
  } else if(sparsity > clf$params$sparsity.max)
  {
    reservoir <- which(child != 0)
    genes2rm  <- sample(reservoir, sparsity - clf$params$sparsity.max) ### Pareil
    child[genes2rm] <- 0
  }
  
  return(individual(X, y, clf, coeffs = child))
}

crossingIndividual_v3 <- function(X, y, clf, parent1, parent2)
{
  # Get the all the features in the parents
  feat_reserv <- c(parent1$indices_, parent2$indices_)
  feat_reserv <- feat_reserv[order(feat_reserv)]
  
  if(length(unique(feat_reserv)) < max(parent1$eval.sparsity,parent2$eval.sparsity))
  {
    if(length(unique(parent1$indices_)) < parent1$eval.sparsity)
    {
      stop(paste("Parent1 has sparsity", parent1$eval.sparsity, "but only", 
                 length(unique(parent1$indices_)), "unique features"))
    } else if(length(unique(parent2$indices_)) < parent2$eval.sparsity)
    {
      stop(paste("Parent2 has sparsity", parent2$eval.sparsity, "but only", 
                 length(unique(parent2$indices_)), "unique features"))
    } else
    {
      stop("One of the parents already has a feature list with a few identical values")
    }
  }
  
  if(length(unique(feat_reserv)) > 1)
  {
    if(length(feat_reserv) != length(unique(feat_reserv)))
    {
      prob <- lapply(unique(feat_reserv), 
                     function(x, feat_reserv) 
                       {sum(x == feat_reserv)/length(feat_reserv)}, feat_reserv)
    } else
    {
      prob <- rep(1/length(feat_reserv), length(feat_reserv))
    }
    child.spar <- sample(c(parent1$eval.sparsity,parent2$eval.sparsity), 1)
    
    if(length(unique(feat_reserv)) < child.spar)
    {
      stop(paste("Error in crossingIndividual_v3, not enought differents variables in both parents. 
                 Only", length(unique(feat_reserv)), "in the parents, but", child.spar, 
                 "required in the child\n", sep = " "))
    }
    
    child.indices <- sample(unique(feat_reserv), child.spar, prob = prob)
    child.indices <- child.indices[order(child.indices)]
    return(individual(X, y, clf, ind = child.indices))
  } else 
  {
    return(parent1)
  }
}

crossing2 <- function(X, y, clf, pop, parents) 
{
  # select the needed number of couples
  number_couples <- length(pop) - length(parents)
  old_pop <- parents
  
  couples <- list()
  for(i in 1:number_couples) 
  {
    if(length(old_pop) < 2) # If there is no more parents available in the list we start again with 
      old_pop <- parents    # the full population of parents
    selectedParents <- sample(1:length(old_pop), 2)
    couples[[i]] <- old_pop[selectedParents]
    old_pop <- old_pop[-selectedParents]
  }
  
  children <- list()
  # parallel computing
  if(clf$params$parallel.local)
  {
    children <- foreach(i = 1:number_couples)  %dorng% 
    {
      child <- clf$functions$crosser(X, y, clf, couples[[i]][[1]], couples[[i]][[2]])
      return(child)
    }
  } else
  {
    for (i in 1:number_couples)
    {
      children[[i]] <- clf$functions$crosser(X, y, clf, couples[[i]][[1]], couples[[i]][[2]])
    }
  }
  return(c(parents, children))
}

#########################
#                       #
#  ###################  #
#  #### SELECTION ####  #
#  ###################  #
#                       #
#########################

#' Does an elite selection on a population
#'
#' @description This function is a template for other selectors, it takes a population and a number of individuals to select. The result is the \code{number} bests element of the population.
selector_v1 <- function(pop, number, clf)
{
  # evaluation <- getFitPopulation(population)
  # ind <- order(evaluation,decreasing = T)[1:number]
  return(sortPopulation(pop, evalToOrder = "fit_")[1:number])
  # return(population[ind])
}

selector_v2 <- function(pop, number, clf = NULL)
{ #TODO: changer ce selecteur
  ind <- sample(x = 1:length(pop), size  =  number, replace = FALSE)
  return(pop[ind])
}

# This function realise the selection by using the selection functions given in the clf
select <- function(clf, pop)
{
  nb2BeMutated <- ceiling(clf$params$select_perc * clf$params$size_pop/100)
  nbByMethod   <- lapply(clf$params$select_percByMethod, 
                         function(x) {x * nb2BeMutated/100})
  res <- rep(NA, nb2BeMutated)
  start <- 1
  end   <- nbByMethod[[1]]
  popTMP <- pop
  for(i in 1:length(nbByMethod))
  {
    popTMP <- pop[which(!(popTMP %in% res))]
    res[start:end] <- clf$functions$selector[[i]](popTMP, nbByMethod[[i]], clf)
    start <- end +1
    end <- end + nbByMethod[[i]]
  }
  return(res)
}


# Select the best individuals for each sparsity to keep them for the next generation
selectIndividualToKeep <- function(clf, pop)
{
  spar <- populationGet_X(element2get = "eval.sparsity", toVec = TRUE, na.rm = TRUE)(pop)
  popBySparsity <- lapply(clf$params$sparsity, function(k, spar, pop) 
    {
      pop[which(spar == k)]
    }, spar, pop)
  popBySparsity <- popBySparsity[sapply(popBySparsity, function(x)
    {
      if(length(x) > 0)
      {
        TRUE
      }else
      {
        FALSE
      }
    })]
  bestBySparsity <- lapply(popBySparsity, function(pop, evalToFit) {getTheBestIndividual(pop, evalToFit)}, "fit_")
  return(bestBySparsity)
}

#####
#####
# A garder dans un coin pour visualiser les fit : 
# plot(sapply(population, function(x) x$fit_)[order(sapply(population, function(x) x$fit_))])
#
#####
#####


individualGetDenseVec <- function(individual, size)
{
  res <- rep(0, size)
  res[individual$indices_] <- individual$coeffs_
  return(res)
}

