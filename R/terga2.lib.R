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

#' Generate and Initialize a Population for Evolutionary Computation
#'
#' This function generates an initial population of individuals for evolutionary
#' algorithms, potentially using predefined parameters and sources, and
#' evaluates their fitness.
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object that includes parameters for generating and
#'   evaluating individuals. Expected to contain `params` and possibly `coeffs_`
#'   elements.
#' @param featEval An optional named vector of feature evaluation metrics. If
#'   \code{NULL}, an NA-filled vector will be created with the same names as the
#'   row names of \code{X}.
#'
#' @details The function initializes an evolutionary population by creating
#' individuals based on signs, injecting predefined individuals (if specified),
#' loading additional individuals from files, and generating new individuals
#' with random or specific feature models. After generation, the function
#' evaluates and sorts the population based on fitness.
#'
#' @return A list representing the initialized and evaluated population, sorted
#'   by fitness. Returns \code{NULL} if the population fails to meet validity
#'   checks.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list(params = list(randomSigns = FALSE, language = "bin", size_pop_random = 5, parallel.local = FALSE, in_pop = "NULL", popSourceFile = "NULL", verbose = TRUE, evolve_k1 = FALSE))
#' pop <- population2(X, y, clf)
#' }
#'
#' @import foreach
#' @import doRNG
#' @export
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


#' Generate a Sparse Individual Vector for Evolutionary Computation
#'
#' This function generates an individual vector representation for evolutionary
#' algorithms. The vector is generated with specified sparsity, based on binary
#' or ternary language, depending on the configuration in the classifier object
#' (`clf`).
#'
#' @param clf A classifier object containing parameters that guide the vector
#'   generation, including `language`, `size_world`, `sparsity.mean`,
#'   `sparsity.min`, and `sparsity.max`.
#' @param signs An optional vector of signs used to set or adjust the values in
#'   the generated individual vector. Defaults to \code{NULL}.
#'
#' @details The function creates a vector of length `clf$params$size_world`,
#' with elements initialized based on a specified probability distribution. For
#' binary language, the vector elements are either 1 or 0, while for ternary
#' language, they can be -1, 0, or 1. The function ensures that the vector
#' respects specified minimum and maximum sparsity levels
#' (`clf$params$sparsity.min` and `clf$params$sparsity.max`), adjusting the
#' number of non-zero elements as needed.
#'
#' @return A numeric vector of length `clf$params$size_world`, representing an
#'   individual with the desired sparsity and language settings.
#'
#' @examples
#' \dontrun{
#' clf <- list(params = list(language = "ternary", size_world = 100, sparsity.mean = 0.1, sparsity.min = 5, sparsity.max = 20))
#' individual <- individual_vec_v1(clf)
#' }
#'
#' @export
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


#' Generate a Sparse Individual Vector with Specified Density
#'
#' This function generates an individual vector with a specific sparsity level,
#' based on parameters defined in a classifier object (`clf`). It supports
#' binary and ternary language configurations and allows for sign customization.
#'
#' @param clf A classifier object that contains parameters for generating the
#'   vector, such as `size_world`, `sparsity.min`, `sparsity.max`,
#'   `sparsity.mean`, and `language`.
#' @param signs An optional vector of signs used to populate the individual
#'   vector at specified indices. If \code{NULL}, the function randomly assigns
#'   values based on the language setting in \code{clf}.
#'
#' @details The function constructs a vector of length `clf$params$size_world`
#' and assigns values based on either binary or ternary language, as specified
#' in `clf$params$language`. The sparsity level is controlled by `sparsity.min`
#' and `sparsity.max`, or by `sparsity.mean` if the minimum and maximum sparsity
#' levels are equal. If `signs` is provided, it will populate the vector at
#' selected indices.
#'
#' If `clf$params$size_world` is not set, a warning is issued (if warnings are
#' enabled in `clf`) and the function returns \code{NULL}.
#'
#' @return A numeric vector of length `clf$params$size_world` with elements set
#'   based on sparsity and language parameters. Returns \code{NULL} if
#'   `size_world` is not properly set.
#'
#' @examples
#' \dontrun{
#' clf <- list(params = list(language = "bin", size_world = 50, sparsity.mean = 10, sparsity.min = 5, sparsity.max = 15, warnings = TRUE))
#' individual <- individual_vec_v2(clf)
#' }
#'
#' @export
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


#' Generate Sparse Individual Vector with Indices and Coefficients
#'
#' This function generates a sparse individual vector representation for
#' evolutionary algorithms, returning a list with the indices and coefficients
#' of non-zero elements. The values and sparsity level are based on the
#' parameters defined in a classifier object (`clf`).
#'
#' @param clf A classifier object containing parameters for generating the
#'   vector, including `size_world`, `sparsity.min`, `sparsity.max`, and
#'   `language`.
#' @param signs An optional vector of signs (not used in this function but
#'   included for consistency with other versions).
#'
#' @details The function randomly selects a number of non-zero elements within
#' the range specified by `sparsity.min` and `sparsity.max`. It then assigns
#' values to these elements based on `clf$params$language`: for binary language,
#' only 1s are used; for ternary language, values are -1 or 1.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{indices}{A vector of indices corresponding to non-zero elements in the generated vector.}
#'   \item{coeffs}{A vector of coefficients (1 or -1) at the specified indices.}
#' }
#'
#' @examples
#' \dontrun{
#' clf <- list(params = list(language = "ternary", size_world = 50, sparsity.min = 5, sparsity.max = 15))
#' individual <- individual_vec_v3(clf)
#' print(individual)
#' }
#'
#' @export
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


#' Select Elite Individuals in a Population
#'
#' This function selects a specified number of elite individuals from a
#' population by marking them as selected. The population is sorted by fitness,
#' and the top individuals are tagged based on the `nbToSelect` parameter.
#'
#' @param clf A classifier object that may contain parameters for selection
#'   (though it is not directly used in this function).
#' @param pop A list representing the population, where each individual is an
#'   element with a `selected` attribute and a `fit_` fitness attribute.
#' @param nbToSelect An integer specifying the number of individuals to mark as
#'   elite in the population.
#'
#' @details The function sorts the population by the `fit_` attribute
#' (presumably representing fitness) in descending order. It then iterates
#' through the sorted list, marking the top `nbToSelect` individuals as
#' `selected` by setting their `selected` attribute to `TRUE`. If the population
#' is empty, a warning is issued and the function returns \code{NULL}.
#'
#' @return The modified population list with the `selected` attribute set to
#'   \code{TRUE} for the top elite individuals.
#'
#' @examples
#' \dontrun{
#' pop <- list(
#'   list(fit_ = 0.9, selected = FALSE),
#'   list(fit_ = 0.8, selected = FALSE),
#'   list(fit_ = 0.7, selected = FALSE)
#' )
#' clf <- list() # Placeholder for classifier object
#' pop <- tag_SelectElite(clf, pop, nbToSelect = 2)
#' # Check which individuals are selected
#' sapply(pop, function(ind) ind$selected)
#' }
#'
#' @export
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


#' Randomly Select Individuals in a Population
#'
#' This function randomly selects a specified number of unselected individuals
#' from a population and marks them as selected. The function can operate in
#' parallel if specified in the classifier parameters.
#'
#' @param clf A classifier object that contains parameters for selection,
#'   including `params$warnings` (to control warnings) and
#'   `params$parallel.local` (for parallel execution).
#' @param pop A list representing the population, where each individual has a
#'   `selected` attribute.
#' @param nbToSelect An integer specifying the number of individuals to randomly
#'   select and mark as selected.
#'
#' @details The function identifies individuals in the population who have not
#' been selected (i.e., their `selected` attribute is `FALSE`). It then randomly
#' selects up to `nbToSelect` of these unselected individuals and sets their
#' `selected` attribute to `TRUE`. If fewer than `nbToSelect` unselected
#' individuals are available, it selects all unselected individuals.
#'
#' If parallel processing is enabled (`clf$params$parallel.local` is `TRUE`),
#' selection is performed in parallel using `foreach` and `doRNG`. If the
#' population is empty, the function issues a warning (if `clf$params$warnings`
#' is `TRUE`) and returns \code{NULL}.
#'
#' @return The modified population list with randomly selected individuals
#'   marked as `selected = TRUE`.
#'
#' @examples
#' \dontrun{
#' pop <- list(
#'   list(fit_ = 0.9, selected = FALSE),
#'   list(fit_ = 0.8, selected = FALSE),
#'   list(fit_ = 0.7, selected = FALSE)
#' )
#' clf <- list(params = list(warnings = TRUE, parallel.local = FALSE))
#' pop <- tag_SelectRandom(clf, pop, nbToSelect = 2)
#' # Check which individuals are selected
#' sapply(pop, function(ind) ind$selected)
#' }
#'
#' @import foreach
#' @import doRNG
#' @export
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


#' Pair Individuals in Population as Mating Couples
#'
#' This function pairs individuals in a population to form mating couples,
#' assigning each pair a `mate` attribute. Each individual can be matched with
#' one other individual from the specified parents list.
#'
#' @param pop A list representing the population, where each individual has a
#'   `mate` attribute to indicate their paired individual.
#' @param parents A vector of indices or identifiers for individuals in the
#'   population that are eligible to be paired.
#'
#' @details The function randomly pairs individuals from the `parents` list. For
#'   each couple, it assigns the `mate` attribute of one individual to the
#'   other’s index, and vice versa. If an individual’s `mate` attribute is
#'   positive, it indicates the index of their partner. If it’s negative, it
#'   indicates their partner’s index with a negative sign.
#'
#'   If the population is empty, the function issues a warning and returns
#'   \code{NULL}.
#'
#' @return The modified population list, with `mate` attributes assigned for
#'   paired individuals.
#'
#' @examples
#' \dontrun{
#' pop <- list(
#'   list(fit_ = 0.9, selected = FALSE, mate = NULL),
#'   list(fit_ = 0.8, selected = FALSE, mate = NULL),
#'   list(fit_ = 0.7, selected = FALSE, mate = NULL),
#'   list(fit_ = 0.6, selected = FALSE, mate = NULL)
#' )
#' parents <- c(1, 2, 3, 4)
#' pop <- tag_Couples(pop, parents)
#' # Check mates for each individual
#' sapply(pop, function(ind) ind$mate)
#' }
#'
#' @export
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


#' Tag Individuals for Mutation in a Population
#'
#' This function tags a specified number of individuals in the population for
#' mutation, while optionally protecting certain individuals from being
#' selected. The tagged individuals are marked with a `toBeMutated` attribute
#' set to `TRUE`.
#'
#' @param pop A list representing the population, where each individual has a
#'   `toBeMutated` attribute.
#' @param mutate_size An integer specifying the number of individuals to mark
#'   for mutation.
#' @param protected An optional vector of indices representing individuals that
#'   should not be tagged for mutation. Defaults to \code{NULL}.
#'
#' @details The function selects `mutate_size` individuals from the population,
#' excluding any indices specified in `protected`. If `mutate_size` exceeds the
#' number of non-protected individuals, the function issues a warning and
#' reduces the mutation count to the maximum available.
#'
#' If the number of protected individuals is too large to allow sufficient
#' mutations, the function removes some individuals from the `protected` list as
#' necessary to meet the mutation quota.
#'
#' If the population is empty, the function issues a warning (if
#' `clf$params$warnings` is `TRUE`) and returns \code{NULL}.
#'
#' @return The modified population list with `toBeMutated` attributes set to
#'   `TRUE` for the selected individuals.
#'
#' @examples
#' \dontrun{
#' pop <- list(
#'   list(fit_ = 0.9, toBeMutated = FALSE),
#'   list(fit_ = 0.8, toBeMutated = FALSE),
#'   list(fit_ = 0.7, toBeMutated = FALSE),
#'   list(fit_ = 0.6, toBeMutated = FALSE)
#' )
#' protected <- c(1, 2)
#' pop <- tag_ToBeMutated(pop, mutate_size = 2, protected = protected)
#' # Check which individuals are tagged for mutation
#' sapply(pop, function(ind) ind$toBeMutated)
#' }
#'
#' @export
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


#' Select and Tag Individuals in a Population for Evolutionary Operations
#'
#' This function selects individuals in a population for evolutionary
#' operations, including elite and random selection, pairing individuals for
#' mating, and tagging some individuals for mutation. The selection process is
#' based on specified percentages in the classifier object (`clf`).
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object containing parameters for selection, including
#'   `select_perc`, `select_percByMethod`, `size_pop`, and `mutate_size`.
#' @param pop A list representing the population, where each individual has
#'   attributes including `mate` and `toBeMutated`.
#'
#' @details The function begins by sorting the population based on fitness
#' (`fit_` attribute). It then calculates the number of individuals to be
#' selected for mutation based on `select_perc` in `clf`. Selection is split
#' into two methods, elite and random, as specified by `select_percByMethod`,
#' and individuals are tagged for mating within these groups.
#'
#' Following selection, individuals are paired for mating, with each selected
#' individual having a `mate` attribute set. The function then calls
#' `tag_ToBeMutated` to mark some individuals for mutation based on
#' `mutate_size` from `clf`.
#'
#' @return The modified population list, where selected individuals are tagged
#'   with `mate` attributes for mating and `toBeMutated` attributes for
#'   mutation.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list(params = list(select_perc = 50, select_percByMethod = c(30, 20), size_pop = 10, mutate_size = 10))
#' pop <- list(
#'   list(fit_ = 0.9, mate = list(), toBeMutated = FALSE),
#'   list(fit_ = 0.8, mate = list(), toBeMutated = FALSE),
#'   # Add additional individuals as needed
#' )
#' pop <- tag_select(X, y, clf, pop)
#' # Check which individuals are selected for mating and mutation
#' sapply(pop, function(ind) ind$mate)
#' sapply(pop, function(ind) ind$toBeMutated)
#' }
#'
#' @export
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



#' Retrieve Selected Individuals as Parents
#'
#' This function identifies and returns the indices of individuals in the
#' population who are marked as selected, which can be used as parents in
#' evolutionary processes.
#'
#' @param pop A list representing the population, where each individual has a
#'   `selected` attribute indicating if they are selected.
#'
#' @details The function uses `populationGet_X` to retrieve individuals with the
#' `selected` attribute set to `TRUE`. It returns the indices of these
#' individuals, allowing them to be used as parents in subsequent evolutionary
#' steps.
#'
#' @return A vector of indices corresponding to individuals in the population
#'   who are marked as selected.
#'
#' @examples
#' \dontrun{
#' pop <- list(
#'   list(fit_ = 0.9, selected = TRUE),
#'   list(fit_ = 0.8, selected = FALSE),
#'   list(fit_ = 0.7, selected = TRUE)
#' )
#' parents <- get_Parents(pop)
#' print(parents) # Should return indices of selected individuals
#' }
#'
#' @export
get_Parents <- function(pop)
{
  return(which(populationGet_X(element2get = "selected", toVec = TRUE, na.rm = TRUE)(pop)))
}


#' Retrieve Individuals Marked for Mutation
#'
#' This function identifies and returns the indices of individuals in the
#' population who are marked for mutation, as indicated by the `toBeMutated`
#' attribute.
#'
#' @param pop A list representing the population, where each individual has a
#'   `toBeMutated` attribute indicating if they are marked for mutation.
#'
#' @details The function uses `populationGet_X` to retrieve individuals with the
#' `toBeMutated` attribute set to `TRUE`. It returns the indices of these
#' individuals, allowing them to be processed for mutation in subsequent
#' evolutionary steps.
#'
#' @return A vector of indices corresponding to individuals in the population
#'   who are marked for mutation.
#'
#' @examples
#' \dontrun{
#' pop <- list(
#'   list(fit_ = 0.9, toBeMutated = TRUE),
#'   list(fit_ = 0.8, toBeMutated = FALSE),
#'   list(fit_ = 0.7, toBeMutated = TRUE)
#' )
#' toBeMutated <- get_IndividualToBeMutated(pop)
#' print(toBeMutated) # Should return indices of individuals marked for mutation
#' }
#'
#' @export
get_IndividualToBeMutated <- function(pop)
{
  return(which(populationGet_X(element2get = "toBeMutated", toVec = TRUE, na.rm = TRUE)(pop)))
}


#' Reset Tags for All Individuals in a Population
#'
#' This function resets specific tags for all individuals in a population,
#' including `selected`, `toBeMutated`, and `mate` attributes. It can be used to
#' clear or initialize these tags before further processing.
#'
#' @param pop A list representing the population, where each individual is
#'   expected to have attributes `selected`, `toBeMutated`, and `mate`.
#' @param selected A logical value to set for the `selected` attribute of each
#'   individual. Defaults to \code{FALSE}.
#' @param toBeMutated A logical value to set for the `toBeMutated` attribute of
#'   each individual. Defaults to \code{FALSE}.
#' @param mate A value to set for the `mate` attribute of each individual.
#'   Defaults to \code{-1}.
#'
#' @details The function iterates through each individual in the population and
#' sets their `selected`, `toBeMutated`, and `mate` attributes to the specified
#' values. This function is useful for resetting the population's status before
#' selection and mutation operations.
#'
#' If `pop` is \code{NULL}, the function issues a warning (if
#' `clf$params$warnings` is enabled) and returns \code{NULL}. If `pop` is not a
#' list, an error is raised.
#'
#' @return The modified population list with updated `selected`, `toBeMutated`,
#'   and `mate` attributes for each individual.
#'
#' @examples
#' \dontrun{
#' pop <- list(
#'   list(fit_ = 0.9, selected = TRUE, toBeMutated = TRUE, mate = 2),
#'   list(fit_ = 0.8, selected = TRUE, toBeMutated = FALSE, mate = 3),
#'   list(fit_ = 0.7, selected = FALSE, toBeMutated = TRUE, mate = -1)
#' )
#' pop <- resetTags(pop, selected = FALSE, toBeMutated = FALSE, mate = -1)
#' # Check updated population
#' str(pop)
#' }
#'
#' @export
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



#' Evolve a Population in an Evolutionary Algorithm
#'
#' This function performs one cycle of evolution on a population, including
#' evaluation, selection, crossover, and mutation steps, based on the specified
#' classifier parameters. It returns the evolved population after these
#' operations.
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object containing parameters for evolution, including
#'   options for crossover, mutation, and verbose output.
#' @param pop A list representing the population, where each individual has
#'   attributes like `fit_` for fitness.
#' @param featEval A vector or list containing feature evaluation metrics, used
#'   in mutation to assess feature significance.
#'
#' @details The function first evaluates the fitness of the current population
#' if it has not been done. It then selects the best individual and proceeds
#' through the evolutionary steps:
#'
#' 1. **Evaluation**: If fitness values are missing, the population is
#' evaluated. 2. **Selection & Crossover**: If crossover is enabled
#' (`clf$params$cross`), parents are selected and crossover is performed to
#' generate new individuals. 3. **Mutation**: If mutation is enabled
#' (`clf$params$mutate`), some individuals are selected for mutation based on
#' `mutate_size`, and mutations are applied using `featEval`.
#'
#' After each step, the population is re-evaluated, and the best individual is
#' updated. If `verbose` is enabled, details about the best individual are
#' printed.
#'
#' @return The evolved population list with updated individuals and fitness
#'   values.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list(params = list(cross = TRUE, mutate = TRUE, mutate_size = 10, verbose = TRUE))
#' pop <- list(
#'   list(fit_ = NA, indices_ = 1:5),
#'   list(fit_ = NA, indices_ = 6:10)
#' )
#' featEval <- runif(10)
#' evolved_pop <- evolve1(X, y, clf, pop, featEval)
#' }
#'
#' @export
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



#' Evolve a Population with Advanced Selection, Crossover, and Mutation
#'
#' This function performs one generation of evolution on a population,
#' implementing selection, crossover, and mutation with customizable methods. It
#' returns the evolved population after these operations.
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object containing parameters for evolutionary
#'   operations, including options for selection percentages, mutation size,
#'   population size, parallelization, and debugging.
#' @param pop A list representing the population, where each individual has
#'   attributes like `fit_`, `selected`, `mate`, and `toBeMutated`.
#' @param featEval A vector or list containing feature evaluation metrics, used
#'   to guide mutation operations.
#' @param generation An integer representing the current generation number in
#'   the evolutionary process.
#'
#' @details The `evolve2m` function implements a detailed evolutionary cycle
#' with the following steps:
#'
#' 1. **Population Sorting**: The population is sorted by fitness, and
#' individuals are evaluated if necessary. 2. **Tagging for Selection**:
#' Individuals are tagged for elite and random selection based on selection
#' percentages in `clf$params$select_percByMethod`. 3. **Tagging for Mutation**:
#' Some individuals are tagged for mutation based on `mutate_size`. 4.
#' **Parallel and Sequential Processing**: If parallel processing is enabled
#' (`clf$params$parallel.local`), crossover and mutation are performed in
#' parallel; otherwise, they are processed sequentially. 5. **Mating and
#' Mutation**: Selected individuals are paired for mating, and some are mutated
#' based on the feature evaluation metrics (`featEval`). 6. **Population Pruning
#' and Expansion**: The evolved population is pruned to retain the top 2/3 of
#' individuals by fitness, then expanded by adding new individuals or retaining
#' some from the previous population if necessary.
#'
#' The population is further cleaned, sorted, and prepared for the next
#' generation. If `debug` mode is enabled in `clf`, various messages are printed
#' to track the evolution process.
#'
#' @return The evolved population list, containing updated individuals with
#'   modified `fit_`, `selected`, `mate`, and `toBeMutated` attributes.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list(params = list(
#'   select_perc = 50, select_percByMethod = c(30, 20), size_pop = 10,
#'   mutate_size = 10, cross = TRUE, mutate = TRUE, parallel.local = FALSE, debug = TRUE
#' ))
#' pop <- list(
#'   list(fit_ = NA, indices_ = 1:5, selected = FALSE, toBeMutated = FALSE, mate = NULL),
#'   list(fit_ = NA, indices_ = 6:10, selected = FALSE, toBeMutated = FALSE, mate = NULL)
#' )
#' featEval <- runif(10)
#' evolved_pop <- evolve2m(X, y, clf, pop, featEval, generation = 1)
#' }
#'
#' @export
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



#' Evolve a Population with Advanced Mating and Mutation
#'
#' This function performs one generation of evolution on a population, including
#' advanced mating and mutation, with optional parallel processing. It returns
#' the evolved population after these operations.
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object containing parameters for evolution, including
#'   options for mutation size, parallel processing, and debugging.
#' @param pop A list representing the population, where each individual has
#'   attributes like `mate` and `toBeMutated`.
#' @param featEval A vector or list containing feature evaluation metrics, used
#'   to guide mutation.
#'
#' @details The `evolve3m` function performs an evolutionary cycle on the
#' population with the following steps:
#'
#' 1. **Population Evaluation**: Each individual in the population is evaluated
#' if not already done, and `mate` and `toBeMutated` attributes are initialized.
#' 2. **Tagging for Selection and Mating**: Individuals are tagged based on
#' selection criteria from `clf`. 3. **Crossover and Mutation**: For individuals
#' with mates, offspring are generated using the `crosser` function. If an
#' individual is marked for mutation, the `mutator` function is applied.
#'
#' The function supports parallel processing, where crossover and mutation steps
#' are performed concurrently if `clf$params$parallel.local` is `TRUE`. After
#' crossover and mutation, children and mutated individuals are combined into a
#' new population. If `debug` mode is enabled, messages track the evolution
#' process.
#'
#' @return The evolved population list, including new offspring and mutated
#'   individuals.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list(params = list(
#'   mutate_size = 10, parallel.local = FALSE, debug = TRUE
#' ))
#' pop <- list(
#'   list(fit_ = NA, indices_ = 1:5, selected = FALSE, toBeMutated = FALSE, mate = NULL),
#'   list(fit_ = NA, indices_ = 6:10, selected = FALSE, toBeMutated = FALSE, mate = NULL)
#' )
#' featEval <- runif(10)
#' evolved_pop <- evolve3m(X, y, clf, pop, featEval)
#' }
#'
#' @export
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


#' Mutate an Individual's Genes in the Population
#'
#' This function mutates a specified individual's genes by altering a subset of
#' them according to a defined mutation rate and gene reservoir.
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object containing parameters for mutation, including
#'   `mutate_rate` and `size_world`.
#' @param pop The current population list, used as context for the mutation but
#'   not directly modified by this function.
#' @param individual_to_be_mutated A list representing the individual to be
#'   mutated, with attributes `eval.sparsity` (sparsity level) and `indices_`
#'   (genes).
#' @param all_genes A vector of all possible gene indices, representing the
#'   complete set of genes that could be selected during mutation.
#' @param featEval A vector or list of feature evaluation metrics, used to guide
#'   the mutation.
#'
#' @details The function calculates the number of genes to mutate (`perc`) based
#' on the individual’s sparsity level and the mutation rate defined in `clf`. If
#' the number of genes to mutate exceeds the available gene reservoir, the
#' mutation rate is adjusted accordingly. Genes for mutation are chosen from
#' `all_genes` and replaced in the individual's genome.
#'
#' The mutation process includes:
#' - Calculating mutation degree based on mutation rate and sparsity.
#' - Selecting genes for mutation and replacing them with new genes from the reservoir.
#'
#' @return A mutated individual list, representing the modified genome of the
#'   input individual.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list(params = list(mutate_rate = 20, size_world = 100))
#' individual_to_be_mutated <- list(eval.sparsity = 5, indices_ = sample(1:100, 5))
#' all_genes <- 1:100
#' featEval <- runif(10)
#' mutated_individual <- mutator_v1(X, y, clf, pop = NULL, individual_to_be_mutated, all_genes, featEval)
#' }
#'
#' @export
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


#' Mutate an Individual's Genes Based on Feature Evaluation
#'
#' This function mutates a specified individual's genes, choosing mutation
#' actions based on feature evaluation probabilities rather than a fixed
#' mutation rate. Actions include modifying, adding, or removing genes,
#' depending on mutation conditions and individual sparsity.
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object containing parameters for mutation, including
#'   `current_seed`, `language`, and `sparsity`.
#' @param pop The current population list, used as context for mutation but not
#'   directly modified by this function.
#' @param individual_to_be_mutated A list representing the individual to mutate,
#'   with attributes `indices_` (genes) and `coeffs_` (gene coefficients).
#' @param all_genes A vector of all possible gene indices, representing the
#'   complete set of genes that could be selected during mutation.
#' @param featEval A vector or list of feature evaluation metrics, used to guide
#'   mutation probabilities for selecting features.
#'
#' @details The mutation process is guided by `featEval`, which assigns a
#' probability for each feature to be selected for mutation based on its
#' evaluation score. The mutation action is chosen based on the individual’s
#' sparsity and classifier settings:
#'
#' - **Modify**: Inverts the coefficient of a randomly selected gene within the individual.
#' - **Add**: Adds a new gene to the individual, drawn from `all_genes` based on feature evaluation probabilities.
#' - **Remove**: Deletes a randomly selected gene from the individual.
#'
#' The mutation action is determined based on the individual’s current sparsity
#' and the classifier's settings. If no suitable genes remain for mutation, the
#' function returns the individual unchanged.
#'
#' @return A mutated individual list with updated `indices_` and `coeffs_`,
#'   representing the modified genome of the input individual.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list(params = list(current_seed = 42, language = "ternary", sparsity = c(5, 10)))
#' individual_to_be_mutated <- list(indices_ = sample(1:100, 5), coeffs_ = rep(1, 5))
#' all_genes <- 1:100
#' featEval <- runif(10)
#' mutated_individual <- mutator_v2(X, y, clf, pop = NULL, individual_to_be_mutated, all_genes, featEval)
#' }
#'
#' @export
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



#' Rank Features Based on Model Fitness
#'
#' This function evaluates and ranks features in the population based on their
#' presence in individual models and their associated fitness scores. It updates
#' a feature evaluation metric (`featEval`) to reflect the average fitness of
#' models in which each feature appears.
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object containing parameters for evaluation,
#'   including `evalToFit`, which specifies the fitness attribute in each model.
#' @param pop A list representing the population, where each individual model
#'   contains `indices_` (genes or features used in the model) and `fit_`
#'   (fitness score).
#' @param featEval A numeric vector of feature evaluation scores, where each
#'   element corresponds to a feature in \code{X}. This vector is updated based
#'   on the average fitness of models that include each feature.
#'
#' @details The function iterates over each model in the population (`pop`) and
#' each feature (gene) in the model. For each feature, it updates `featEval` by
#' averaging the current feature score with the fitness score of the model in
#' which it appears. The result is an updated ranking of features based on their
#' presence in high-fitness models.
#'
#' The function checks for negative values in `featEval` and raises an error if
#' any are found, as negative values indicate an issue in the evaluation
#' process.
#'
#' @return A numeric vector `featEval`, with updated values reflecting the
#'   average fitness scores of models that include each feature.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list(params = list(evalToFit = "fit_"))
#' pop <- list(
#'   list(indices_ = c(1, 3, 5), fit_ = 0.8),
#'   list(indices_ = c(2, 4, 6), fit_ = 0.7)
#' )
#' featEval <- rep(0, 10)
#' updatedFeatEval <- rankFeatures(X, y, clf, pop, featEval)
#' }
#'
#' @export
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


#' Apply Mutation to a Selected Subset of the Population
#'
#' This function performs mutation on a selected subset of individuals in a
#' population, using feature evaluations to guide mutation decisions. It returns
#' the mutated population.
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object containing parameters for mutation, including
#'   `size_world` (total number of genes) and `parallel.local` (to enable
#'   parallel computation).
#' @param pop A list representing the population, where each individual has
#'   attributes relevant to mutation.
#' @param selection A vector of indices representing the individuals in `pop`
#'   that should be mutated.
#' @param featEval A numeric vector of feature evaluation scores, which will be
#'   updated to reflect the mutation's effect on feature importance.
#'
#' @details The `mutate2` function first ranks features using `rankFeatures`,
#' updating `featEval` based on model fitness in the population. It then
#' iterates over the selected individuals (`selection`), applying mutation to
#' each. If parallel processing is enabled (`clf$params$parallel.local`),
#' mutations are performed concurrently.
#'
#' The mutation process uses the classifier's `mutator` function to alter gene
#' indices or coefficients in each selected individual based on feature
#' evaluation probabilities.
#'
#' @return A list representing the mutated population, with the individuals in
#'   `selection` modified according to mutation operations.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list(params = list(size_world = 10, parallel.local = FALSE),
#'             functions = list(mutator = mutator_v1))
#' pop <- list(
#'   list(indices_ = c(1, 3, 5), fit_ = 0.8),
#'   list(indices_ = c(2, 4, 6), fit_ = 0.7)
#' )
#' selection <- 1:2
#' featEval <- rep(0, 10)
#' mutated_pop <- mutate2(X, y, clf, pop, selection, featEval)
#' }
#'
#' @export
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


#' Create a New Individual by Crossing Genes from Two Parents
#'
#' This function generates a new individual by combining genes from two parent
#' individuals. It samples a subset of genes from the parents' combined gene
#' pool, with the child’s gene sparsity determined based on the parents'
#' sparsity levels.
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object containing parameters for generating a new
#'   individual.
#' @param parent1 A list representing the first parent, which includes
#'   `indices_` (genes) and `eval.sparsity` (gene sparsity).
#' @param parent2 A list representing the second parent, also with `indices_`
#'   and `eval.sparsity`.
#'
#' @details The function combines the genes of `parent1` and `parent2` into a
#' unique pool and selects a subset for the child. The child’s sparsity level is
#' randomly determined based on the range between the two parents'
#' `eval.sparsity`. If the parents have identical sparsity, a subset is sampled
#' within the parent's sparsity limit.
#'
#' **Note**: The function currently has limitations and may need further development to handle edge cases and differing gene pools effectively.
#'
#' @return A list representing the new individual, with a subset of genes
#'   (`indices`) inherited from the parents.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list() # Placeholder for classifier settings
#' parent1 <- list(indices_ = c(1, 3, 5, 7), eval.sparsity = 4)
#' parent2 <- list(indices_ = c(2, 4, 6, 8), eval.sparsity = 4)
#' child <- crossingIndividual_v1(X, y, clf, parent1, parent2)
#' print(child)
#' }
#'
#' @export
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



#' Create a New Individual by Combining Gene Vectors from Two Parents
#'
#' This function generates a new individual by combining dense gene vectors from
#' two parent individuals. It applies sparsity control to ensure the child’s
#' gene vector meets specified minimum and maximum sparsity constraints.
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object containing parameters for generating a new
#'   individual, including `sparsity.min` and `sparsity.max` to control
#'   sparsity.
#' @param parent1 A list representing the first parent, containing attributes
#'   necessary for obtaining a dense gene vector.
#' @param parent2 A list representing the second parent, also containing
#'   attributes necessary for obtaining a dense gene vector.
#'
#' @details The function creates a dense gene vector (`child`) by taking half of
#' the gene vector from `parent1` and half from `parent2`. After combining, it
#' applies sparsity control:
#'
#' - **If the child’s sparsity is below the minimum** (`sparsity.min`), additional genes are randomly added from the parents' combined gene reservoir until the minimum sparsity is reached.
#' - **If the child’s sparsity exceeds the maximum** (`sparsity.max`), excess genes are randomly removed until the sparsity meets the maximum threshold.
#'
#' This approach ensures that the child’s gene vector adheres to the specified
#' sparsity constraints, creating a balanced and controlled gene distribution.
#'
#' @return A list representing the new individual, with a gene vector created by
#'   combining and adjusting genes from both parents.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list(params = list(sparsity.min = 3, sparsity.max = 6))
#' parent1 <- list() # Placeholder, assuming `individualGetDenseVec` can process this
#' parent2 <- list() # Placeholder, assuming `individualGetDenseVec` can process this
#' child <- crossingIndividual_v2(X, y, clf, parent1, parent2)
#' print(child)
#' }
#'
#' @export
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



#' Create a New Individual by Sampling Genes from Two Parents with Controlled
#' Sparsity
#'
#' This function generates a new individual by sampling genes from the combined
#' gene pool of two parents, with controlled sparsity based on each parent's
#' sparsity level. It ensures sufficient unique features in the gene pool to
#' meet the child's sparsity requirements.
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object containing parameters for generating a new
#'   individual.
#' @param parent1 A list representing the first parent, containing attributes
#'   `indices_` (gene indices) and `eval.sparsity` (sparsity level).
#' @param parent2 A list representing the second parent, also containing
#'   `indices_` and `eval.sparsity`.
#'
#' @details The function starts by creating a combined feature reservoir
#' (`feat_reserv`) from the gene indices of both parents. It checks for adequate
#' unique features to satisfy the maximum sparsity of either parent. If either
#' parent lacks enough unique genes to meet its declared sparsity, the function
#' raises an error.
#'
#' The child’s sparsity level is randomly chosen from the sparsity of one of the
#' parents. If the feature reservoir has duplicates, probabilities are assigned
#' based on feature frequency. The function then samples a subset of genes from
#' the reservoir for the child, ensuring the child’s sparsity matches the chosen
#' level.
#'
#' @return A list representing the new individual, with a subset of genes
#'   sampled from both parents. If insufficient unique genes are available, the
#'   function raises an error.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list() # Placeholder for classifier settings
#' parent1 <- list(indices_ = c(1, 3, 5, 7), eval.sparsity = 4)
#' parent2 <- list(indices_ = c(2, 4, 6, 8), eval.sparsity = 4)
#' child <- crossingIndividual_v3(X, y, clf, parent1, parent2)
#' print(child)
#' }
#'
#' @export
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


#' Generate Offspring by Pairing Parents in a Population
#'
#' This function generates a new set of offspring by pairing parents from a
#' population. It randomly selects pairs of parents to create offspring through
#' a specified crossover function. The result is a combined list of the original
#' parents and their offspring.
#'
#' @param X A matrix or data frame of feature values, where each row represents
#'   a feature and each column represents a sample.
#' @param y A response vector or target variable for supervised learning. The
#'   length of \code{y} should match the number of columns in \code{X}.
#' @param clf A classifier object containing parameters for crossover, including
#'   `parallel.local` (to enable parallel computation) and a `crosser` function.
#' @param pop A list representing the current population, where each individual
#'   can be selected for pairing.
#' @param parents A list of selected parents (indices or individuals) from the
#'   population to be used for pairing and offspring generation.
#'
#' @details The function creates pairs of parents from `parents` and uses the
#' `crosser` function defined in `clf$functions` to generate offspring for each
#' pair. If there aren’t enough unique parents, it reuses the full list to
#' ensure the required number of couples. Offspring generation can be
#' parallelized if `clf$params$parallel.local` is set to `TRUE`.
#'
#' The process includes:
#' - **Pairing**: Randomly pairs parents from `parents`.
#' - **Offspring Generation**: Calls the `crosser` function to create offspring for each pair.
#' - **Parallelization**: Uses parallel computation if enabled, otherwise processes sequentially.
#'
#' @return A combined list of the original parents and their offspring.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10) # Random features
#' y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
#' clf <- list(params = list(parallel.local = FALSE),
#'             functions = list(crosser = crossingIndividual_v3)) # Assume `crossingIndividual_v3` is defined
#' pop <- list(
#'   list(indices_ = c(1, 3, 5), fit_ = 0.8),
#'   list(indices_ = c(2, 4, 6), fit_ = 0.7)
#' )
#' parents <- 1:2
#' combined_population <- crossing2(X, y, clf, pop, parents)
#' print(combined_population)
#' }
#'
#' @export
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

#' Select Top Individuals from a Population Based on Fitness
#'
#' This function selects a specified number of top individuals from a population
#' based on their fitness scores. It sorts the population in descending order of
#' fitness and returns the highest-ranking individuals.
#'
#' @param pop A list representing the population, where each individual has a
#'   `fit_` attribute that represents its fitness score.
#' @param number An integer specifying the number of individuals to select from
#'   the sorted population.
#' @param clf A classifier object containing any additional parameters needed
#'   for selection (though not directly used in this function).
#'
#' @details The function sorts the population in descending order by the `fit_`
#' attribute and returns the top `number` individuals. This selection strategy
#' prioritizes individuals with the highest fitness scores, making it useful for
#' selecting elites in evolutionary algorithms.
#'
#' @return A list of the top `number` individuals from the population, sorted by
#'   fitness.
#'
#' @examples
#' \dontrun{
#' pop <- list(
#'   list(fit_ = 0.9, indices_ = 1:5),
#'   list(fit_ = 0.8, indices_ = 6:10),
#'   list(fit_ = 0.7, indices_ = 11:15)
#' )
#' clf <- list() # Placeholder for classifier settings
#' selected_individuals <- selector_v1(pop, number = 2, clf)
#' print(selected_individuals)
#' }
#'
#' @export
selector_v1 <- function(pop, number, clf)
{
  # evaluation <- getFitPopulation(population)
  # ind <- order(evaluation,decreasing = T)[1:number]
  return(sortPopulation(pop, evalToOrder = "fit_")[1:number])
  # return(population[ind])
}



#' Randomly Select Individuals from a Population
#'
#' This function randomly selects a specified number of individuals from a
#' population without considering their fitness scores. It is a simple selection
#' strategy that introduces diversity by choosing individuals at random.
#'
#' @param pop A list representing the population from which individuals are
#'   selected.
#' @param number An integer specifying the number of individuals to select
#'   randomly from the population.
#' @param clf An optional classifier object that can contain additional
#'   parameters for selection (not directly used in this function).
#'
#' @details The function selects `number` individuals randomly from the
#' population `pop`, without replacement. This approach does not prioritize
#' fitness and can be useful for introducing genetic diversity in evolutionary
#' algorithms.
#'
#' @return A list of randomly selected individuals from the population.
#'
#' @examples
#' \dontrun{
#' pop <- list(
#'   list(fit_ = 0.9, indices_ = 1:5),
#'   list(fit_ = 0.8, indices_ = 6:10),
#'   list(fit_ = 0.7, indices_ = 11:15)
#' )
#' selected_individuals <- selector_v2(pop, number = 2)
#' print(selected_individuals)
#' }
#'
#' @export
selector_v2 <- function(pop, number, clf = NULL)
{ #TODO: changer ce selecteur
  ind <- sample(x = 1:length(pop), size  =  number, replace = FALSE)
  return(pop[ind])
}


#' Select Individuals from a Population Using Multiple Selection Methods
#'
#' This function selects individuals from a population using multiple selection
#' methods specified in the classifier. It divides the selection based on
#' specified percentages for each method, allowing a combination of selection
#' strategies.
#'
#' @param clf A classifier object containing parameters and functions for
#'   selection, including:
#' - `select_perc`: The percentage of the population to be selected.
#' - `size_pop`: The total size of the population.
#' - `select_percByMethod`: A vector of percentages defining the proportion of individuals to select by each method.
#' - `functions$selector`: A list of selection functions to apply, corresponding to each method in `select_percByMethod`.
#' @param pop A list representing the population from which individuals are
#'   selected.
#'
#' @details The function calculates the number of individuals to select
#' (`nb2BeMutated`) based on `select_perc` and the population size. It then
#' splits the selection across methods, as defined in `select_percByMethod`. For
#' each method, the respective selector function (from `clf$functions$selector`)
#' is applied to a subset of the population, excluding individuals already
#' selected by previous methods.
#'
#' The result is a list of selected individuals, combining the results of each
#' selection method.
#'
#' @return A list of selected individuals from the population, chosen by
#'   multiple selection methods.
#'
#' @examples
#' \dontrun{
#' pop <- list(
#'   list(fit_ = 0.9, indices_ = 1:5),
#'   list(fit_ = 0.8, indices_ = 6:10),
#'   list(fit_ = 0.7, indices_ = 11:15)
#' )
#' clf <- list(
#'   params = list(select_perc = 50, size_pop = 10, select_percByMethod = c(60, 40)),
#'   functions = list(selector = list(selector_v1, selector_v2))
#' )
#' selected_individuals <- select(clf, pop)
#' print(selected_individuals)
#' }
#'
#' @export
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


#' Select the Best Individual by Sparsity Level from a Population
#'
#' This function selects the best individual from a population for each
#' specified sparsity level, based on fitness. It groups individuals by their
#' sparsity, then identifies the highest-ranking individual in each group.
#'
#' @param clf A classifier object containing parameters for selection, including
#'   `sparsity`, which defines the target sparsity levels.
#' @param pop A list representing the population, where each individual has
#'   attributes like `eval.sparsity` (sparsity level) and `fit_` (fitness
#'   score).
#'
#' @details The function first groups individuals in `pop` by their sparsity
#' levels, as specified in `clf$params$sparsity`. For each group, it selects the
#' individual with the highest fitness score (`fit_`). Only sparsity levels with
#' at least one individual in the population are considered.
#'
#' This approach allows for maintaining a diverse set of individuals across
#' different sparsity levels, focusing on the best-performing individual within
#' each sparsity category.
#'
#' @return A list of the best individuals for each sparsity level specified in
#'   `clf$params$sparsity`, selected based on fitness.
#'
#' @examples
#' \dontrun{
#' pop <- list(
#'   list(eval.sparsity = 3, fit_ = 0.9),
#'   list(eval.sparsity = 4, fit_ = 0.8),
#'   list(eval.sparsity = 3, fit_ = 0.85)
#' )
#' clf <- list(params = list(sparsity = c(3, 4, 5)))
#' best_individuals <- selectIndividualToKeep(clf, pop)
#' print(best_individuals)
#' }
#'
#' @export
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

#' Convert an Individual's Sparse Gene Representation to a Dense Vector
#'
#' This function converts an individual's sparse gene representation (where only
#' specific indices have non-zero values) into a dense vector of a specified
#' size. It creates a vector with values from the individual's coefficients at
#' the specified indices and zeros elsewhere.
#'
#' @param individual A list representing an individual, with attributes:
#'   - `indices_`: A vector of indices where the individual has non-zero gene values.
#'   - `coeffs_`: A vector of coefficients corresponding to the genes at `indices_`.
#' @param size An integer specifying the size of the dense vector to be
#'   returned.
#'
#' @details The function initializes a dense vector (`res`) of length `size`
#' with all elements set to zero. It then populates this vector with the values
#' from `individual$coeffs_` at the positions specified in
#' `individual$indices_`. The result is a dense vector representation of the
#' individual's gene information.
#'
#' @return A numeric vector of length `size`, with non-zero values at positions
#'   specified by `individual$indices_` and zeroes elsewhere.
#'
#' @examples
#' \dontrun{
#' individual <- list(indices_ = c(2, 4, 6), coeffs_ = c(1.5, -0.5, 2.0))
#' size <- 10
#' dense_vector <- individualGetDenseVec(individual, size)
#' print(dense_vector) # Should show a vector of length 10 with specified values at indices 2, 4, and 6
#' }
#'
#' @export
individualGetDenseVec <- function(individual, size)
{
  res <- rep(0, size)
  res[individual$indices_] <- individual$coeffs_
  return(res)
}

