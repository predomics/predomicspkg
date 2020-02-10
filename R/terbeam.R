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


#' terbeam: ternary beam searching algorithm
#'
#' @title terbeam
#' @description terbeam is a model search algorithm on a beam search approach.
#' @param sparsity: number of features in a given model. This is a vector with multiple lengths.
#' @param maxNbOfModels: number of models to be explored for a given k_sparsity. This is equivalent to a population size in terga.
#' @param nbVeryBest: is the number of features to be kept that appear in the very best models. They will be kept even if they are not frequent in the best models (default: 1 percent of maxNbOfModels).
#' @param nbBest: is the number of features that will be used to build the k+1 sparsity combinations (default: 10 percent of maxNbOfModels).
#' @param final.pop.perc: a percentage of nbVeryBest translates in a number of models to be kept for k_sparsity.
#' @param popSaveFile: (??)
#' @param saveFiles: ??
#' @param language is the language that is used by the different algorithms {bin, bininter, ter, terinter, ratio}, (default:"terinter")
#' @param scoreFormula: a Function that contains the ratio Formula or other specific ones
#' @param epsilon: a small value to be used with the ratio language (useCustomLanguage) (default: NULL). When null it is going to be calculated by the minimum value of X divided by 10.
#' @param objective: this can be auc, cor or aic. Terga can also predict regression, other than class prediction. (default:auc)
#' @param max.nb.features: focuses only on the subset of top most significant features (default:1000)
#' @param estimate_coefs: non ternary solution for the aic objective (default:FALSE)
#' @param evalToFit: The model performance attribute to use as fitting score (default:"fit_"). Other choices are c("auc_","accuracy_","precision_","recall_","f_score_")
#' @param k_penalty: Penalization of the fit by the k_sparsity (default: 0)
#' @param intercept: (??) (default:NULL)
#' @param testAllSigns: ??
#' @param plot: Plot different graphics (default:FALSE).
#' @param verbose: print out information on the progress of the algorithm (default:TRUE)
#' @param warnings: Print out warnings when runnig (default:FALSE).
#' @param debug: print debug information (default:FALSE)
#' @param print_ind_method: One of c("short","graphical") indicates how to print a model and subsequently a population during the run (default:"short").
#' @param nCores: the number of cores to execute the program. If nCores=1 than the program runs in a non parallel mode
#' @param parallelize.folds: parallelize folds when cross-validating (default:TRUE)
#' @param seed: the seed to be used for reproductibility. If seed=NULL than it is not taken into account (default:NULL).
#### TODO check
#' @param experiment.id: The id of the experiment that is to be used in the plots and comparitive analyses (default is the learner's name, when not specified)
#' @param experiment.description: A longer description of the experiment. This is important when many experiments are run and can also be printed in by the printExperiment function.
#' @param experiment.save: Data from an experiment can be saved with different levels of completness, with options to be selected from c("nothing", "minimal", "full"), default is "minimal"
#' @param parallel: parallel
#' @return an object containing a list of parameters for this classifier
#' @export
terBeam <- function(sparsity = 1:5, max.nb.features = 1000, 
                    # maxNbOfModels corresponds to the width of the beam search in terms of models
                    maxNbOfModels = 10000, # population size
                    # nbBest is the number of models to keep such that the most frequent feature in the best models are kept
                    nbBest = round(maxNbOfModels/10), 
                    # nbVeryBest 
                    nbVeryBest = round(maxNbOfModels/100), 
                    final.pop.perc = 100,  # percentage of models in the returned results
                    # population options
                    popSaveFile = "NULL", saveFiles = FALSE, 
                    # language in {bin, bininter, ter, terinter, ratio}
                    language = "terinter", 
                    # language options
                    scoreFormula=scoreRatio, epsilon = "NULL",
                    # evaluation options
                    objective = "auc", k_penalty=0, evalToFit = 'auc_', estimate_coefs=FALSE, intercept = "NULL", testAllSigns = FALSE, 
                    # output options
                    plot = FALSE, verbose = TRUE, warnings = FALSE, debug = FALSE, print_ind_method = "short", parallelize.folds = TRUE,
                    # computing options
                    nCores = 4, seed = "NULL", #maxTime = Inf,
                    # experiment options
                    experiment.id = "NULL", experiment.description = "NULL", experiment.save = "nothing")
{
  # standard means that we use the standard heuristics
  clf <- list()
  clf$learner <- "terBeam"
  clf$params  <- list()
  clf$experiment <- list()                        # information about the experiment
  clf$params$objective          <- objective
  clf$params$estimate_coefs     <- estimate_coefs
  clf$params$sparsity           <- sparsity
  clf$params$max.nb.features    <- max.nb.features
  #  clf$params$maxBeam           <- maxBeam 
  #  clf$params$FILENAME          <- FILENAME
  #  clf$params$PREFIX            <- PREFIX    
  clf$params$saveFiles          <- saveFiles      # It would be interesting to add this in the future
  #  clf$params$pathSave          <- pathSave    
  
  #  clf$params$size_pop          <- size_pop
  clf$params$maxNbOfModels      <- maxNbOfModels
  clf$params$nbBest             <- nbBest
  clf$params$nbVeryBest         <- nbVeryBest
  
  # print out intermediary results
  clf$params$plot               <- plot           # print out logs.
  clf$params$verbose            <- verbose        # print out logs.
  clf$params$warnings           <- warnings       # print out warnings
  clf$params$debug              <- debug          # print out debugging information.
  clf$params$print_ind_method   <- print_ind_method # method to print individual
  
  # Computing options
  clf$params$nCores             <- nCores         # parallel computing
  clf$params$parallel           <- nCores > 1     # parallel computing
  clf$params$parallelize.folds  <- parallelize.folds
  clf$params$parallel.local     <- FALSE
  clf$params$seed               <- seed           # fix the seed to be able to reproduce results
  
  clf$params$testAllSigns       <- testAllSigns
  clf$params$final.pop.perc     <- final.pop.perc
  clf$params$evalToFit          <- evalToFit
  clf$params$k_penalty          <- k_penalty
  
  clf$params$language           <- language
  clf$params$popSaveFile        <- popSaveFile
  clf$params$epsilon            <- epsilon
  clf$params$scoreFormula       <- scoreFormula
  clf$params$intercept          <- intercept
  
  # Experiment information
  if(!(experiment.id=="NULL"))
  {
    clf$experiment$id           <- experiment.id
  }else
  {
    clf$experiment$id           <- clf$learner
  }
  
  if(!(experiment.description=="NULL"))
  {
    clf$experiment$description  <- experiment.description
  } else 
  {
    clf$experiment$description  <- paste(clf$learner, date() , sep = " ")
  }
  #match.arg(experiment.save) 
  clf$experiment$save           <- experiment.save
  
  return(clf)
}

terBeam_fit <- function(X, y, clf)
{
  
  # Setting the language environment
  switch(clf$params$language, 
         ter= 
         {
           # ternary language without intercept (maximize the auc)
           if(clf$params$verbose){print("Setting environment for the language 'ter'")}
           if(clf$params$objective == "cor")
           {
             clf$params$evalToFit <- "cor_"
           }else
           {
             # note that here with the ter language we could not use auc to fit since the intercept should be 0
             clf$params$intercept = 0
             if(clf$params$evalToFit == "auc_")
             {
               clf$params$evalToFit <- "accuracy_"
               warning("terga1_fit: changing evalToFit from auc_ to accuracy_ because of the language.")
             }
           }
         },
         terinter=
         {
           # ternary language without intercept (maximize the accuracy)
           if(clf$params$verbose){print("Setting environment for the language 'terinter'")}
           if(clf$params$objective == "cor")
           {
             clf$params$evalToFit <- "cor_"
           }
         },
         bin=
         {
           # ternary language without intercept (maximize the auc)
           if(clf$params$verbose){print("Setting environment for the language 'bin'")}
           if(clf$params$objective == "cor")
           {
             clf$params$evalToFit <- "cor_"
           }else
           {
             # note that here with the ter language we could not use auc to fit since the intercept should be 0
             clf$params$intercept = 0
             if(clf$params$evalToFit == "auc_")
             {
               clf$params$evalToFit <- "accuracy_"
               warning("terga1_fit: changing evalToFit from auc_ to accuracy_ because of the language.")
             }
           }
         },
         bininter=
         {
           # ternary language without intercept (maximize the auc)
           if(clf$params$verbose){print("Setting environment for the language 'bininter'")}
           if(clf$params$objective == "cor")
           {
             clf$params$evalToFit <- "cor_"
           }
         },
         ratio=
         {
           # ternary language without intercept (maximize the auc)
           if(clf$params$verbose){print("Setting environment for the language 'ratio'")}
           if(clf$params$objective == "cor")
           {
             clf$params$evalToFit <- "cor_"
           }
         },
         {
           stop(paste("The language",clf$params$language, "is not implemented !"))
         }
  )
  
  if(clf$params$verbose) print(paste("... ... parameters are checked and set"))
  # Print the experiment configuration
  if(clf$params$verbose) printClassifier(obj = clf)
  # Rprof("Profiling_terbeam", line.profiling = TRUE)
  
  # Rprof(NULL)
  # summaryRprof("Profiling_terbeam", lines = "show")
  # store all the features to keep
  features.to.keep <- allFeatures <- rownames(X)
  fullPop <- list()
  res.mod.coll <- list()
  
  # for each sparsity
  for(k in clf$params$sparsity)
  {
    if(k == 1)
    {
      # For k = 1 we generate every possible Model with only one feature
      if(!is.null(clf$feature.cor))
      {
        # for the ration language, force to have the same number of negative and positive features selected
        if((clf$params$language == "ratio" | clf$params$language == "ter" | clf$params$language == "terinter") & clf$params$objective == "auc")
        {
          # select the best features here no need to compute all
          nb.selected.features <- min(nrow(clf$feature.cor),clf$params$maxNbOfModels)
          nb.selected.features.neg <- nb.selected.features.pos <- round(nb.selected.features/2)
          # get the pool of features
          features.pool <- as.character(rownames(clf$feature.cor))[order(clf$feature.cor$p)]
          features.pool.coeffs <- clf$coeffs_[features.pool]
          # negative
          neg.ind <- features.pool.coeffs == "-1" & !is.na(features.pool.coeffs)
          selected.features.neg <- features.pool[neg.ind][1:min(sum(neg.ind), nb.selected.features.neg)]
          # positive
          pos.ind <- features.pool.coeffs == "1" & !is.na(features.pool.coeffs)
          selected.features.pos <- features.pool[pos.ind][1:min(sum(pos.ind), nb.selected.features.neg)]
          
          selected.features <- c(selected.features.neg, selected.features.pos)
        }else
        {
          nb.selected.features <- min(nrow(clf$feature.cor), clf$params$maxNbOfModels)
          # get the pool of features
          features.pool <- rownames(clf$feature.cor)[order(clf$feature.cor$p)]
          # select the best features here no need to compute all
          selected.features <- features.pool[1:nb.selected.features]
        }
        # get the index in the rownames
        ind.features    <- which(rownames(X) %in% selected.features)
        if(clf$params$verbose) print(paste("... ... generating only best single feature models"))
      }else
      {
        if(clf$params$verbose) print(paste("... ... generating all single feature models"))
        stop("terBeam_fit: clf$feature.cor is missing")
        ind.features    <- seq_len(nrow(X))
      }
      
      pop               <- generateAllSingleFeatureModel(X, y, clf, ind.sub = ind.features)
      
    } else
    {
      ### Get the features to keep for next k
      nbCombinaisons    <- choose(n = c(1:length(features.to.keep)), k = k)
      
      
      # for the ration language, force to have the same number of negative and positive features selected
      if((clf$params$language == "ratio" | clf$params$language == "ter" | clf$params$language == "terinter") & clf$params$objective == "auc")
      {
        # select the best features here no need to compute all
        nb.selected.features <- max(which(nbCombinaisons < clf$params$maxNbOfModels))
        nb.selected.features.neg <- nb.selected.features.pos <- round(nb.selected.features/2)
        # get the pool of features
        features.pool <- intersect(as.character(rownames(clf$feature.cor))[order(clf$feature.cor$p)], features.to.keep)
        features.pool.coeffs <- clf$coeffs_[features.pool]
        # negative
        neg.ind <- features.pool.coeffs == "-1" & !is.na(features.pool.coeffs)
        selected.features.neg <- features.pool[neg.ind][1:min(sum(neg.ind), nb.selected.features.neg)]
        # positive
        pos.ind <- features.pool.coeffs == "1" & !is.na(features.pool.coeffs)
        selected.features.pos <- features.pool[pos.ind][1:min(sum(pos.ind), nb.selected.features.neg)]
        features.to.keep <- selected.features <- c(selected.features.neg, selected.features.pos)
        
      }else
      {
        nb.selected.features <- max(which(nbCombinaisons < clf$params$maxNbOfModels))
        # get the pool of features
        features.pool <- intersect(as.character(rownames(clf$feature.cor))[order(clf$feature.cor$p)], features.to.keep)
        
        # select the best features here no need to compute all
        features.to.keep <- selected.features <- features.pool[1:min(nb.selected.features, length(features.pool))]
      }
      
      # # nbCombinaisons contains the maximum number of models that would be generated for a given number of features
      # features.to.keep  <- features.to.keep[1:max(which(nbCombinaisons < clf$params$maxNbOfModels))] 
      
      ### For k > 1 we generate every possible combinations of features of size k
      ind.features.to.keep <- which(allFeatures %in% features.to.keep)
      
      if(length(ind.features.to.keep) >= k)
      {
        pop               <- generateAllCombinations(X = X, y = y, clf = clf, 
                                                     ind.features.to.keep = ind.features.to.keep, 
                                                     sparsity = k, 
                                                     allFeatures = allFeatures)
      }else
      {
        break
      }
      
    }
    
    # Evaluate the population
    pop                 <- evaluatePopulation(X = X, y = y, clf = clf, pop = pop, 
                                              eval.all = TRUE, 
                                              force.re.evaluation = TRUE, 
                                              estim.feat.importance = FALSE,
                                              mode = "train",
                                              delete.null.models = TRUE)
    
    # Sort the population according to the clf$params$evalToFit attribute
    pop                 <- sortPopulation(pop, evalToOrder = "fit_")
    
    # it may happen that the population is empty, in this case we break and return
    if(is.null(pop))
    {
      next
    }
    
    # Sample the best and veryBest population
    best                <- pop[1:min(clf$params$nbBest,length(pop))]
    veryBest            <- pop[1:min(clf$params$nbVeryBest,length(pop))]
    
    # features.to.keep
    ### Evaluate the apperance of every features in the best and veryBest models
    featuresApperance   <- countEachFeatureApperance(clf, allFeatures, pop, best, veryBest)
    features.to.keep    <- getFeatures2Keep(clf, featuresApperance)
    
    if(clf$params$verbose)
    {
      if(isModel(pop[[1]]))
      {
        try(printModel(mod = pop[[1]], method = clf$params$print_ind_method, score = "fit_"), silent = TRUE)
      }
    }
    
    #EP: keep only the verybest
    #fullPop[(length(fullPop) +1):(length(fullPop) + length(pop))] <- pop
    minsize <- min(length(pop), clf$params$nbVeryBest)
    fullPop[(length(fullPop) +1):(length(fullPop) + minsize)] <- pop[1:minsize]
    
    # save populatio in a file
    if(!(clf$params$popSaveFile=="NULL"))
    {
      savePopulation(fullPop, paste("resulstsForSparsity", k, clf$params$popSaveFile, sep = "_"))
    }
    
    # stopping testing
    if((length(features.to.keep) < k + 2) & (k != 1)) # If we exhausted all the combinations
    {
      break      
    } # end if stopping test
  } # end loop sparsity
  
  if(clf$params$verbose) print(paste("... ... models are created"))
  
  return.perc       <- clf$params$final.pop.perc
  if(return.perc > 100)
  {
    return.perc = 100 # upper bound
    warning("terBeam_fit: clf$params$final.pop.perc can not be greater than 100")
  }
  
  #fullPop           <- sortPopulation(fullPop, evalToOrder = "fit_")
  fullPop           <- unique(fullPop) # keep only unique models
  
  if(return.perc == 100)
  {
    # transform the population onto a model collection
    res.mod.coll    <- listOfModels2ModelCollection(pop = fullPop)
  } else # if smaller percentage
  {
    #if(clf$params$final.pop.perc>100)
    nBest           <- round(return.perc * clf$params$nbVeryBest / 100)
    res.mod.coll    <- listOfModels2ModelCollection(pop = fullPop, nBest = nBest)
  }
  
  if(clf$params$verbose) print(paste("... ... models are coverted onto a model collection"))
  return(res.mod.coll) 
}

