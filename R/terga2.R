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
# @script: terga2.R                                          
# @author: Edi Prifti
# @author: Lucas Robin
# @date: August 2016                                                    
################################################################

#' Model search algorithm based on genetic algorithms (GA).
#' @description TerGA is a model search algorithm based on genetic algorithms (GA). 
#' An “individual” (i.e. genome) in this context is a combination of features that 
#' will be associated together using a selected "language" to compute a score that 
#' will constitute the prediction model. Depending on the type of fitting (i.e. evaluation)
#' function that is maximized, the fatures are weighed by specific coefficients. 
#' In short the algorithm is based on different operations such as crossing, mutating 
#' and evolving different “individuals” and evaluating their fitness to the “environment” 
#' which is represented by the variable to be predicted.
#' @param sparsity: number of features in a given model (default:1:10). 
#' This is a vector with the model-size range (number of features used by a model).
#' @param language is the language that is used by the different algorithms 
#' {bin, bininter, ter, terinter}, (default:"terinter")
#' @param objective: This is the task that is to be learned and can be either classification 
#' (auc) or can be a regression (cor) (default:auc).
#' @param evalToFit: The model performance attribute to use as fitting score (default:"accuracy_"). 
#' Other choices are c("accuracy_", "auc_", "precision_","recall_","f_score_") for the 
#' classification task. It can be either rho, rho-squared or minimizing the 
#' standar error of the regression for the regression task.
#' @param k_penalty: Model-size penalization effect applied on the fit scpre (default: 0).
#' @param estimate_coefs: _deprecated_ A particular option for the regression mode
#'  with the aic objective (default:FALSE)
#' @param max.nb.features: If this number is smaller than the number of variables in the
#' dataset, the max.nb.features most significant features will be selected and the 
#' dataset will be restricted (default:1000).
#' @param size_pop: the number of individuals in a population to be evolved (default:100)
#' @param size_pop_random the number of individuals initialized randomly. This is used 
#' by the metal algorithm (i.e. aggregator method).
#' @param final.pop.perc: What percentage of the final population should be returned (default:100)
#' @param in_pop: a specific population of models that can be evolved. This is particulary
#' useful for the metal algorithm
#' @param popSourceFile: It is possible to load a population of models that has been
#' already learned before. With this option we can specify such file (default:NULL).
#' @param popSaveFile: Once the population of models evolved, we can store it in 
#' another file (default:NULL).
#' @param scoreFormula: a Function that contains the ratio Formula or other specific ones
#' @param epsilon: a very small value to be used with the ratio language 
#' (useCustomLanguage) (default: NULL). When null it is going to be calculated by the 
#' minimum value of X divided by 10.
#' @param individual_vec: The function that is used to generate an individual 
#' (default:individual_vec_v2).
#' @param randomSigns: When generating an individual composed of a set of features, we 
#' can set the coefficients of the variables from -1 or 1 randomly (default:FALSE).
#' @param unique_vars: When performing operations on multiple individuals it can be 
#' that in an individual we have multiple time the same feature. If set to TRUE this 
#' individual will be destroyed (default:FALSE)
#' @param select_perc: The percentage of the population to be selected for crossing/mutation 
#' (default:50)
#' @param selector: During the selection process, the parent population can be
#' selected using different strategies. For instance the default process is performed
#' using both elite and random selection (default:list(selector_v1, selector_v2)).
#' @param select_percByMethod: A list contaning the percentage of individuals that
#' each of the methods specified in selector should get.
#' @param cross: A swithch, which activates the crossing operator (default:TRUE).
#' @param crosser: The method that should be applied to cross individuals
#' together (default:crossingIndividual_v4).
#' @param mutate: A swithch, which activates the mutation operator (default:TRUE).
#' @param mutate_size: The percentage of individuals in the population to be mutated (default:70).
#' @param mutate_rate: The percentage of features in an individual to be mutated (default:50).
#' @param mutator: The method that should be applied to mutate individuals (default:mutator_v2).
#' The operations can be, deletion, insertion or changing the coeffiecient (from -1 to 1
#' and vice-versa).
#' @param evolver: The method that will be used to evolve the individuals together.
#' This is the core of the algorithm and can be one of different implementations 
#' c("v1", "v2", "v3","v4") where the default one is "v4".
#' @param nb_generations: The maximum number of generations to evolve the population.
#' @param convergence: A switch which activates the automatic convergence of the algorithm
#' when the best individual is not improving (default:TRUE).
#' @param convergence_steps: The number of generations after which we consider 
#' convergence (default:10).
#' @param evolve_k1: Whether or not to evaluate exhaustively the features for 
#' model size = 1. This will take a lot of time if the dataset is large, thus the
#' possibility to evolve this using the GA is interesting. (default:TRUE)
#' @param plot: Plot graphics indicating the evolution of the simulation (default:FALSE)
#' @param verbose: Print out information on the progress of the algorithm (default:FALSE).
#' @param warnings: Print out warnings when runnig (default:FALSE).
#' @param debug: Print out detailed information on the progress of the algorithm 
#' (default:FALSE)
#' @param print_ind_method: One of c("short","graphical") indicates how to print 
#' a model and subsequently a population during the run (default:"short").
#' @param parallelize.folds: parallelize folds when cross-validating (default:TRUE).
#' @param nCores: The number of cores to execute the program. If nCores = 1 than 
#' the program runs in a non parallel mode
#' @param seed: The seed to be used for reproductibility. If seed=NULL than it is 
#' not taken into account (default:NULL).
#' @param maxTime: We can use a time limit to evolve a population (default:Inf).
#' @param experiment.id: The id of the experiment that is to be used in the plots 
#' and comparitive analyses (default is the learner's name, when not specified)
#' @param experiment.description: A longer description of the experiment. This is 
#' important when many experiments are run and can also be printed in by the 
#' printExperiment function.
#' @param experiment.save: Data from an experiment can be saved with different 
#' levels of completness, with options to be selected from 
#' c("nothing", "minimal", "full"), default is "minimal"
#' @return an object of the classifier class, containing a list of parameters
#' @export
terga2 <- function(sparsity = c(1:10), max.nb.features = 1000, 
                   # language in {bin, bininter, ter, terinter, ratio}
                   language = "terinter",
                   # evaluation options
                   objective = "auc", evalToFit = "accuracy_", k_penalty = 0, estimate_coefs = FALSE, 
                   # ratio formula (or other future specific)
                   scoreFormula = scoreRatio, epsilon = "NULL",
                   # population options
                   size_pop = 100, size_pop_random = size_pop, final.pop.perc = 100, 
                   in_pop = "NULL", popSourceFile = "NULL", popSaveFile = "NULL", 
                   # individual options
                   individual_vec = individual_vec_v2, randomSigns = FALSE, unique_vars = FALSE,
                   # selection options
                   select_perc = 25, selector = list(selector_v1, selector_v2), select_percByMethod = list(50, 50),
                   # crossing options
                   cross = TRUE, crosser = crossingIndividual_v3,
                   # mutation options
                   mutate = TRUE,  mutate_size = 75, mutate_rate = 50, mutator = mutator_v2, 
                   # evoluion options
                   evolver = "v2m", nb_generations = 100, convergence = TRUE, convergence_steps = 10, evolve_k1 = TRUE,
                   # output options
                   plot = FALSE, verbose = FALSE, warnings = FALSE, debug = FALSE, print_ind_method = "short", parallelize.folds = TRUE,
                   # computing options
                   nCores = 4, seed = "NULL", maxTime = Inf,
                   # experiment options
                   experiment.id = "NULL", experiment.description = "NULL", experiment.save = "nothing")
{
  # some paramters need to be added to the description

  clf                           <- list()         # create a classifier object
  clf$learner                   <- "terga2"       # name of the method
  clf$params                    <- list()         # parameter list
  clf$experiment                <- list()         # information about the experiment
  
  # POPULATION
  clf$params$sparsity           <- sparsity       # number of non zero variables in the model
  clf$params$sparsity.min       <- min(sparsity)
  clf$params$sparsity.max       <- max(sparsity)
  clf$params$sparsity.mean      <- mean(sparsity)
  clf$params$current_sparsity   <- NA             # number of non zero variables in the model
  
  # Feature space information
  clf$params$max.nb.features    <- max.nb.features
  clf$params$size_world         <- "NULL"         # total number of variables
  clf$params$unique_vars        <- unique_vars    # weather in a model we can have one variable more than once
  clf$params$testAllSigns       <- FALSE          # _deprecated_
  clf$params$randomSigns        <- randomSigns
  
  # FITTING
  clf$params$objective          <- objective      # prediction (roc) or regression (cor)
  clf$params$estimate_coefs     <- estimate_coefs # integer or real estimated coefficients.
  clf$params$language           <- language
  clf$params$evalToFit          <- evalToFit
  clf$params$k_penalty          <- k_penalty
  clf$params$epsilon            <- epsilon  
  clf$params$scoreFormula       <- scoreFormula
  clf$params$intercept          <- "NULL"
  
  # POPULATION
  clf$functions$individual_vec  <- individual_vec # The function that will be used to make individuals
  
  clf$params$in_pop             <- in_pop
  clf$params$popSourceFile      <- popSourceFile
  clf$params$popSaveFile        <- popSaveFile
  clf$params$size_pop           <- size_pop       # how many models in the population to evolve
  clf$params$size_pop_random    <- size_pop_random
  clf$params$final.pop.perc     <- final.pop.perc
  
  # SELECTION
  clf$functions$selector        <- selector      # il faudra vérifier que cette liste est de la même longueur que select_ByMethod
  clf$params$select_perc        <- select_perc   # the percentage of the population to be selected
  clf$params$select_percByMethod<- select_percByMethod # how much each selector should select
  if(length(selector) != length(clf$params$select_percByMethod))
  {
    stop("terga2: please verify that select_percByMethod is a list indicating a percentage for each selector function provided")
  }
  
  if(sum(unlist(clf$params$select_percByMethod)) != 100)
  {
    stop("terga2: please verify that select_percByMethod is a list containing percentages whose sum should be 100 percent")
  }
    
  # CROSSING
  clf$params$cross              <- cross          # do we activate the crossing process
  clf$functions$crosser         <- crosser        # The function that will be applied to cross
  
  # MUTATION
  clf$params$mutate             <- mutate         # do we activate the mutation process
  clf$params$mutate_size        <- mutate_size    # what percentage of models in the population are mutated
  clf$params$mutate_rate        <- mutate_rate    # what percentage of the variables in the model are mutated
  clf$functions$mutator         <- mutator        # The function that will be applied
  
  # CONVERGENCE
  clf$params$evolver            <- evolver        # The method be used during the evolution process
# TODO make this a function as for the rest
  clf$params$nb_generations     <- nb_generations # number of generation to evolve
  clf$params$convergence        <- convergence    # what should the simulation stop when convergence ?
# TODO this needs to be tested
  clf$params$convergence_steps  <- convergence_steps # after how many steps without improvement do we consider convergence?
  clf$params$evolve_k1          <- evolve_k1      # weather to evolve models with k_1 or to search them exhaustively.
  
  # print out intermediary results
  clf$params$plot               <- plot           # plot results? 
  clf$params$verbose            <- verbose        # print out logs.
  clf$params$warnings           <- warnings       # print out warnings
  clf$params$debug              <- debug          # print out logs.
  clf$params$print_ind_method   <- print_ind_method # method to print individual
  
  # Computing options
  clf$params$nCores             <- nCores         # parallel computing
  clf$params$parallel           <- nCores > 1     # parallel computing
  clf$params$parallelize_folds  <- parallelize.folds
  clf$params$parallel.local     <- FALSE
  clf$params$seed               <- seed           # fix the seed to be able to reproduce results
  clf$params$maxTime            <- maxTime
  
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
    clf$experiment$description <- experiment.description
  } else 
  {
    clf$experiment$description <- paste(clf$learner, date() , sep = " ")
  }
  #match.arg(experiment.save)
  clf$experiment$save          <- experiment.save
  
  # set the class
  class(clf) <- c("classifier", "predomics")
  
  return(clf)
}



#######################################################
##      New Version of the terga_fit function :      ##
#######################################################
terga2_fit <- function(X, y, clf)
{
  startingTime <- Sys.time()
  
  # NOTE: for speed and simplicty purposes we decided to use accuracy_ as the 
  # variable to fit the models. We noticed that when you optimize with auc_ the 
  # best models are not necessarely the best in accuracy. We hypothesize that 
  # when optimizing in accuracy the models will be rather good also in auc_ 
  # based on how the intercept evaluation is performed is very similar to AUC.
  # For ter and bin languages the threshold will be set to zero.
  
  switch(clf$params$language, 
         ter= 
         {
           # ternary language without intercept (maximize the accuracy)
           if(clf$params$verbose){cat("... ... Setting environment for the language 'ter'\n")}
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
               if(clf$params$warnings) warning("terga1_fit: changing evalToFit from auc_ to accuracy_ because of the language.")
             }
           }
         },
         terinter=
         {
           # ternary language without intercept (maximize the auc)
           if(clf$params$verbose){cat("... ... Setting environment for the language 'terinter'\n")}
           if(clf$params$objective == "cor")
           {
             clf$params$evalToFit <- "cor_"
           }
         },
         bin=
         {
           # ternary language without intercept (maximize the accuracy)
           if(clf$params$verbose){cat("... ... Setting environment for the language 'bin'\n")}
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
               if(clf$params$warnings) warning("terga1_fit: changing evalToFit from auc_ to accuracy_ because of the language.")
             }
           }
         },
         bininter=
         {
           # ternary language without intercept (maximize the auc)
           if(clf$params$verbose){cat("... ... Setting environment for the language 'bininter'\n")}
           if(clf$params$objective == "cor")
           {
             clf$params$evalToFit <- "cor_"
           }
         },
         ratio=
         {
           # ternary language without intercept (maximize the auc)
           if(clf$params$verbose){cat("... ... Setting environment for the language 'ratio'\n")}
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
  
  if(clf$params$max.nb.features < nrow(X))
  {
    if(clf$params$verbose) print(paste("... ... restricting X to the",clf$params$max.nb.features,"most significant variables"))  
    selected.features <- rownames(clf$feature.cor[order(clf$feature.cor$p),][1:min(clf$params$max.nb.features, nrow(X)),])
    X <- X[selected.features,]
    clf$coeffs_ <- clf$coeffs_[selected.features]
    clf$feature.cor <- clf$feature.cor[selected.features,]
  }
  
  # set the size of the world
  clf$params$size_world <- nrow(X)
  
  # Print the experiment configuration
  if(clf$params$verbose) printClassifier(obj = clf)
  
  featEval <-rep(NA, length(rownames(X)))
  names(featEval) <- rownames(X)
  
  # parallel switch
  if(clf$params$parallel) 
  {
    parfold <- TRUE
  }else 
  {
    parfold <- FALSE
  }

  if(isPopulation(clf$params$in_pop))
  {
    # starting population
    pop <- clf$params$in_pop
    if(clf$params$verbose) cat(paste("... ... using seeding population",length(pop),"\n"))
  }else
  {
    # otherwise start fresh
    pop <- population2(X, y, clf, featEval) # We generate randomly a first population
    if(clf$params$verbose) cat(paste("... ... using generated population",length(pop),"\n"))
  }
  
  if(clf$params$debug) 
  {
    print(paste("population:", length(pop)))
  }
    
  
  # it may happen that the population is empty, in this case return NULL
  if(is.null(pop))
  {
    return(NULL)
  }
  
  if(clf$params$debug)
  {
    v <- populationGet_X(element2get = "fit_", toVec = TRUE, na.rm = TRUE)(pop)
    summary(v)
  }
  featEval <- rankFeatures(X, y, clf, pop, featEval)

  stop <- FALSE # En utilisant une variable pour la condition d'arrêt on peut calculer l'arrêt du fit
                # via une fonction qui elle même pourra être un paramêtre donné dans le clf
  
  generation <- 0
  
  doWeStopHere <- function(X, y, clf, pop, generation, startingTime) #{ifelse(generation == clf$params$nb_generations, TRUE, FALSE)}
  {
    if(is.finite(clf$params$maxTime))
    {
      execTime <- Sys.time() - startingTime
      print(execTime)
      if(execTime > clf$params$maxTime)
        return(TRUE)
      return(FALSE)
    } else
    {
      ifelse(generation == clf$params$nb_generations, TRUE, FALSE)
    }
  }
  
  # initiate variables
  trace_evolution <- list()
  evaluation <- c()
  
  while(!stop)
  {
    if(clf$params$verbose) 
    {
      print(paste("Generation", generation))
    }else
    {
      # print a heartbeat
      cat(".")
    }
    
    pop <- cleanPopulation(pop, clf)
    
    if(clf$params$debug) 
    {
      print(paste("population after clean:", length(pop)))
    }
    
    pop <- evaluatePopulation(X, y, clf, pop, force.re.evaluation = TRUE, eval.all = FALSE) 
    
    # Let's keep the best models of each sparsity in a buffer structured as modelCollection
    pop2keep.mc <- listOfModels2ModelCollection(pop = pop, nBest = 5)
    
    # keep valid populations and omit those in the mc that are no ok
    valid.pops <- unlist(lapply(pop2keep.mc, isPopulation))
    if(any(!valid.pops))
    {
      pop2keep.mc <-pop2keep.mc[valid.pops]
    }
    
    # restrain if bigger than the given size but by keeping the best for each sparsity
    # sometimes the best models will be in high sparsity setting
    if(length(pop) > clf$params$size_pop)
    {
      nBest           <- round(length(pop)/length(clf$params$sparsity))
      pop.restrict    <- listOfModels2ModelCollection(pop = pop, nBest = nBest)
      pop             <- modelCollectionToPopulation(pop.restrict)
    }
    
    # print the best for each sparsity
    if(clf$params$verbose) 
    {
      print("Best individuals kept for the next generation :")
      lapply(pop2keep.mc, function(x)
      {
        if(!is.null(x))
        {
          try(cat("\t", printModel(mod = x[[1]], method = clf$params$print_ind_method, score = "fit_"), "\n"), silent = TRUE)
        } else 
        {
          cat("No individual found for this sparsity")
        }
      })
      spar <- populationGet_X(element2get = "eval.sparsity", toVec = TRUE, na.rm = TRUE)(pop)
      print(table(spar))
    }
    
    # Evolve methods
    switch (clf$params$evolver,
      v1 = {pop       <- evolve1(X, y, clf, pop, featEval)},
      v2m = {pop      <- evolve2m(X, y, clf, pop, featEval, generation)},
      v3m = {pop      <- evolve3m(X, y, clf, pop, featEval)},
      {pop            <- evolve2m(X, y, clf, pop, featEval)} # default is v2
    )
    
    if(clf$params$debug) 
    {
      print(paste("population after evolver:", length(pop)))
    }
    
    
    if(clf$params$evolver %in% list('v2m') && generation == 0)
    {
      clf$params$size_pop <- round(clf$params$size_pop*3/2)
    }
    
    # transform the bugger onto a population
    pop2keep <- modelCollectionToPopulation(pop2keep.mc)
    
    # and add it to the evolved pop
    #pop[(length(pop) +1):(length(pop) + length(pop2keep))] <- pop2keep
    pop <- c(pop, pop2keep)
    pop <- unique(pop) # delete dupplicates for efficiency
    # sort the population
    pop <- sortPopulation(pop, evalToOrder = "fit_")
    
    # DEBUG and trace evolution
    best <- pop[[1]]
    if(clf$params$debug)
    {
      printy(best)
    }
    trace_evolution[[generation+1]] <- best ## A utiliser pour tester la convergence
    evaluation <- populationGet_X(element2get = "fit_", toVec = TRUE, na.rm = TRUE)(pop = trace_evolution)
    
    
    if(clf$params$debug) 
    {
      print(paste("population end generation:", summary(evaluation)))
    }
    
    
    # # to be updated by the app
    # if(clf$params$plot)
    # {
    #   plot(evaluation, xlim=c(1,clf$params$nb_generations), type = "b", col = "red", ylim = c(0.5,1), cex = 0.5, pch = 19)
    #   abline(h = max(table(y)/sum(table(y))), lty=2)
    # }
    
    featEval <- rankFeatures(X, y, clf, pop, featEval)
    
#========> TODO add convergence early stopping and also print out evolution through epochs
    # # Test convergence and exit if atteint
    # if (clf$params$convergence)
    # {  # if we are testing convergence
    #   trace.nona <- trace_evolution[!is.na(trace_evolution)]
    #   if (length(trace.nona) > clf$params$convergence_steps)
    #   {  # run at least enough steps to test the convergence
    #     if (convergence.test(x = trace.nona, steps = clf$params$convergence_steps))
    #     {  # if convergence atteint
    #       break("")
    #     }
    #   }
    # }
#========> TODO     
 
    # #add again the best individuels in case that they are changed
    # if(!(length(clf$params$in_pop)==1 && clf$params$in_pop=="NULL"))
    # {
    #   pop_add <- clf$params$in_pop
    #   pop[(length(pop) +1):(length(pop) + length(pop_add))] <- pop_add
    #   pop <- unique(pop)
    # }
    
    # if(generation == 47)
    # {
    #   print("*********************************HERE*********************************")
    # }
    
    # CAUTION: we revaluate the population forcing the score after all the modifications of the indices during crossing, mutation etc...
    # clean the population
    #pop <- cleanPopulation(pop, clf)
    #pop <- evaluatePopulation(X, y, clf, pop, force.re.evaluation = TRUE, eval.all = TRUE) 
    
    # if(!(clf$params$popSaveFile=="NULL"))
    # {
    #   savePopulation(pop, paste("generation", generation, clf$params$popSaveFile, sep = "_"))
    # }
    
    # increase generation
    generation <- generation + 1
    stop <- doWeStopHere(X, y, clf, pop, generation, startingTime) # temporary fonction name, this could be latter something like
    # clf$params$stopCondition(X, y, clf, pop, generation)
    
  } # end while
  # end the heartbeat
  if(!clf$params$verbose) cat("\n")
  
  # to be updated by the app
  if(clf$params$plot)
  {
    plot(evaluation, xlim=c(1,clf$params$nb_generations), type = "b", col = "red", ylim = c(0.5,1), cex = 0.5, pch = 19)
    abline(h = max(table(y)/sum(table(y))), lty=2)
  }
  
  # We revaluate the population preparing for final result
  # clean the population
  pop <- cleanPopulation(pop, clf)
  pop <- evaluatePopulation(X, y, clf, pop, force.re.evaluation = TRUE, eval.all = TRUE)
  pop <- unique(pop) # keep only unique models
  
  # At the end we revaluate the population forcing the scores and everything else
  # pop <- evaluatePopulation(X, y, clf, pop, force.re.evaluation = TRUE, eval.all = FALSE) 
  if(clf$params$final.pop.perc == 100)
  {
    res <- listOfModels2ModelCollection(pop)
  } else 
  {
    res <- listOfModels2ModelCollection(pop, nBest = 10)
  }
  
  if(clf$params$verbose) print(paste("... ... models are coverted onto a model collection"))
  return(res)
}
