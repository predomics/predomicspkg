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
# @script: terga1.R                                          
# @author: Edi Prifti
# @author: Lucas Robin
# @date: August 2016                                                    
################################################################


#' terga1: Model search algorithm based on genetic algorithms (GA)
#'
#' @title terga1
#' @description terga1 is a model search algorithm based on genetic algorithms (GA). A “genome” or “individual” in this context is a combination of features that will be associated together to compute a score that will be the prediction model. Depending on the type of fitting function that is maximized the fatures are weighed by specific coefficients. In short the algorithm is based on different operations such as crossing, mutating and evolving different “individuals” and evaluating their fitness to the “environment” which is represented by the variable to be predicted.
#' @param sparsity: number of features in a given model. This is a vector with multiple lengths.
#' @param size_pop: the number of individuals in a population to be evolved.
#' @param size_world: this is the number of features in the dataset.
#' @param max.nb.features: focuses only on the subset of top most significant features (default:1000)
#' @param popSourceFile: A population of models that can start as a first generation to be evolved (default:NULL).
#' @param popSaveFile: (??)
#' @param language is the language that is used by the different algorithms {bin, bininter, ter, terinter, ratio}, (default:"terinter")
#' @param scoreFormula: a Function that contains the ratio Formula or other specific ones
#' @param epsilon: a small value to be used with the ratio language (default: NULL). When null it is going to be calculated by the minimum value of X divided by 10.
#' @param unique_vars: logical (default: FALSE) indicates weather unique variables can be used in a model or population.
#' @param objective: this can be auc, cor or aic. Terga can also predict regression, other than class prediction. (default:auc)
#### TODO UPDATE THIS AS c(classification, regression)
#' @param estimate_coefs: non ternary solution for the aic objective (default:FALSE)
#' @param intercept: (Interceot for the a given model) (default:NULL)
#' @param evalToFit: The model performance attribute to use as fitting score (default:"fit_"). Other choices are c("auc_","accuracy_","precision_","recall_","f_score_")
#' @param k_penalty: Penalization of the fit by the k_sparsity (default: 0)
#' @param select_type: the selection operator type. can be mixed, elite or tournoi (default: mixed)
#' @param select_perc1: percentage of individuals to be selected with elite
#' @param select_perc2: percentage of individuals to be selected with tournoi
#' @param perc_best_ancestor: percentage of best ancentors as seeding in the new population
#' @param mutate_size: percentage of individuals in the population to be mutated
#' @param mutate_rate: percentage of features in an individual to be mutated
#' @param plot: plot graphics indicating the evolution of the simulation (default:FALSE)
#### TOCHECK is this still needed for tergaV2
#' @param convergence: should the algorithm converge when the best individual is not improving (default:TRUE).
#' @param convergence_steps: the number of generations after which we consider convergence (default:10).
#' @param evolve_k1: weather or not to evaluate exhaustively the features for k_sparse=1. This will take a lot of time if the dataset is large, thus the possibility to evolve this using the GA. (default:TRUE)
#' @param verbose: print out information on the progress of the algorithm (default:TRUE)
#' @param warnings: Print out warnings when runnig (default:FALSE).
#' @param debug: print debug information (default:FALSE)
#' @param print_ind_method: One of c("short","graphical") indicates how to print a model and subsequently a population during the run (default:"short").
#' @param parallelize.folds: parallelize folds when cross-validating (default:TRUE)
#' @param nb_generations: maximum number of generations to evolve the population.
#' @param nCores: the number of cores to execute the program. If nCores=1 than the program runs in a non parallel mode
#' @param seed: the seed to be used for reproductibility. If seed=NULL than it is not taken into account (default:NULL).
#### TODO check
#' @param experiment.id: The id of the experiment that is to be used in the plots and comparitive analyses (default is the learner's name, when not specified)
#' @param experiment.description: A longer description of the experiment. This is important when many experiments are run and can also be printed in by the printExperiment function.
#' @param experiment.save: Data from an experiment can be saved with different levels of completness, with options to be selected from c("nothing", "minimal", "full"), default is "minimal"
#' @return an object containing a list of parameters for this classifier
#' @export
terga1 <- function(sparsity = c(1:10), 
                   # population options
                   size_pop = 100, size_world = "NULL", max.nb.features = 1000, popSourceFile = "NULL", popSaveFile = "NULL",
                   # language in {bin, bininter, ter, terinter, ratio}
                   language = "terinter",
                   # language options
                   scoreFormula=scoreRatio, epsilon = "NULL",
                   # individual options
                   unique_vars = FALSE,
                   # evaluation options
                   objective = "auc", k_penalty=0, evalToFit = "fit_", estimate_coefs = FALSE, intercept = "NULL",
                   # selection options
                   select_type = "mixed", select_perc1 = 20, select_perc2 = 30, perc_best_ancestor = 10,
                   # mutation options
                   mutate_size = 70, mutate_rate = 50,  
                   # evolution options
                   nb_generations = 100, convergence = TRUE, convergence_steps = 10, evolve_k1 = TRUE,
                   # output options
                   plot = FALSE, verbose = TRUE, warnings = FALSE, debug = FALSE, print_ind_method = "short", parallelize.folds = TRUE,
                   # computing options
                   nCores = 4, seed = "NULL", 
                   # experiment options
                   experiment.id = "NULL", experiment.description = "NULL", experiment.save = "nothing") 
{
  clf <- list() # create a classifier object
  clf$learner <- "terga1" # name of the method
  clf$params <- list() # parameter list
  clf$experiment <- list() # information about the experiment
  
  # POPULATION
  clf$params$sparsity           <- sparsity # number of non zero variables in the model
  clf$params$current_sparsity   <- NA # number of non zero variables in the model
  clf$params$size_pop           <- size_pop # how many models in the population to evolve
  clf$params$size_world         <- size_world # total number of variables
  clf$params$max.nb.features    <- max.nb.features
  clf$params$nb_generations     <- nb_generations # number of generation to evolve
  clf$params$unique_vars        <- unique_vars # weather in a model we can have one variable more than once
  clf$params$popSourceFile      <- popSourceFile
  clf$params$popSaveFile        <- popSaveFile
  
  # FITTING
  clf$params$objective          <- objective # prediction (roc) or regression (cor)
  clf$params$estimate_coefs     <- estimate_coefs # integer or real estimated coefficients.
  # SELECTION
  clf$params$select_type        <- select_type # selection method (mixte, elitist, tournoi)
  clf$params$select_perc1       <- select_perc1 # selection percentage for elitist
  clf$params$select_perc2       <- select_perc2 # selection percentage for tournoi
  clf$params$perc_best_ancestor <- perc_best_ancestor # percentage of best k-1 models to be used as seed for the k simulation
  # MUTATION
  clf$params$mutate_size        <- mutate_size # what percentage of models in the population are mutated
  clf$params$mutate_rate        <- mutate_rate # what percentage of the variables in the model are mutated
  # CONVERGENCE
  clf$params$convergence        <- convergence # what should the simulation stop when convergence ?
  clf$params$convergence_steps  <- convergence_steps # after how many steps without improvement do we consider convergence?
  clf$params$evolve_k1          <- evolve_k1 # weather to evolve models with k_1 or to search them exhaustively.
  # print out intermediary results
  clf$params$plot               <- plot           # plot results? 
  clf$params$verbose            <- verbose        # print out logs.
  clf$params$warnings           <- warnings       # print out warnings
  clf$params$debug              <- debug          # print out debugging information.
  clf$params$print_ind_method   <- print_ind_method # method to print individual
  
  # Computing options
  clf$params$nCores             <- nCores # parallel computing
  clf$params$parallel           <- nCores > 1 # parallel computing
  clf$params$parallelize.folds  <- parallelize.folds
  clf$params$parallel.local     <- FALSE
  clf$params$seed               <- seed
  
  # Evaluation options
  clf$params$evalToFit          <- evalToFit
  clf$params$k_penalty          <- k_penalty
  
  clf$params$language           <- language
  #clf$params$useCustomLanguage  <- useCustomLanguage
  clf$params$epsilon            <- epsilon
  clf$params$scoreFormula       <- scoreFormula
  clf$params$intercept          <- intercept

  # Experiment information
  if(!(experiment.id=="NULL"))
  {
    clf$experiment$id          <- experiment.id
  }else
  {
    clf$experiment$id          <- clf$learner
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
  
  return(clf)
}


# Launch the fit classifier routine
terga1_fit <- function(X, y, clf) {
  
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
  
  # if(clf$params$max.nb.features < nrow(X))
  # {
  #   if(clf$params$verbose) print(paste("... ... restricting X to the",clf$params$max.nb.features,"most significant variables"))  
  #   selected.features <- rownames(clf$feature.cor[order(clf$feature.cor$p),][1:min(clf$params$max.nb.features, nrow(X)),])
  #   X <- X[selected.features,]
  #   clf$coeffs_ <- clf$coeffs_[selected.features]
  #   clf$feature.cor <- clf$feature.cor[selected.features,]
  # }
  
  # set the size of the world
  clf$params$size_world <- nrow(X)
  
  # Print the experiment configuration
  if(clf$params$verbose) printClassifier(obj = clf)
  
  res <- list() # the final object containing the evolved models
  for(i in clf$params$sparsity) # sparsity is = k, i.e. the number of features in a model
  { 
    cat("\t\tResolving problem with\t", i, "\tvariables ...\n")
    
    # set the current_sparsity
    clf$params$current_sparsity <- i
    
    
    # test for
    if (clf$params$current_sparsity == 1 & !clf$params$evolve_k1) # if we want to evolve features for k_sparse=1 we create a normal population
    {  
      pop_last      <- as.list(1:nrow(X)) # create population with k_sparse = 1
      
    }else 
    {
      best_ancestor = NULL
      
      if (exists("pop_last")) 
      {
        if(length(evaluation) != 0)
        {
          best_ancestor = pop_last[[which.max(evaluation)]]  
        }
        
        # build a new population seeded by the last one
        pop         <- population(clf = clf, 
                                  size_ind = i, 
                                  size_world = nrow(X), 
                                  best_ancestor = best_ancestor, 
                                  size_pop = clf$params$size_pop,
                                  seed = clf$params$current_seed)
      }else 
      {
        # build a new population from scratch
        pop         <- population(clf = clf, 
                                  size_ind = i, 
                                  size_world = nrow(X), 
                                  best_ancestor = best_ancestor, 
                                  size_pop = clf$params$size_pop,
                                  seed = clf$params$current_seed)
      }
      
      # if this is population with model objects we transform in an index population
      if(isPopulation(obj = pop))
      {
        pop <- listOfModelsToListOfSparseVec(list.models = pop)
      }
      
      # then we evolve
      pop_last      <- evolve(X, y, clf, pop, seed = clf$params$current_seed)
    }
    
    # evaluate the fitting function for all the models of the populaion
    # transform to a population of model objects
    pop_last.mod <- listOfSparseVecToListOfModels(X = X, y = y, clf = clf, v = pop_last)
    # evaluate the population
    pop.last.eval <- evaluatePopulation(X, y, clf, pop_last.mod, force.re.evaluation = TRUE, eval.all = TRUE)
    # get the evaluation vector      
    evaluation    <- as.numeric(populationGet_X(element2get = "fit_", toVec = TRUE, na.rm = TRUE)(pop = pop.last.eval))
    # get the evaluation
    evaluation.ord <- order(abs(evaluation), decreasing = TRUE)
    # order by best in front
    evaluation <- evaluation[evaluation.ord]
    pop_ordered_mod <- pop.last.eval[evaluation.ord] # we gain speed
    
    # get the best individuals (should be the first)
    best_individual_index <- which.max(abs(evaluation))
    best_individual <- pop_ordered_mod[[best_individual_index]]
    
    # print it out
    if(clf$params$verbose) 
    {
      if(isModel(best_individual))
      {
        try(cat(paste("gen =",i,"\t", printModel(mod = best_individual, method = clf$params$print_ind_method, score = "fit_"),"\n")), silent = TRUE)
      }
    }
    
    # transform the indexes into models
    if(!isPopulation(obj = pop_ordered_mod))
    {
      pop_ordered_mod <- evaluatePopulation(X, y, clf, pop_ordered_mod, force.re.evaluation = TRUE, eval.all = TRUE)
    }
    
    # keep only models that are unique
    pop_ordered_mod <- unique(pop_ordered_mod)
    
    if(!(clf$params$popSaveFile=="NULL"))
    {
      #pop2Save      <- evaluatePopulation(X, y, clf, pop_ordered_mod, eval.all = TRUE)
      pop2Save      <- pop_ordered_mod
      savePopulation(pop2Save, paste("generation", i, clf$params$popSaveFile, sep = "_"))
    }
    
    # save the whole list of models ordered by fitting score. The best is the first
    res[[paste("k",i,sep = "_")]] <- pop_ordered_mod
  }
  return(res)
}
## End of the function terga_fit
#######################################################
