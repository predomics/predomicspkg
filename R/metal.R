


#' metal: metal searching algorithm
#' @description metal is a model search algorithm on a list of beam search approach and get the populations into GA.
#' @param sparsity: number of features in a given model. This is a vector with multiple lengths.
#' @param max.nb.features: focuses only on the subset of top most significant features (default:1000)
#' @param popSaveFile: (??)
#' @param saveFiles: ??
#' @param language is the language that is used by the different algorithms {bin, bininter, ter, terinter, ratio}, (default:"terinter")
#' @param scoreFormula: a Function that contains the ratio Formula or other specific ones
#' @param epsilon: a small value to be used with the ratio language (useCustomLanguage) (default: NULL). When null it is going to be calculated by the minimum value of X divided by 10.
#' @param objective: this can be auc, cor or aic. Terga can also predict regression, other than class prediction. (default:auc)
#' @param estimate_coefs: non ternary solution for the aic objective (default:FALSE)
#' @param evalToFit: The model performance attribute to use as fitting score (default:"fit_"). Other choices are c("auc_","accuracy_","precision_","recall_","f_score_")
#' @param k_penalty: Penalization of the fit by the k_sparsity (default: 0)
#' @param intercept: (??) (default:NULL)
#' @param testAllSigns: ??
#' @param plot: plot graphics indicating the evolution of the simulation (default:FALSE)
#' @param verbose: print out information on the progress of the algorithm (default:TRUE)
#' @param warnings: Print out warnings when runnig (default:FALSE).
#' @param debug: print debug information (default:FALSE)
#' @param print_ind_method: One of c("short","graphical") indicates how to print a model and subsequently a population during the run (default:"short").
#' @param parallelize.folds: parallelize folds when cross-validating (default:TRUE)
#' @param nCores: the number of cores to execute the program. If nCores=1 than the program runs in a non parallel mode
#' @param seed: the seed to be used for reproductibility. If seed=NULL than it is not taken into account (default:NULL).
#' @param experiment.id: The id of the experiment that is to be used in the plots and comparitive analyses (default is the learner's name, when not specified)
#' @param experiment.description: A longer description of the experiment. This is important when many experiments are run and can also be printed in by the printExperiment function.
#' @param experiment.save: Data from an experiment can be saved with different levels of completness, with options to be selected from c("nothing", "minimal", "full"), default is "minimal"
#' @param list.clfs: list of Genetor and Unificator
#' @param unificator.method: the default unificator is a terga2. Other one specified will yield a stop of the program.
#' @param unificator.evolver: the default evolve method used by the unificator which is by default a terga2.
#' @return an object containing a list of parameters for this classifier
#' @export
metal <- function(sparsity = 1:10, max.nb.features = 1000, 
                  #final.pop.perc = 100,  # percentage of models in the returned results
                  # population options
                  popSaveFile = "NULL", saveFiles = FALSE, pathSave="NULL",
                  # language in {bin, bininter, ter, terinter, ratio, mix}
                  language = "mix", 
                  # language options
                  scoreFormula = scoreRatio, epsilon = "NULL", 
                  # evaluation options
                  objective = "auc", k_penalty = 0, evalToFit = "accuracy_", estimate_coefs = FALSE, intercept = "NULL", testAllSigns = FALSE, 
                  # output options
                  plot = FALSE, verbose = TRUE, warnings = FALSE, debug = FALSE, print_ind_method = "short", 
                  parallelize.folds = TRUE,
                  # computing options
                  nCores = 10, seed = "NULL", #maxTime = Inf,
                  # experiment options
                  experiment.id = "NULL", experiment.description = "NULL", experiment.save = "nothing", 
                  # list of clfs
                  list.clfs="NULL", unificator.method = "terga2", unificator.evolver = "v2m_new")
{
  # standard means that we use the standard heuristics
  clf <- list()
  clf$learner <- "metal"
  clf$params  <- list()
  clf$experiment <- list() # information about the experiment
  clf$params$objective          <- objective
  clf$params$sparsity           <- sparsity
  clf$params$max.nb.features    <- max.nb.features
  
  clf$params$saveFiles          <- saveFiles   # It would be interesting to add this in the future
  clf$params$pathSave           <- pathSave   
  
  
  # print out intermediary results
  clf$params$plot               <- plot           # plot results? 
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
  #clf$params$final.pop.perc     <- final.pop.perc
  clf$params$evalToFit          <- evalToFit
  clf$params$k_penalty          <- k_penalty
  
  clf$params$language           <- language
  clf$params$popSaveFile        <- popSaveFile
  clf$params$epsilon            <- epsilon
  clf$params$scoreFormula       <- scoreFormula
  clf$params$intercept          <- intercept
  clf$params$list.clfs          <- list.clfs
  
  
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
    clf$experiment$description      <- experiment.description
  } else 
  {
    clf$experiment$description      <- paste(clf$learner, date(), sep = " ")
  }
  #match.arg(experiment.save) 
  clf$experiment$save               <- experiment.save
  
  # unificator information
  clf$params$unificator.method      <- unificator.method
  clf$params$unificator.evolver     <- unificator.evolver
  
  # add here the list of classifiers
  if(is.null(list.clfs))
  {
    languages                       <- c("terinter","bininter","ratio")
    learners                        <- c("terbeam","terda","terga2")
    generator_matrix2run            <- matrix(1, nrow = length(languages), ncol = length(learners))
    rownames(generator_matrix2run)  <- languages
    colnames(generator_matrix2run)  <- learners
    generator_matrix2run["ratio","terda"]       <- 0 # for terda
    # generate the default generators and unificators
    clf$params$list.clfs            <- generator_metal(mat = generator_matrix2run, clf=clf)
  }else
  {
    if(!is.list(list.clfs))
    {
      languages                       <- c("terinter","bininter","ratio")
      learners                        <- c("terbeam","terda","terga2")
      generator_matrix2run            <- matrix(1, nrow = length(languages), ncol = length(learners))
      rownames(generator_matrix2run)  <- languages
      colnames(generator_matrix2run)  <- learners
      generator_matrix2run["ratio","terda"]       <- 0 # for terda
      # generate the default generators and unificators
      clf$params$list.clfs            <- generator_metal(mat = generator_matrix2run, clf=clf)
    }else
    {
      clf$params$list.clfs            <- list.clfs  
    }
  }
  
  return(clf)
}


metal_fit <- function(X, y, clf)
{
  
  # TODO: sanity check
  list.clfs <- clf$params$list.clfs

  # get the number of clfs
  n <- length(list.clfs) # length of the attribute attribute
  if((n-1) != length(list.clfs[[n]]))
  {
    stop("metal_fit: number of clf objects  list does not correspond to the expected number")
  }

  if(clf$params$verbose) print("... ... sanity check of the list of clf")
  
  if(clf$params$verbose) printClassifier(obj = clf)
  # parallel switch
  #parfold <- FALSE
  # if(!clf$params$parallelize.folds & clf$params$parallel)
  # {
  #   parfold <- TRUE
  # }

  # use the attributes (attribute) to determine which clf is the unifier.
  ind.uni <- which(list.clfs[[n]]=="U")
  # Test wether this is a terga family algorithm
  if(list.clfs[[ind.uni]]$learner!='terga2')
  {
    stop("metal_fit: check the unificator to be of terga2 family")
  }
  u.clf <- list.clfs[[ind.uni]]
  # initiate the current sparsity
  u.clf$params$current_seed <- clf$params$current_seed
  
  size_pop <- u.clf$params$size_pop # the population size
  #u.clf$params$cluster <- clf$params$cluster

  #res_list        <- list() # generator results storage
  #size_pop_add    <- list() # a vector containing the number of models from each generator
  #popSourceFile   <- list() # a list of files to save populations of generators

  in_pop         <- list() #population of generator terbeams input for terga

  # for all the generators launch the learning to generate models
  for(i in c(1:(n-2)))
  {
    if(i != ind.uni)
    {
      g.clf                           <- list.clfs[[i]]
      g.clf$params$cluster            <- clf$params$cluster
      g.clf$params$parallelize.folds  <- FALSE
      # initiate the current sparsity
      g.clf$params$current_seed       <- clf$params$current_seed
      g.clf$params$objective          <- clf$params$objective
      
      # fix the max number of features
      g.clf$params$max.nb.features    <- clf$params$max.nb.features
      # don't compute importance for generators
      g.clf$params$compute.importance <- FALSE
      
      # copy the feature.cor from the main clf to the generators
      if(!is.null(clf$feature.cor)) 
      {
        g.clf$feature.cor <- clf$feature.cor
        if(clf$params$verbose) print("... ... found feature.cor adding to generator clf")
      }
      # copy the coeffs_ from the main clf to the generators
      if(!is.null(clf$coeffs_)) 
      {
        g.clf$coeffs_ <- clf$coeffs_
        if(clf$params$verbose) print("... ... found coeffs_ adding to generator clf")
      }
      
      switch(g.clf$learner,
             terBeam={  res       <- terBeam_fit(X, y, g.clf)},
             terda={    res       <- terda_fit(X, y, g.clf)},
             terga2={   res       <- terga2_fit(X, y, g.clf)}
      )

      if((g.clf$params$epsilon=="NULL"))
      {
        g.clf$params$epsilon      <- 0
      }

      NBest                       <- min(floor(size_pop / (n-2) ), length(modelCollectionToPopulation(res)))
      #bestIndividuals             <- getNBestIndividuals(res, g.clf, NBest, equal.sparsity = TRUE)
      
      bestIndividuals             <- getNBestModels(obj = res, 
                                                    significance = FALSE, 
                                                    by.k.sparsity = TRUE,
                                                    k.penalty = 0,
                                                    n.best = NBest,
                                                    single.best = FALSE,
                                                    single.best.cv = TRUE,
                                                    single.best.k = NULL,
                                                    max.min.prevalence = FALSE,
                                                    X = NULL,
                                                    verbose = FALSE, 
                                                    evalToOrder = "fit_",
                                                    return.population = TRUE,
                                                    unique.control = TRUE
      )
      
      # when the generators did not generate anything
      if(!is.null(bestIndividuals))
      {
        bestIndividuals           <- evaluatePopulation(X, y, g.clf, bestIndividuals, force.re.evaluation = TRUE, eval.all = TRUE)
        bestIndividuals           <- sortPopulation(pop = bestIndividuals, evalToOrder = "fit_")
        #pop <- sortPopulation(modelCollectionToPopulation(res$models))[1:floor(size_pop/2/(n-2))]
        in_pop[(length(in_pop) +1):(length(in_pop) + length(bestIndividuals))] <- bestIndividuals
      }
      
      if(clf$params$verbose) print(paste("... ... generator nb", i,"ended"))
      
      # test save if needed
      #savePopulation(bestIndividuals, paste(i,"saved_NBest",id,sep = '.'),compress = FALSE)
      #popSourceFile[i] <- paste(paste(i,"saved_NBest",id,sep = '.'),'yml',sep = '.')
      #size_pop_add[i] <- length(bestIndividuals)#floor(size_pop/2/n/length(clf$params$sparsity))*length(clf$params$sparsity)
    } # end for all generators
  } # end for

  # compute the size of the random population, knowing the size of the generated ones to be merged
  #u.clf$params$size_pop_random <- size_pop - Reduce("+", size_pop_add)
  #u.clf$params$popSourceFile   <- popSourceFile
  #in_pop <- evaluatePopulation(X, y, u.clf, in_pop, eval.all = TRUE,force.re.evaluation = TRUE)

  u.clf$params$size_pop           <- length(in_pop)
  u.clf$params$in_pop             <- in_pop
  u.clf$params$cluster            <- clf$params$cluster
  u.clf$params$objective          <- clf$params$objective
  
  # fix the max number of features
  u.clf$params$max.nb.features    <- clf$params$max.nb.features
  
  if(!is.null(clf$feature.cor)) 
  {
    u.clf$feature.cor <- clf$feature.cor
    if(clf$params$verbose) print("... ... found feature.cor adding to unificator clf")
  }
  if(!is.null(clf$coeffs_)) 
  {
    u.clf$coeffs_ <- clf$coeffs_
    if(clf$params$verbose) print("... ... found coeffs_ adding to unificator clf")
  }

  # run learner Unificator
  if(clf$params$verbose) print(paste("... ... launching the unificator"))
  res.cv   <- terga2_fit(X, y, u.clf)
  if(clf$params$verbose) print(paste("... ... unificator has finished"))
  
  if(!is.null(res.cv))
  {
    return(res.cv)
  }else
  {
    return(NULL)
  }
}
