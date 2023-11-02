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
# @script:  predomics.R
# @author:  Edi Prifti
# @author:  Lucas Robin
# @author:  Shasha Cui
# @author:  Blaise Hanczar
# @author:  Yann Chevaleyre
# @author:  Jean-Daniel Zucker
# @date:    May 2016
# @date:    October 2023
################################################################

#' fit: runs the classifier on a dataset
#'
#' @import doSNOW
#' @import snow
#' @import doRNG
#' @description This function runs a learning experiment based on the classifier
#' object and the given dataset.
#' @param X: Dataset to classify
#' @param y: Variable to predict
#' @param clf: The classifier object object containing the settings of the classifier
#' @param cross.validate: Whether or not the classification should be run in
#' cross-validation mode (default:TRUE)
#' @param lfolds: The folds to be used for the cross-validation
#' @param nfolds: The number of folds to use in the cross-validation. If lfolds
#' are not specified this option allows to set them up (default:10)
#' @param parallelize.folds: Switch setting the parallelization mode based on
#' cross-validation folds and nothing else in the algorithm (default:TRUE).
#' This is much more efficient.
#' @param compute.importance: The importance of variables in the learning process
#' during cross-validation can be computed. This is based on data perturbation
#' similar to the mean decrease accuracy in the random forest algorithm. Moreover,
#' this gives feature prevalence in models during CV (default:TRUE)
#' @param return.all: Should all results from the cross-validation steps be
#' returned. This is usually needed when testing stability of the models
#' (default:FALSE)
#' @param log.file: The output file for parallel logs (default:'parallel.log')
#' @param path: The path where to save temporary data
#' @return An experiment object containing the classifier along with the
#' classification results as a sub-element
#' @export
fit <- function(X,
                y,
                clf,
                cross.validate = FALSE,
                lfolds = NULL,
                nfolds = 10,
                parallelize.folds = TRUE,
                compute.importance = TRUE,
                return.all = FALSE,
                log.file = "parallel.log",
                path = NULL)
{
  # test the classifier object
  if (!isClf(clf))
  {
    stop("fit: please provide a valid classifier object!")
  }
  
  clf$params$compute.importance <- compute.importance
  check.X_y_w(X, y, w = NULL)
  
  if (!is.matrix(X))
  {
    cat("... Database X is not a matrix! Converting ...\n")
    
    # log_info("... Database X is not a matrix! Converting ...")
    
    # transforming to dataframe first will make sure to keep the dimnames such as for instance when it is sparse table or otu_table (phyloseq)
    X <- as.matrix(as.data.frame(X))
  }
  
  # if no rownames add some
  if (is.null(rownames(X)))
  {
    rownames(X) <- paste("F", c(1:nrow(X)), sep = "_")
    cat("... No names for X variables: renaming by default\n")
  }
  
  if (ncol(X) != length(y))
  {
    stop("fit: the number of observations in 'X' is diffrent from the one in 'y' ")
  }
  
  if (clf$params$objective == "auc")
  {
    cat("... Classification mode, computing factor(y) for speedup and robustness\n")
    y <- factor(as.character(y))
  }
  
  # The possibility to select and focus on the top best features of X
  if (!is.null(clf$params$max.nb.features))
  {
    max.nb.features <-
      min(nrow(X), clf$params$max.nb.features, na.rm = TRUE)
  } else
  {
    max.nb.features <-  nrow(X)
  }
  
  # set a flag
  path.feature.cor <- paste(path, "feature.cor.rda", sep = "")
  
  if (file.exists(path.feature.cor))
  {
    cat("... Loading feature correlation for speedup\n")
    # restore precomputed object from the hard drive
    load(path.feature.cor)
    cat("... Correlation file loaded\n")
    
    if (exists("feature.cor"))
    {
      if (any(is.na(match(
        rownames(feature.cor), rownames(X)
      ))))
      {
        stop(
          paste(
            "... feature.cor does not match X... needs to be recomputed. You can delete the following file",
            path.feature.cor
          )
        )
      }
    } else
    {
      stop("feature.cor object does not exist. Please check the name.")
    }
  } else
  {
    cat("... Computing feature correlation for speedup\n")
    if (clf$params$objective == "cor")
      # correlation
    {
      # if any NA in y omit them
      if (any(is.na(y)))
      {
        ina <- is.na(y)
        cat(
          paste(
            "... y contains ",
            sum(ina),
            "NA values ... ommiting the observations in both y and X\n"
          )
        )
        # transforming to dataframe first will make sure to keep the dimnames such as for instance when it is sparse table or otu_table (phyloseq)
        y.nona <- y[!ina]
        X.nona <- X[, !ina]
      } else
      {
        y.nona <- y
        X.nona <- X
      }
      
      feature.cor     <-
        filterfeaturesK(
          data = t(apply(X.nona, 1, rank)),
          # for speedup
          trait = rank(y.nona),
          # for speedup
          k = max.nb.features,
          type = "pearson",
          # for speedup
          sort = TRUE,
          verbose = clf$params$verbose,
          return.data = FALSE
        ) # to avoid having to recompute this all the time
    } else
      # classification
    {
      feature.cor     <- filterfeaturesK(
        data = X,
        trait = y,
        k = max.nb.features,
        type = "wilcoxon",
        sort = TRUE,
        verbose = clf$params$verbose,
        return.data = FALSE
      ) # to avoid having to recompute this all the time
    }
    
    if (!is.null(path))
    {
      save(feature.cor, file = path.feature.cor)
      cat("... Correlation file saved\n")
    }
  }
  
  cat("... Storing data in the classifier object for speedup\n")
  clf$feature.cor <- feature.cor # add them to the clf
  
  if (all(is.na(clf$feature.cor$p)))
  {
    warning("runClassifier: does not seem to have produced a pvalue")
  }
  
  # store the initial order and indexes
  clf$data          <- list()
  clf$data$features <- as.character(rownames(X))
  names(clf$data$features) <- clf$data$features
  
  # select the top "min(nrow(X), clf$params$max.nb.features)" X features
  X <- X[rownames(clf$feature.cor)[1:max.nb.features], ]
  
  clf$data$X        <- X
  clf$data$X.min    <- min(X, na.rm = TRUE)
  clf$data$X.max    <- max(X, na.rm = TRUE)
  clf$data$y        <- y
  
  # compute the coefficients once for all to improve performance
  cat("... Computing ternary coefficients for speedup\n")
  coeffs          <-
    getSign(
      X = X,
      y = y,
      clf = clf,
      parallel.local = FALSE
    )
  clf$coeffs_     <- coeffs # add them to the clf
  
  
  # check sparsity not to be larger than variables in X
  if (any(is.na(clf$params$sparsity > nrow(X))))
  {
    # adding the maximum number of featuers
    clf$params$sparsity <- c(clf$params$sparsity, nrow(X))
  }
  # mark NAs the bigger ones
  clf$params$sparsity[clf$params$sparsity > nrow(X)] <- NA
  # delete them
  clf$params$sparsity <-
    clf$params$sparsity[!is.na(clf$params$sparsity)]
  
  # set the seed while sanity checking
  if (!(any(clf$params$seed == "NULL")))
  {
    # convert from list to a vector
    if (is.list(clf$params$seed))
    {
      clf$params$seed <- unlist(clf$params$seed)
      if (any(class(clf$params$seed) != "numeric"))
      {
        stop("fit: convertion of seed from list to numeric vector failed.")
      }
    }
    
    # we can have multiple k-folds per experiment if the seed is vectors of seeds
    if (length(clf$params$seed) == 1)
    {
      set.seed(clf$params$seed)
      cat("... One seed found, setting by default\n")
    } else
    {
      if (length(clf$params$seed) == 0)
      {
        stop("fit: the seed should have at least one value.")
      }
      
      # if we are here this means that it everything is as expected seed is a
      # vector of numeric values.
      set.seed(clf$params$seed[1])
      cat("... Multiple seeds found, setting the default\n")
    }
  }
  
  # check and set the number of folds
  if (!is.null(lfolds))
  {
    cat("... Custom folds are provided\n")
    if (!is.list(lfolds))
    {
      lfolds = NULL
    } else
      #if lfolds exists
    {
      nfolds = length(lfolds)
    }
  }
  
  if (nfolds == 1)
  {
    cat(
      "... The number of folds is set to 1. In this case I'm deactivating the cross-validation process\n"
    )
    cross.validate = FALSE
  }
  
  # set the parallelize.folds parameter. If no crossval than it is deactivated
  clf$params$parallelize.folds <-
    parallelize.folds & cross.validate & clf$params$parallel
  
  # add a parallel.local parameter if we wish to speed out some local steps
  clf$params$parallel.local <- FALSE
  
  # START CLUSTER if parallel computing set the cluster
  if (clf$params$parallel)
  {
    cat(
      paste(
        "... Running the classifier",
        clf$learner,
        "in parallel mode with ",
        clf$params$nCores + 1,
        "CPUs\n"
      )
    )
    # adding another core for the whole dataset. If it is taken into account
    # during the launch this will be a sleeping thread so no harm, if not it
    # will allow to run faster as we won't forget to increment it
    registerDoSNOW(
      clf$params$cluster <-
        makeCluster(clf$params$nCores + 1, type = "SOCK", outfile = log.file)
    )
    if (!clf$params$parallelize.folds)
      # if folds are not parallelized
    {
      clf$params$parallel.local <-
        TRUE # switch the local parallel to TRUE
    }
  } else
  {
    cat(paste(
      "... Running the classifier",
      clf$learner,
      "with a single CPU\n"
    ))
  }
  
  # for all the predomics learners
  if (!isLearnerSota(clf))
    # if runing the BTR algorithms
  {
    # set the epsilon
    if (!is.null(clf$params$epsilon))
    {
      if (clf$params$epsilon == "NULL")
      {
        #clf$params$epsilon        <- .Machine$double.xmin
        clf$params$epsilon        <- 0
        if (clf$params$verbose)
          cat("... Setting epsilon for predomics learners\n")
      }
    }
    
    # # for regression rank y since the beginning to avoid doing it for every model
    # if(clf$params$objective == "cor")
    # {
    #   if(clf$params$verbose) cat("... Correlation mode, computing rank(y) for speedup\n")
    #   # for spearman implementation to scale-up. Deactivatet at the current moment
    #   y <- rank(y)
    # }
  }
  
  # save the data step by step to be able to resume
  if (clf$experiment$save != "nothing")
  {
    if (clf$params$verbose)
      cat("... Saving experiments\n")
    fileNames                 <- gsub(" ", "_", clf$experiment$id)
    dir.create(fileNames)
    setwd(fileNames)
    saveResults(X, paste("X", fileNames, "yml", sep = "."))
    saveResults(y, paste("Y", fileNames, "yml", sep = "."))
    experiment <- list()
    experiment$desc           <- clf$experiment$description
    experiment$params         <-
      list(clf = list(
        learner = clf$learner,
        params = clf$params,
        experiment = clf$experiment
      ))
    if (clf$experiment$save == "full")
    {
      if ((clf$params$popSaveFile == "NULL"))
      {
        clf$params$popSaveFile <- fileNames
      }
    }
  }
  
  switch(
    clf$learner,
    terda =
      {
        # Here we handle also the sota.glmnet as well since this is a derivate of terda
        cat('... terda fitting based on Linear programming relaxation ...\n')
      },
    terga1 =
      {
        cat('... First version of terga fitting based on Genetic Algorithm heuristics ...\n')
      },
    terga2 =
      {
        cat(
          '... Second and faster version of terga fitting based on Genetic Algorithm heuristics ...\n'
        )
      },
    terBeam =
      {
        cat('... terbeam fitting based on Exhaustive Heuristic beam search ...\n')
      },
    metal =
      {
        cat('... model fitting based on aggregating different Heuristics ...\n')
      },
    sota.svm =
      {
        cat('... SOTA: state of the art SVM fitting ...\n')
      },
    sota.rf =
      {
        cat('... SOTA: state of the art Ranfom Forest fitting ...\n')
      },
    {
      warning('This method does not exist !')
    }
  )
  
  # open plot file
  if (clf$params$plot)
  {
    pdf(file = paste0("graphics_from_", clf$experiment$description, ".pdf"))
  }
  
  # lancer le classifier
  if (!cross.validate)
  {
    cat("... No cross validation mode\n")
    res.clf               <- list()
    res.clf$classifier    <- runClassifier(X, y, clf)
    
  } else
  {
    cat("... Cross validation mode\n")
    res.clf               <- list()
    res.clf$classifier    <- list()
    
    res.clf$crossVal      <-
      runCrossval(
        X,
        y,
        clf,
        lfolds = lfolds,
        nfolds = nfolds,
        return.all = return.all
      )
    
    # store the whole dataset results in the right place
    res.clf$classifier    <- res.clf$crossVal$whole
    # unweight the crossVal object that played the role of the transporter
    res.clf$crossVal      <-
      res.clf$crossVal[-match("whole", names(res.clf$crossVal))]
    res.clf$lfolds        <- lfolds
    
    # add FEATURE IMPORTANCE to the model objects based on the crossval feature importance
    if (compute.importance & cross.validate)
    {
      if (!is.null(res.clf$crossVal$fip))
      {
        feature.importance.cv <-
          rowMeans(res.clf$crossVal$fip$mda, na.rm = TRUE)
        feature.prevalence.cv <-
          res.clf$crossVal$fip$fpf / ncol(res.clf$crossVal$fip$mda) # normalize %
        if (!is.null(res.clf$classifier$fip$mda))
        {
          feature.importance  <- res.clf$classifier$fip$mda
        } else
        {
          feature.importance  <- rep(NA, length(clf$data$features))
          names(feature.importance) <- names(clf$data$features)
        }
        
        if (isModelCollection(res.clf$classifier$models))
        {
          pop <- modelCollectionToPopulation(res.clf$classifier$models)
        } else
        {
          pop <- NULL
        }
        
        if (!isPopulation(pop))
        {
          cat("... Bad news! The population of models seems to be empty. No models were found\n")
          res.clf$classifier$models <-
            listOfModels2ModelCollection(pop)
          cat("... Thank you for using Predomics\n")
          return(res.clf)
        }
        
        # for each model add a vector with feature importance
        for (i in 1:length(pop))
        {
          # MDA generalization
          pop[[i]]$mda.cv_ <- rep(0, length(pop[[i]]$names_))
          names(pop[[i]]$mda.cv_) <- pop[[i]]$names_
          ind.features <-
            intersect(pop[[i]]$names_, names(feature.importance.cv))
          pop[[i]]$mda.cv_[ind.features] <-
            feature.importance.cv[ind.features]
          
          # prevalence in top models in the folds
          pop[[i]]$prev.cv_ <- rep(0, length(pop[[i]]$names_))
          names(pop[[i]]$prev.cv_) <- pop[[i]]$names_
          ind.features <-
            intersect(pop[[i]]$names_, names(feature.prevalence.cv))
          pop[[i]]$prev.cv_[ind.features] <-
            feature.prevalence.cv[ind.features]
          
          # MDA empirical
          pop[[i]]$mda_ <- rep(0, length(pop[[i]]$names_))
          names(pop[[i]]$mda_) <- pop[[i]]$names_
          ind.features <-
            intersect(pop[[i]]$names_, names(feature.importance))
          pop[[i]]$mda_[ind.features] <-
            feature.importance[ind.features]
        }
        
        # convert to model collection and put back
        res.clf$classifier$models <-
          listOfModels2ModelCollection(pop)
      }
    }
    
    if (clf$experiment$save != "nothing")
    {
      experiment$params$lfolds <- lfolds
    }
  } # end cross.validate test
  
  cat("... Learning process is finished succesfuly\n")
  
  # STOP CLUSTER
  if (clf$params$parallel)
  {
    stopCluster(clf$params$cluster)
    cat("... Stopping the parallel cluster\n")
  }
  
  # SAVE EXPERIMENT
  if (clf$experiment$save != "nothing")
  {
    experiment$results    <- res.clf
    saveResults(experiment, paste(fileNames, "yml", sep = "."))
    setwd("..")
    cat("... Saving exmeriment finished\n")
  }
  
  # closing graphics
  if (clf$params$plot)
  {
    dev.off()
  }
  
  cat("... Thank you for using Predomics. Don't forget to digest the results now.\n")
  
  class(res.clf) <- c("experiment", "predomics")
  
  return(res.clf)
}



#' Runs the learning on a dataset
#' @export
#' @import doSNOW
#' @import snow
#' @description This function runs a classifier in a given dataset
#' @param X: The dataset to classify
#' @param y: The variable to predict
#' @param clf: The classifier object containing the different settings of the classifier.
#' @param x_test: if not NULL (default) this dataset will be used to evaluate the models in a subset for the feature importance
#' @param y_test: if not NULL (default) this dataset will be used to evaluate the models in a subset for the feature importance
#' @return the classifier along with the classification results as a sub-element
runClassifier <- function(X,
                          y,
                          clf,
                          x_test = NULL,
                          y_test = NULL)
{
  if (clf$params$verbose)
    cat("... Entering runClassifier\n")
  # test the classifier object
  if (!isClf(clf))
  {
    stop("fit: please provide a valid classifier object!")
  }
  
  check.X_y_w(X, y, w = NULL)
  
  if (clf$params$verbose)
    cat("... Storing data in classifier object for speedup\n")
  
  # set the epsion
  if (!is.null(clf$params$epsilon))
  {
    if (clf$params$epsilon == "NULL")
    {
      #clf$params$epsilon        <- min(X[X!=0], na.rm = T)/10
      #clf$params$epsilon        <- .Machine$double.xmin
      clf$params$epsilon        <- 0
    }
  }
  
  # when not in crossval this will be null so we need to initiate this
  if (is.null(clf$params$current_seed))
  {
    clf$params$current_seed <- clf$params$seed[1]
  }
  
  if (clf$params$debug)
    cat("=> DBG: before fit\n")
  
  startingTime <- Sys.time()
  switch(
    clf$learner,
    terda =
      {
        # Here we handle also the sota.glmnet as well since this is a derivate of terda
        if (clf$params$verbose)
          cat('... terda fitting based on Linear programming relaxation ...\n')
        res <- terda_fit(X, y, clf)
      },
    terga1 =
      {
        if (clf$params$verbose)
          cat('... First version of terga fitting based on Genetic Algorithm heuristics ...\n')
        res <- terga1_fit(X, y, clf)
      },
    terga2 =
      {
        if (clf$params$verbose)
          cat(
            '... Second and faster version of terga fitting based on Genetic Algorithm heuristics ...\n'
          )
        res <- terga2_fit(X, y, clf)
      },
    terBeam =
      {
        if (clf$params$verbose)
          cat('... terbeam fitting based on Exhaustive Heuristic beam search ...\n')
        res <- terBeam_fit(X, y, clf)
      },
    metal =
      {
        if (clf$params$verbose)
          cat('... model fitting based on aggregating different Heuristics ...\n')
        res <- metal_fit(X, y, clf)
      },
    sota.svm =
      {
        if (clf$params$verbose)
          cat('... SOTA: state of the art SVM fitting ...\n')
        res <- sota.svm_fit(X, y, clf)
      },
    sota.rf =
      {
        if (clf$params$verbose)
          cat('... SOTA: state of the art Ranfom Forest fitting ...\n')
        res <- sota.rf_fit(X, y, clf)
      },
    {
      warning('This method does not exist !')
    }
  )
  
  if (clf$params$debug)
  {
    cat("=> DBG: after fit\n")
    printy(res)
  }
  
  if (isModelCollection(res))
  {
    clf$models <- res
  } else
  {
    clf$models <- NULL
    warning("runClassifier: major issue - no models produced ... stoping")
  }
  
  # FEATURE IMPORTANCE
  if (clf$params$compute.importance & !isLearnerSota(clf))
  {
    # the population of models that is learned in this fold
    pop.fold <-
      modelCollectionToPopulation(mod.collection = clf$models)
    
    # In the case where we are in the whole dataset and x_test is null as is y_test we use X as a whole to learn the feature importance.
    # In this case we will have an empirical importance
    if (is.null(x_test) | is.null(y_test))
    {
      warning(
        "runClassifier: setting test data to train. We shouldn't be here... please investigate"
      )
      x_test <- X
      y_test <- y
    }
    
    # for each of the best models in the population compute the importance in the test population
    efip.fold <-
      evaluateFeatureImportanceInPopulation(
        pop = pop.fold,
        X = x_test,
        y = y_test,
        clf = clf,
        score = "fit_",
        filter.ci = TRUE,
        method = "extensive",
        seed = c(1:10),
        # 10 times the perturbation for more accurate importance
        aggregation = "mean",
        verbose = clf$params$verbose
      )
    clf$fip <- efip.fold
  }
  
  if (clf$params$debug)
  {
    cat("=> DBG: after efip\n")
  }
  
  if (isModelCollection(clf$models))
  {
    # update the final indexes as the input X
    clf$models <-
      updateObjectIndex(obj = clf$models,
                        features = clf$data$features)
  }
  
  clf$execTime <-
    as.numeric(Sys.time() - startingTime, units = "mins")
  
  if (clf$params$debug)
  {
    cat("=> DBG: after end runclassifier\n")
  }
  
  return(clf)
}


#' Compute the cross-validation emprirical and generalization scores
#'
#' @description Compute the cross-validation emprirical and generalization scores.
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the response vector
#' @param clf: the classifier parameter object
#' @param nfolds: the number of folds for the cross-validation
#' @param return.all: return all results from the crossvalidation for feature stability testing
#' @return a list containing empirical, generalisation scores for each fold as well as a matrix with the mean values.
#' @export
runCrossval <- function(X,
                        y,
                        clf,
                        lfolds = NULL,
                        nfolds = 10,
                        return.all = FALSE)
{
    # test the classifier object
    if (!isClf(clf))
    {
      stop("fit: please provide a valid classifier object!")
    }
    
    check.X_y_w(X, y, w = NULL)
    
    # if folds exist
    if (!is.null(lfolds))
    {
      # add the whole dataset
      lfolds            <- c(NA, lfolds)
      names(lfolds)[1]  <- "whole"
      # count
      nfolds <- length(lfolds)
      if (clf$params$verbose)
        cat("... Folds are provided\n")
    } else
    {
      # if they don't exist create them
      if (!myAssertNotNullNorNa(nfolds))
      {
        stop("runCrossval: lfolds or nfolds does not exist, please provide!")
      }
      
      if (length(clf$params$seed) > 1)
      {
        lfolds <- list()
        for (i in 1:length(clf$params$seed))
        {
          lfolds.tmp      <-
            create.folds(
              y[!is.na(y)],
              k = nfolds,
              list = TRUE,
              returnTrain = FALSE,
              seed = clf$params$seed[i]
            ) # Index of 10-fold
          names(lfolds.tmp) <-
            paste(names(lfolds.tmp), "_r", i, sep = "") # modify names
          lfolds          <- c(lfolds, lfolds.tmp)
        }
      } else
        # if only once
      {
        lfolds            <-
          create.folds(
            y[!is.na(y)],
            k = nfolds,
            list = TRUE,
            returnTrain = FALSE,
            seed = clf$params$seed
          ) # Index of 10-fold
      }
      
      clf$lfolds <- lfolds
      
      # add the whole dataset as a fold to process at the same time in //
      lfolds              <- c(NA, lfolds)
      names(lfolds)[1]    <- "whole"
      # count
      nfolds              <- length(lfolds)
      
      if (clf$params$verbose)
        cat(paste("... All the ", nfolds, " CV folds are set\n"))
    }
    
    # create the empty object structure that will be returned
    res.crossval                                <- list()
    res.crossval$nfold                          <- list()
    res.crossval$k                              <- list()
    res.crossval$scores                         <- list()
    res.crossval$scores$empirical.auc           <-
      as.data.frame(matrix(
        nrow = max(clf$params$sparsity),
        ncol = (nfolds - 1)
      ))
    rownames(res.crossval$scores$empirical.auc) <-
      c(paste("k", c(1:max(
        clf$params$sparsity
      )), sep = "_"))
    colnames(res.crossval$scores$empirical.auc) <-
      paste("fold", 1:(nfolds - 1), sep = "_") #-1 since the whole won't be taken into account
    
    # add others using the same model
    res.crossval$scores$generalization.auc      <-
      res.crossval$scores$empirical.auc # auc
    # accuracy
    res.crossval$scores$empirical.acc           <-
      res.crossval$scores$empirical.auc # accuracy
    res.crossval$scores$generalization.acc      <-
      res.crossval$scores$empirical.auc # accuracy
    # recall
    res.crossval$scores$empirical.rec           <-
      res.crossval$scores$empirical.auc # recall
    res.crossval$scores$generalization.rec      <-
      res.crossval$scores$empirical.auc # recall
    # precision
    res.crossval$scores$empirical.prc           <-
      res.crossval$scores$empirical.auc # precision
    res.crossval$scores$generalization.prc      <-
      res.crossval$scores$empirical.auc # precision
    # f1-score
    res.crossval$scores$empirical.f1s           <-
      res.crossval$scores$empirical.auc # f1-score
    res.crossval$scores$generalization.f1s      <-
      res.crossval$scores$empirical.auc # f1-score
    # correlation
    res.crossval$scores$empirical.cor           <-
      res.crossval$scores$empirical.auc # cor
    res.crossval$scores$generalization.cor      <-
      res.crossval$scores$empirical.auc # cor
    
    # for each fold compute the best models
    if (clf$params$parallelize.folds)
    {
      cat("... Starting CV in parallel\n")
      # execute each crossVal in //
      
      res.all <- foreach(i = 1:nfolds) %dorng%
        {
          # prepare the datasets
          if (all(is.na(lfolds[[i]])))
          {
            # Here we are in the whole dataset
            # training dataset
            x_train = X
            y_train = y
            # # testing dataset, needed for checking extreme cases
            x_test = X
            y_test = y
            
          } else
          {
            # Here we are in the whole dataset
            # training dataset
            x_train = X[, -lfolds[[i]]]
            y_train = y[-lfolds[[i]]]
            # # testing dataset needed for checking extreme cases
            if (length(lfolds[[i]]) == 1)
              # for leave one out
            {
              x_test = as.matrix(X[, lfolds[[i]]])
            } else
            {
              x_test = X[, lfolds[[i]]]
            }
            y_test = y[lfolds[[i]]]
            
          } # end else other folds
          
          # omit some particular cases
          if (any(table(y_train) == 0) |
              length(table(y_train)) == 1)
            # to take into account the leve one out case
          {
            NULL
          } else
          {
            # TODO: test saveing all the k-folds during execution
            if (!(clf$params$popSaveFile == "NULL"))
            {
              #dirName <- paste(clf$params$popSaveFile,paste("crossVal.fold", i, sep = ""),sep="/")
              dirName                         <-
                paste("crossVal.fold", i, sep = "")
              dir.create(dirName)
              setwd(dirName)
            }
            
            # Launch the classifier in the train dataset and digest results
            clf$params$current_seed <-
              clf$params$seed[1] + i # set the current seed
            if (clf$params$debug)
            {
              cat("=> DBG: before runclassifier\n")
              
              runClassifier(
                X = x_train,
                y =  y_train,
                clf = clf,
                x_test = x_test,
                y_test = y_test
              )
            } else
            {
              try({
                runClassifier(
                  X = x_train,
                  y =  y_train,
                  clf = clf,
                  x_test = x_test,
                  y_test = y_test
                )
              }, silent = TRUE)
            }
          }
        } # end of folds loop (foreach)
      
    } else
      # no parallel
    {
      # execute each crossval in serial
      cat("... Starting cross validation not in parallel\n")
      
      res.all <- list()
      for (i in 1:nfolds)
      {
        if (clf$params$verbose) {
          cat("===> k-fold\t", names(lfolds)[i], "\n")
        }
        
        # prepare the datasets
        if (all(is.na(lfolds[[i]])))
        {
          # Here we are in the whole dataset
          # training dataset
          x_train = X
          y_train = y
          # # testing dataset, needed for checking extreme cases
          x_test = X
          y_test = y
          
        } else
        {
          # Here we are in the whole dataset
          # training dataset
          x_train = X[, -lfolds[[i]]]
          y_train = y[-lfolds[[i]]]
          # # testing dataset needed for checking extreme cases
          if (length(lfolds[[i]]) == 1)
            # for leave one out
          {
            x_test = as.matrix(X[, lfolds[[i]]])
          } else
          {
            x_test = X[, lfolds[[i]]]
          }
          y_test = y[lfolds[[i]]]
        } # end other folds
        
        # omit some particular cases
        if (any(table(y_train) == 0) |
            length(table(y_train)) == 1)
          # to take into account the leve one out case
        {
          warning("runCrossval: only one level in the class impossible to compute fitness")
          next
        }
        
        # TODO: test saveing all the k-folds during execution
        if (!(clf$params$popSaveFile == "NULL"))
        {
          #dirName <- paste(clf$params$popSaveFile,paste("crossVal.fold", i, sep = ""),sep="/")
          dirName                         <-
            paste("crossVal.fold", i, sep = "")
          dir.create(dirName)
          setwd(dirName)
        }
        
        # Launch the classifier in the train dataset and digest results
        clf$params$current_seed           <- clf$params$seed[1] + i
        
        if (clf$params$debug)
        {
          res.all[[i]]                    <-
            runClassifier(
              X = x_train,
              y =  y_train,
              clf = clf,
              x_test = x_test,
              y_test = y_test
            )
        } else
        {
          res.all[[i]]                    <- try({
            runClassifier(
              X = x_train,
              y =  y_train,
              clf = clf,
              x_test = x_test,
              y_test = y_test
            )
          }, silent = FALSE)
          # end try
        } # end else debug
      } # end of folds loop (for)
    } # end else parallelize.folds
    
    if (clf$params$verbose)
      cat("... Cross validation finished\n")
    
    # store the empirical result separately
    res.crossval$whole <- res.all[[1]]
    # omit it from the empirical results as they will be extracted separately and
    res.all <- res.all[-1]
    # also clean the lfolds object
    lfolds  <- lfolds[-1]
    
    # results for FEATURE IMPORTANCE
    # the MDA for each fold
    mda.all <-
      as.data.frame(matrix(NA, nrow = nrow(X), ncol = length(res.all)))
    rownames(mda.all) <-
      rownames(X)
    colnames(mda.all) <- names(res.all)
    # create also results for the standard deviation and the prevalence in the folds
    pda.all <- sda.all <- mda.all
    
    # DISPATCH the results in the custom output structure
    for (i in 1:length(res.all))
    {
      # prepare the datasets
      # training dataset
      x_train = X[, -lfolds[[i]]]
      y_train = y[-lfolds[[i]]]
      # testing dataset
      if (length(lfolds[[i]]) == 1)
        # for leave one out
      {
        x_test = as.matrix(X[, lfolds[[i]]])
      } else
      {
        x_test = X[, lfolds[[i]]]
      }
      y_test = y[lfolds[[i]]]
      
      res_train <- res.all[[i]]
      
      if (is.list(res_train))
        # if result object exist
      {
        if (is.null(res_train$models))
          # and is a model collection
        {
          if (!isModelCollection(res_train))
          {
            res_train.digest              <- NULL
          } else
          {
            # digest
            res_train.digest              <-
              digestModelCollection(obj = res_train,
                                    X = x_train,
                                    clf = clf)
          }
        } else
          # is a crossval result
        {
          res_train.digest                <-
            digestModelCollection(obj = res_train$models,
                                  X = x_train,
                                  clf = clf)
        } # end if/else models exist
      } else
        # return nothing
      {
        res_train.digest                  <- NULL
      }
      
      # if the results could be digested
      if (!is.null(res_train.digest))
      {
        # for all the best models of each k-sparse (create empty matrix) for auc
        res.crossval$k$auc                <-
          as.data.frame(matrix(nrow = max(clf$params$sparsity), ncol = 2))
        rownames(res.crossval$k$auc)      <-
          c(paste("k", c(1:max(
            clf$params$sparsity
          )), sep = "_"))
        colnames(res.crossval$k$auc)      <-
          c("empirical", "generalization")
        # add another table for accuracy
        res.crossval$k$acc                <- res.crossval$k$auc
        # recall
        res.crossval$k$rec                <- res.crossval$k$auc
        # precision
        res.crossval$k$prc                <- res.crossval$k$auc
        # f1-score
        res.crossval$k$f1s                <- res.crossval$k$auc
        # add another table for correlation
        res.crossval$k$cor                <- res.crossval$k$auc
        
        # for all k-sparse BEST models
        for (k in 1:length(res_train.digest$best_models))
        {
          k_sparse.name     <- names(res_train.digest$best_models)[k]
          mod               <-
            res_train.digest$best_models[[k_sparse.name]]
          # update the final indexes as the input X
          mod               <-
            updateObjectIndex(obj = mod, features = rownames(X))
          #mod.train         <- evaluateModel(mod = mod, X=x_train, y=y_train, clf=clf, eval.all = TRUE, force.re.evaluation = TRUE, mode='train')
          mod.train         <-
            mod # since this is the same as computed above
          mod.test          <-
            evaluateModel(
              mod = mod,
              X = x_test,
              y = y_test,
              clf = clf,
              eval.all = TRUE,
              force.re.evaluation = TRUE,
              mode = 'test'
            )
          
          if (!is.null(mod.train) & !is.null(mod.test))
          {
            # Empirical fitting score
            res.crossval$scores$empirical.auc[k_sparse.name, i]            <-
              mod.train$auc_
            res.crossval$scores$empirical.acc[k_sparse.name, i]            <-
              mod.train$accuracy_
            res.crossval$scores$empirical.rec[k_sparse.name, i]            <-
              mod.train$recall_
            res.crossval$scores$empirical.prc[k_sparse.name, i]            <-
              mod.train$precision_
            res.crossval$scores$empirical.f1s[k_sparse.name, i]            <-
              mod.train$f1_
            res.crossval$scores$empirical.cor[k_sparse.name, i]            <-
              mod.train$cor_
            
            # Generalization fitting score
            res.crossval$scores$generalization.auc[k_sparse.name, i]       <-
              mod.test$auc_
            res.crossval$scores$generalization.acc[k_sparse.name, i]       <-
              mod.test$accuracy_
            res.crossval$scores$generalization.rec[k_sparse.name, i]       <-
              mod.test$recall_
            res.crossval$scores$generalization.prc[k_sparse.name, i]       <-
              mod.test$precision_
            res.crossval$scores$generalization.f1s[k_sparse.name, i]       <-
              mod.test$f1_
            res.crossval$scores$generalization.cor[k_sparse.name, i]       <-
              mod.test$cor_
            
            # store by k
            # AUC
            res.crossval$k$auc[k_sparse.name, "empirical"]                 <-
              mod.train$auc_
            res.crossval$k$auc[k_sparse.name, "generalization"]            <-
              mod.test$auc_
            # Accuracy
            res.crossval$k$acc[k_sparse.name, "empirical"]                 <-
              mod.train$accuracy_
            res.crossval$k$acc[k_sparse.name, "generalization"]            <-
              mod.test$accuracy_
            # Recall
            res.crossval$k$rec[k_sparse.name, "empirical"]                 <-
              mod.train$recall_
            res.crossval$k$rec[k_sparse.name, "generalization"]            <-
              mod.test$recall_
            # Precision
            res.crossval$k$prc[k_sparse.name, "empirical"]                 <-
              mod.train$precision_
            res.crossval$k$prc[k_sparse.name, "generalization"]            <-
              mod.test$precision_
            # F1-Score
            res.crossval$k$f1s[k_sparse.name, "empirical"]                 <-
              mod.train$f1_
            res.crossval$k$f1s[k_sparse.name, "generalization"]            <-
              mod.test$f1_
            
            # Regression
            res.crossval$k$cor[k_sparse.name, "empirical"]                 <-
              mod.train$cor_
            res.crossval$k$cor[k_sparse.name, "generalization"]            <-
              mod.test$cor_
            
          } # if training and testing results exist
        } # end of k_sparse loop
        
        # if saving move one level up
        if (!(clf$params$popSaveFile == "NULL"))
        {
          setwd("..")
        }
      } # end if null digest
      
      # Compute FEATURE IMPORTANCE for each classifier in GENERALIZATION
      if (clf$params$compute.importance & !isLearnerSota(clf))
      {
        # we compute the feature importance and BTR languages and algorithms
        if (isClf(res.all[[i]]))
        {
          # if results are valid
          if (!is.null(res.all[[i]]$fip))
          {
            # if object exist
            if (!is.null(res.all[[i]]$fip$mda))
            {
              mda.all[res.all[[i]]$fip$feat.catalogue, i] <- res.all[[i]]$fip$mda
            }
            # if object exist
            if (!is.null(res.all[[i]]$fip$sda))
            {
              sda.all[res.all[[i]]$fip$feat.catalogue, i] <- res.all[[i]]$fip$sda
            }
            # if object exist
            if (!is.null(res.all[[i]]$fip$pda))
            {
              pda.all[res.all[[i]]$fip$feat.catalogue, i] <- res.all[[i]]$fip$pda
            }
          }
        } else
        {
          # print out information (might be errors)
          print(res.all[[i]])
          next
        }
      } # end importance
      
      # Do we need to return everything back ?
      if (return.all)
      {
        res.crossval$nfold[[i]]             <- list(results = res_train,
                                                    resultsDigest = res_train.digest)
      }
      
    } # end for (dispatching results)
    
    if (clf$params$verbose)
      cat("... All results from cross validation are dispatched\n")
    
    # reorder results function
    reorderByRownamesNumeric <- function(mat)
    {
      ind <- as.numeric(gsub("k_", "", rownames(mat)))
      mat <- mat[order(ind), ]
    }
    
    # reorder results
    res.crossval$scores$empirical.auc       <-
      reorderByRownamesNumeric(res.crossval$scores$empirical.auc)
    res.crossval$scores$empirical.acc       <-
      reorderByRownamesNumeric(res.crossval$scores$empirical.acc)
    res.crossval$scores$empirical.rec       <-
      reorderByRownamesNumeric(res.crossval$scores$empirical.rec)
    res.crossval$scores$empirical.prc       <-
      reorderByRownamesNumeric(res.crossval$scores$empirical.prc)
    res.crossval$scores$empirical.f1s       <-
      reorderByRownamesNumeric(res.crossval$scores$empirical.f1s)
    res.crossval$scores$empirical.cor       <-
      reorderByRownamesNumeric(res.crossval$scores$empirical.cor)
    res.crossval$scores$generalization.auc  <-
      reorderByRownamesNumeric(res.crossval$scores$generalization.auc)
    res.crossval$scores$generalization.acc  <-
      reorderByRownamesNumeric(res.crossval$scores$generalization.acc)
    res.crossval$scores$generalization.rec  <-
      reorderByRownamesNumeric(res.crossval$scores$generalization.rec)
    res.crossval$scores$generalization.prc  <-
      reorderByRownamesNumeric(res.crossval$scores$generalization.prc)
    res.crossval$scores$generalization.f1s  <-
      reorderByRownamesNumeric(res.crossval$scores$generalization.f1s)
    res.crossval$scores$generalization.cor  <-
      reorderByRownamesNumeric(res.crossval$scores$generalization.cor)
    
    
    # auc
    if (!is.null(dim(res.crossval$scores$empirical.auc)))
    {
      res.crossval$scores$mean.auc            <-
        data.frame(cbind(
          rowMeans(res.crossval$scores$empirical.auc, na.rm = TRUE),
          rowMeans(res.crossval$scores$generalization.auc, na.rm = TRUE)
        ))
      colnames(res.crossval$scores$mean.auc)  <-
        c("empirical", "generalization")
    }
    
    # accuracy
    if (!is.null(dim(res.crossval$scores$empirical.acc)))
    {
      res.crossval$scores$mean.acc            <-
        data.frame(cbind(
          rowMeans(res.crossval$scores$empirical.acc, na.rm = TRUE),
          rowMeans(res.crossval$scores$generalization.acc, na.rm = TRUE)
        ))
      colnames(res.crossval$scores$mean.acc)  <-
        c("empirical", "generalization")
    }
    
    # recall
    if (!is.null(dim(res.crossval$scores$empirical.rec)))
    {
      res.crossval$scores$mean.rec            <-
        data.frame(cbind(
          rowMeans(res.crossval$scores$empirical.rec, na.rm = TRUE),
          rowMeans(res.crossval$scores$generalization.rec, na.rm = TRUE)
        ))
      colnames(res.crossval$scores$mean.rec)  <-
        c("empirical", "generalization")
    }
    
    # precision
    if (!is.null(dim(res.crossval$scores$empirical.prc)))
    {
      res.crossval$scores$mean.prc            <-
        data.frame(cbind(
          rowMeans(res.crossval$scores$empirical.prc, na.rm = TRUE),
          rowMeans(res.crossval$scores$generalization.prc, na.rm = TRUE)
        ))
      colnames(res.crossval$scores$mean.prc)  <-
        c("empirical", "generalization")
      
      
    }
    
    # f1-score
    if (!is.null(dim(res.crossval$scores$empirical.f1s)))
    {
      res.crossval$scores$mean.f1s            <-
        data.frame(cbind(
          rowMeans(res.crossval$scores$empirical.f1s, na.rm = TRUE),
          rowMeans(res.crossval$scores$generalization.f1s, na.rm = TRUE)
        ))
      colnames(res.crossval$scores$mean.f1s)  <-
        c("empirical", "generalization")
    }
    
    # correlation
    if (!is.null(dim(res.crossval$scores$empirical.cor)))
    {
      res.crossval$scores$mean.cor            <-
        data.frame(cbind(
          rowMeans(res.crossval$scores$empirical.cor, na.rm = TRUE),
          rowMeans(res.crossval$scores$generalization.cor, na.rm = TRUE)
        ))
      colnames(res.crossval$scores$mean.cor)  <-
        c("empirical", "generalization")
    }
    
    # adding results for feature importance
    if (clf$params$compute.importance & !isLearnerSota(clf))
    {
      res.crossval$fip <- list(
        mda = mda.all,
        sda = sda.all,
        pda = pda.all,
        fpf = rowSums(!is.na(mda.all)) # feature prevalenc in folds
      )
    }
    
    return(res.crossval)
  }


# ################################################################
# # COMPUTING PREDICTIONS
# ################################################################
#
# #' Computes the predected classification using a given model
# #'
# #' @description This function computes the predicted classification
# #' @export
# #' @param X: dataset to classify
# #' @param clf: an object containing the different parameters of the classifier
# #' @param mod: a model object to be used in the class prediction
# #' @param intercept: not needed for the ter family
# #' @param lev: the class level to be used for the new class
# #' @return a vector with the predicted classification of the samples
# prediction <- function(X, clf, mod=NULL, intercept=NULL, lev = c(-1,1))
# {
#   # if we don't send a specific model to test we'll use the default best model of the learner
#   if(is.null(mod)){
#     # Make sure we have a model for terLearner family
#     if(clf$learner != "sota.svm" & clf$learner != "sota.rf")
#     {
#       # if the fitting process has not been launched
#       if(is.null(clf$models))
#       {
#         stop("Please provide a model to use in the prediction process")
#       }else
#       {
#         mod <- getTheBestModel(clf)
#       }
#       # # Setting the intercept
#       # if(is.null(intercept))
#       # {
#       #   stop("Setting the default intercept = 0!")
#       # }
#     }
#   }
#
#   # Setting the intercept
#   if(is.null(intercept))
#   {
#     if(!is.null(mod)) # mod exists
#     {
#       if(!is.null(mod$intercept_)) # mod contains an intercept
#       {
#         intercept <- mod$intercept_
#       }else # mod does not contain an intercept
#       {
#         stop("Setting the default intercept = 0!")
#       }
#     }else # mod does not exist and intercept
#     {
#       intercept <- 0
#       warning("Setting intercept to default 0!")
#     }
#     if(clf$learner=="sota.rf") # output of random forest is posterior probability of class +1
#     {
#       intercept = 0.5
#     }
#   }
#
#   switch(clf$learner,
#          terda=
#          {
#            # case 'terda' here...
#            score <- computeModelScore(X = X, mod = mod, clf) # compute the score of the model
#            yhat <- lev[factor(score - intercept > 0)]
#          },
#          terga1=
#          {
#            # case 'terga_v1' here...
#            score <- computeModelScore(X = X, mod = mod, clf) # compute the score of the model
#            yhat <- lev[factor(score - intercept > 0)]
#          },
#          terga2=
#          {
#            # case 'terga' here...
#            score <- computeModelScore(X = X, mod = mod, clf) # compute the score of the model
#            yhat <- lev[factor(score - intercept > 0)]
#          },
#          terBeam=
#          {
#            # case 'terbeam' here...
#            score <- computeModelScore(X = X, mod = mod, clf) # compute the score of the model
#            yhat <- lev[factor(score - intercept > 0)]
#          },
#          sota.svm=
#          {
#            # case 'stoa.svm' here...
#            model <- clf$models[[clf$best.model.id]]
#            X.reduced <- t(X[model$selected.features,])
#            score <- as.numeric(predict(model$obj, X.reduced, type="decision"))
#            # yhat <- rep(1,ncol(X))
#            # yhat[score < intercept] <- -1
#            yhat <- lev[factor(score - intercept > 0)]
#          },
#          sota.rf=
#          {
#            # case 'stoa.rf' here...
#            model <- clf$models[[clf$best.model.id]]
#            X.reduced <- t(X[model$selected.features,])
#            score <- predict(model$obj,X.reduced,type="prob")[,"1"]
#            yhat <- lev[factor(score - intercept > 0)]
#          },
#          {
#            stop('This method does not exist !')
#          })
#         names(yhat) = colnames(X)
#         names(score) = colnames(X)
#   return(list(yhat=yhat, score=score))
# }
