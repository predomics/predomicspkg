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
# @script: sota.svm.R                                          
# @author: Edi Prifti
# @author: Blaise Hanczar
# @date: September 2016
# @date: December 2016 (a working integrated version with predomics)
################################################################


#' sota.svm: launching svm classifier
#'
#' @title sota.svm
#' @importFrom kernlab ksvm
#' @importFrom kernlab predict
#' @description sota.svm is a wrapper that executes svm using the same framework as for the predomics package.
#' @param sparsity: number of features in a given model. This is a vector with multiple lengths.
#' @param objective: prediction mode (default: auc)
#' @param max.nb.features: create the glmnet object using only the top most significant features (default:1000)
#' @param language is the language that is used by the different algorithms {bin, bininter, ter, terinter, ratio}, (default:"sota")
#' @param intercept: (Interceot for the a given model) (default:NULL)
#' @param evalToFit: Which model property will be used to select the best model among different k_sparsities (default: auc_)
#' @param k_penalty: Penalization of the fit by the k_sparsity (default: 0)
#' @param scaled: ??
#' @param type: ??
#' @param kernel: ??
#' @param kpar: ??
#' @param C: (??)
#' @param nu: ??
#' @param epsilon.hp: (??) (for the SVM)
#' @param prob.model: ??
#' @param class.weights: ??
#' @param fit: ??
#' @param cache: (??)
#' @param tol: ??
#' @param shrinking: ??
#' @param na.action: ??
#' @param popSaveFile: (??)
#' @param seed: the seed to be used for reproductibility. If seed=NULL than it is not taken into account (default:NULL).
#' @param nCores: the number of CPUs to run the program in parallel
#' @param plot: Plot graphics indicating the evolution of the simulation (default:FALSE)
#' @param verbose: print out information on the progress of the algorithm (default:TRUE)
#' @param warnings: Print out warnings when runnig (default:FALSE).
#' @param debug: print out information on the progress of the algorithm (default:FALSE)
#' @param print_ind_method: One of c("short","graphical") indicates how to print a model and subsequently a population during the run (default:"short").
#' @param experiment.id: The id of the experiment that is to be used in the plots and comparitive analyses (default is the learner's name, when not specified)
#' @param experiment.description: A longer description of the experiment. This is important when many experiments are run and can also be printed in by the printExperiment function.
#' @param experiment.save: Data from an experiment can be saved with different levels of completness, with options to be selected from c("nothing", "minimal", "full"), default is "minimal"
#' @return an object containing a list of parameters for this classifier
#' @export 
sota.svm <- function(sparsity = c(1:30), # when sparsity == 0 it means that we can not fix it.
                     objective = "auc",
                     max.nb.features = 1000,
                     intercept = 0, # ksvm incorporates the threshold in the score.
                     language = "svm",
                     evalToFit = "auc_",
                     k_penalty=0, 
                     scaled = TRUE,
                     type = NULL,
                     kernel ="rbfdot",
                     kpar = "automatic",
                     C = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000),
                     nu = 0.2,
                     epsilon.hp = 0.1, # hyper parameter
                     prob.model = FALSE,
                     class.weights = NULL,
                     # cross = 0,  done by the fit function of predomics
                     fit = TRUE,
                     cache = 40,
                     tol = 0.001,
                     shrinking = TRUE,
                     na.action = na.omit,
                     popSaveFile = "NULL",
                     seed = "NULL",
                     nCores = 4,
                     verbose = TRUE,
                     plot = FALSE,
                     warnings = FALSE,
                     debug = FALSE, 
                     print_ind_method = "short", 
                     experiment.id = NULL, 
                     experiment.description = NULL, 
                     experiment.save = "nothing") 
{
  
  clf <- list() # create a classifier object
  clf$learner                   <- "sota.svm" # name of the method
  clf$params                    <- list() # parameter list
  clf$params$scoreFormula       <- "NULL"
  clf$params$parallel           <- FALSE
  clf$params$sparsity           <- sparsity # number of non zero variables in the model
  clf$params$scaled             <- scaled
  clf$params$objective          <- objective
  clf$params$max.nb.features    <- max.nb.features
  clf$params$language           <- language
  clf$params$intercept          <- intercept
  clf$params$type               <- type # A logical vector indicating the variables to be scaled. If scaled is of length 1, the value is recycled as many times as needed and all non-binary variables are scaled. Per default, data are scaled internally (both x and y variables) to zero mean and unit variance. The center and scale values are returned and used for later predictions.
  #   C-svc C classification
  #   nu-svc nu classification
  #   C-bsvc bound-constraint svm classification
  #   spoc-svc Crammer, Singer native multi-class
  #   kbb-svc Weston, Watkins native multi-class
  #   one-svc novelty detection
  #   eps-svr epsilon regression
  #   nu-svr nu regression
  #   eps-bsvr bound-constraint svm regression
  clf$params$kernel             <- kernel
  clf$params$kpar               <- kpar
  clf$params$C                  <- C # cost of constraints violation (default: 1) this is the ‘C’-constant of the regularization term in the Lagrange formulation
  clf$params$nu                 <- nu # parameter needed for nu-svc, one-svc, and nu-svr. The nu parameter sets the upper bound on the training error and the lower bound on the fraction of data points to become Support Vectors (default: 0.2).
  clf$params$epsilon.hp         <- epsilon.hp # epsilon in the insensitive-loss function used for eps-svr, nu-svr and eps-bsvm (default: 0.1)
  clf$params$prob.model         <- prob.model # if set to TRUE builds a model for calculating class probabilities or in case of regression, calculates the scaling parameter of the Laplacian distribution fitted on the residuals. Fitting is done on output data created by performing a 3-fold cross-validation on the training data. For details see references. (default: FALSE)
  clf$params$class.weights      <- class.weights # a named vector of weights for the different classes, used for asymmetric class sizes. Not all factor levels have to be supplied (default weight: 1). All components have to be named.
  #clf$params$cross              <- cross # if a integer value k>0 is specified, a k-fold cross validation on the training data is performed to assess the quality of the model: the accuracy rate for classification and the Mean Squared Error for regression
  clf$params$fit                <- fit # indicates whether the fitted values should be computed and included in the model or not (default: TRUE)
  clf$params$cache              <- cache # cache memory in MB (default 40)
  clf$params$tol                <- tol # tolerance of termination criterion (default: 0.001)
  clf$params$shrinking          <- shrinking # option whether to use the shrinking-heuristics (default: TRUE)
  clf$params$na.action          <- na.action # A function to specify the action to be taken if NAs are found. The default action is na.omit, which leads to rejection of cases with missing values on any required variable. An alternative is na.fail, which causes an error if NA cases are found. (NOTE: If given, this argument must be named.)
  clf$params$popSaveFile        <- popSaveFile
  
  clf$params$seed               <- seed
  clf$params$nCores             <- nCores         # parallel computing
  clf$params$parallel           <- nCores > 1     # parallel computing
  clf$params$plot               <- plot           # plot out logs.
  clf$params$verbose            <- verbose        # print out graphics
  clf$params$warnings           <- warnings       # print out warnings
  clf$params$debug              <- debug          # print out logs.
  clf$params$print_ind_method   <- print_ind_method # method to print individual
  
  clf$params$evalToFit          <- evalToFit
  clf$params$k_penalty          <- k_penalty
  
  # Experiment information
  if(!is.null(experiment.id))
  {
    clf$experiment$id          <- experiment.id
  }else
  {
    clf$experiment$id          <- clf$learner
  }
  
  if(!is.null(experiment.description))
  {
    clf$experiment$description <- experiment.description
  } else 
  {
    clf$experiment$description <- paste(clf$learner, date() , sep = " ")
  }
  
  clf$experiment$save          <- experiment.save
  
  return(clf)
}


sota.svm_fit <- function(X, y, clf) 
{
  # transpose (in the predomics package we have observations in the columns)
  x <- as.matrix(t(X))
  
  model_collection <- list()
  for (i in 1:length(clf$params$sparsity)) # sparsity is = k, i.e. the number of features in a model
  { 
    k <- clf$params$sparsity[i]
    if(clf$params$verbose) cat("... ... Resolving problem with\t", k, "\tvariables ...\n")

    selected.features <- rownames(clf$feature.cor)[order(clf$feature.cor$p)][1:k]
    x.reduced <- x[,selected.features]
    
    # Optimization of the hyper-parameter C
    best.C <- NA
    if(length(clf$params$C)>1) # find the best C
    {
      optimize.C <- list() # unify the results
      if(clf$params$parallel) # If parallel computing
      {
        optimize.C <- foreach (j = 1:length(clf$params$C)) %dopar%
        {
          set.seed(clf$params$seed) # otherwise each time this will be different
          svm <- ksvm(x = x.reduced, y=y,
                      type = clf$params$type,
                      scaled = clf$params$scaled,
                      kernel = clf$params$kernel, 
                      kpar = clf$params$kpar,
                      C = clf$params$C[j], 
                      nu = clf$params$nu,
                      epsilon = clf$params$epsilon.hp, 
                      prob.model = clf$params$prob.model,
                      class.weights = clf$params$class.weights, 
                      cross = 5, 
                      fit = clf$params$fit, 
                      cache = clf$params$cache,
                      tol = clf$params$tol, 
                      shrinking = clf$params$shrinking,
                      na.action = clf$params$na.action
          )
          attr(svm,"cross")
        } # end foreach loop
        names(optimize.C) <- clf$params$C
      }else
      {
        for (j in 1:length(clf$params$C)) 
        {
          set.seed(clf$params$seed) # otherwise each time this will be different
          svm <- ksvm(x = x.reduced, y=y,
                      type = clf$params$type,
                      scaled = clf$params$scaled,
                      kernel = clf$params$kernel, 
                      kpar = clf$params$kpar,
                      C = clf$params$C[j], 
                      nu = clf$params$nu,
                      epsilon = clf$params$epsilon, 
                      prob.model = clf$params$prob.model,
                      class.weights = clf$params$class.weights, 
                      cross = 5, 
                      fit = clf$params$fit, 
                      cache = clf$params$cache,
                      tol = clf$params$tol, 
                      shrinking = clf$params$shrinking,
                      na.action = clf$params$na.action
          )
          optimize.C[[j]] <- c(attr(svm,"cross")) 
        } # end for loop
        names(optimize.C) <- clf$params$C
      } # end else //
      
      best.C <- clf$params$C[which.min(unlist(optimize.C))]
      
      if(clf$params$verbose)
      {
        print(paste("best.C", best.C))
      }
    }else # if C is fixed use it as the best.C
    {
      best.C <- clf$params$C
    }
    
    
    # Launch ksvm with the parameters from the clf
    set.seed(clf$params$seed)
    svm <- ksvm(x = x.reduced, y=y,
                type = clf$params$type,
                scaled = clf$params$scaled,
                kernel = clf$params$kernel, 
                kpar = clf$params$kpar,
                C = best.C, 
                nu = clf$params$nu,
                epsilon = clf$params$epsilon, 
                prob.model = clf$params$prob.model,
                class.weights = clf$params$class.weights, 
                cross = 0, 
                fit = clf$params$fit, 
                cache = clf$params$cache,
                tol = clf$params$tol, 
                shrinking = clf$params$shrinking,
                na.action = clf$params$na.action
    )
    
    # build the model for this given sparsity
    mod.res                  <- list() # to send the results
    # Fill the object up with the rest of the values
    mod.res$learner          <- clf$learner
    mod.res$language         <- clf$params$language
    mod.res$names_           <- selected.features
    # match the index
    mod.res$indices_         <- match(selected.features, rownames(X))
    mod.res$coeffs_          <- NA
    
    #mod.res$indices_         <- which(rownames(X) %in% selected.features)
    mod.res$eval.sparsity    <- length(unique(mod.res$names_))
    # Add the svm object
    mod.res$obj              <- svm
    mod.res$auc_             <- NA
    mod.res$cor_             <- NA
    mod.res$aic_             <- NA
    mod.res$score_           <- NA
    mod.res$fit_             <- NA
    mod.res$unpenalized_fit_ <- NA
    mod.res$intercept_       <- NA
    mod.res$sign_            <- NA
    # evaluate all
    mod.res                  <- computeCoeffSVMLin(X, y, clf=clf, mod=mod.res) # compute the coefficients for the linear SVM
    mod.res                  <- evaluateModel(mod = mod.res, 
                                              X = X, 
                                              y = y, 
                                              clf = clf, 
                                              eval.all = TRUE, 
                                              force.re.evaluation = TRUE, 
                                              estim.feat.importance = TRUE)
    
    if(clf$params$verbose)
    {
      try(printModel(mod = mod.res, method = clf$params$print_ind_method, score = "fit_"), silent = TRUE)
    }
    
    # Create a resulting model collection object
    model_collection[[i]]    <- list(mod.res)
    
  } # end of sparsity loop
  
  names(model_collection)    <- paste("k", clf$params$sparsity, sep="_")
  return(model_collection)
}
