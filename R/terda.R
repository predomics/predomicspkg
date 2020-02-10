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
# @script: terda.R                                          
# @author: Yann Chevaleyre
# @author: Lucas Robin
# @date: August 2016                                                    
################################################################


#' terda: terda classifier parameter function
#'
#' @title terda
#' @description terbeam is a model search algorithm.
#' @param sparsity: number of features in a given model. This is a vector with multiple lengths.
#' @param nIterations: ??
#' @param max.nb.features: create the glmnet object using only the top most significant features (default:1000)
#' @param kBest: ??
#' @param method: ??
#' @param kStep: ??
#' @param vartype: (??)
#' @param gamma: ??
#' @param nRR: (??) (default:FALSE)
#' @param lb: ??
#' @param ub: ??
#' @param language is the language that is used by the different algorithms {bin, bininter, ter, terinter, ratio}, (default:"terinter")
#' @param scoreFormula: a Function that contains the ratio Formula or other specific ones
#' @param epsilon: a small value to be used with the ratio language (useCustomLanguage) (default: NULL). When null it is going to be calculated by the minimum value of X divided by 10.
#' @param objective: this can be auc, cor or aic. Terga can also predict regression, other than class prediction. (default:auc)
#' @param evalToFit: The model performance attribute to use as fitting score (default:"fit_"). Other choices are c("auc_","accuracy_","precision_","recall_","f_score_")
#' @param k_penalty: Penalization of the fit by the k_sparsity (default: 0)
#' @param intercept: (??) (default:NULL)
#' @param popSaveFile: (??)
#' @param final.pop.perc: ??
#' @param plot: Plot different graphics (default:FALSE).
#' @param verbose: print out information on the progress of the algorithm (default:TRUE)
#' @param warnings: Print out warnings when runnig (default:FALSE).
#' @param debug: print out debug infotmation when activated (default: FALSE)
#' @param print_ind_method: One of c("short","graphical") indicates how to print a model and subsequently a population during the run (default:"short").
#' @param parallelize.folds: parallelize folds when cross-validating (default:TRUE)
#' @param nCores: the number of cores to execute the program. If nCores=1 than the program runs in a non parallel mode
#' @param seed: the seed to be used for reproductibility. If seed=NULL than it is not taken into account (default:NULL).
#### TODO check
#' @param experiment.id: The id of the experiment that is to be used in the plots and comparitive analyses (default is the learner's name, when not specified)
#' @param experiment.description: A longer description of the experiment. This is important when many experiments are run and can also be printed in by the printExperiment function.
#' @param experiment.save: Data from an experiment can be saved with different levels of completness, with options to be selected from c("nothing", "minimal", "full"), default is "minimal"
#' @return an object containing a list of parameters for this classifier
#' @export
terda <- function(sparsity = 5, nIterations = 5, max.nb.features = 1000, kBest = "NULL",
                  # terda core options
                  method = "glmnetRR", kStep = "NULL", vartype = "real", gamma = 0.7, nRR = 1, lb = -1, ub = 1,
                  # language in {bin, bininter, ter, terinter}
                  language = "terinter",
                  scoreFormula = scoreRatio, epsilon = "NULL",
                  # glmnet parameters
                  nblambdas = 1000, #number of lambdas in the glmnet model
                  # evaluation options
                  objective = "auc", evalToFit = "auc_", k_penalty=0, intercept = "NULL",
                  # population options
                  popSaveFile = "NULL", final.pop.perc = 100, alpha = 0.5,
                  # output options
                  plot = FALSE, verbose = TRUE, warnings = FALSE, debug = FALSE, print_ind_method = "short", parallelize.folds = TRUE,
                  # computing options
                  nCores = 4, seed = "NULL", #maxTime = Inf,
                  # experiment options
                  experiment.id = "NULL", experiment.description = "NULL", experiment.save = "nothing"
)
{
  clf <- list()
  clf$learner                   <- "terda"
  clf$params                    <- list()
  clf$experiment                <- list()         # information about the experiment
  clf$params$method             <- method
  clf$params$objective          <- objective
  clf$params$kStep              <- kStep
  clf$params$kBest              <- kBest
  clf$params$gamma              <- gamma
  clf$params$sparsity           <- sparsity
  clf$params$nRR                <- nRR
  clf$params$lb                 <- lb
  clf$params$ub                 <- ub
  
  # glmnet model parameters
  clf$params$nblambdas          <- nblambdas
  
  clf$params$vartype            <- vartype
  clf$params$alpha              <- alpha          # glmnet's elastcinet parameter
  clf$params$nIterations        <- nIterations
  clf$params$max.nb.features    <- max.nb.features
  
  # print out intermediary results
  clf$params$plot               <- plot           # save plots duting the process
  clf$params$verbose            <- verbose        # print out logs.
  clf$params$warnings           <- warnings       # print out warnings
  clf$params$debug              <- debug          # print out logs.
  clf$params$print_ind_method   <- print_ind_method # method to print individual
  
  # Computing options
  clf$params$nCores             <- nCores         # parallel computing
  clf$params$parallel           <- (nCores > 1)   # parallel computing
  clf$params$parallelize.folds  <- parallelize.folds
  clf$params$parallel.local     <- FALSE
  clf$params$seed               <- seed           # fix the seed to be able to reproduce results
  
  clf$params$evalToFit          <- evalToFit
  clf$params$k_penalty          <- k_penalty
  clf$params$final.pop.perc     <- final.pop.perc
  clf$params$popSaveFile        <- popSaveFile
  
  clf$params$language           <- language
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



terda_fit <- function(X, y, clf) 
{
  
  if(clf$params$objective == "cor")
  {
    cat("... ... terda: setting regression mode.")
    clf$params$language <- "regression"
  }

  # Setting the language environment
  switch(clf$params$language,
         regression=
         {
           # unconstrained logistic regression
           if(clf$params$verbose){ print("Setting environment for regression not implemented for the moment returning null") }
           
           # if(clf$params$evalToFit != "auc_" & clf$params$evalToFit != "accuracy_")
           # {
           #   clf$params$evalToFit <- "auc_"
           #   warning("terga1_fit: not a valid evalToFit. Changing to auc_.")
           # }
           return(NULL)
         },
         logreg=
         {
           # unconstrained logistic regression
           if(clf$params$verbose){ print("Setting environment for unconstrained logistic regression") }
           if(clf$params$evalToFit != "auc_" & clf$params$evalToFit != "accuracy_")
           {
             clf$params$evalToFit <- "auc_"
             warning("terga1_fit: not a valid evalToFit. Changing to auc_.")
           }
         },
         ter= 
         {
           # ternary language without intercept (maximize the accuracy)
           if(clf$params$verbose){print("Setting environment for the language 'ter'")}
           # note that here with the ter language we could not use auc to fit since the intercept should be 0
           if(clf$params$evalToFit != "auc_" & clf$params$evalToFit != "accuracy_")
           {
             clf$params$intercept = 0
             clf$params$evalToFit <- "accuracy_"
             warning("terga1_fit: not a valid evalToFit. Changing to accuracy_.")
           }else if(clf$params$evalToFit == "auc_")
           {
             clf$params$intercept = 0
             clf$params$evalToFit <- "accuracy_"
             warning("terga1_fit: changing evalToFit from auc_ to accuracy_ because of the language.")
           }
         },
         terinter=
         {
           # ternary language with intercept (maximize the accuracy)
           if(clf$params$verbose){print("Setting environment for the language 'terinter'")}
           if(clf$params$evalToFit != "auc_" & clf$params$evalToFit != "accuracy_")
           {
             clf$params$evalToFit <- "auc_"
             warning("terga1_fit: not a valid evalToFit. Changing to auc_.")
           }
         },
         bin=
         {
           # ternary language without intercept (maximize the auc)
           if(clf$params$verbose){print("Setting environment for the language 'bin'")}
           # note that here with the bin language we could not use auc to fit since the intercept should be 0
           if(clf$params$evalToFit != "auc_" & clf$params$evalToFit != "accuracy_")
           {
             clf$params$intercept = 0
             clf$params$evalToFit <- "accuracy_"
             warning("terga1_fit: not a valid evalToFit. Changing to accuracy_.")
           }else if(clf$params$evalToFit == "auc_")
           {
             clf$params$intercept = 0
             clf$params$evalToFit <- "accuracy_"
             warning("terga1_fit: changing evalToFit from auc_ to accuracy_ because of the language.")
           }
         },
         bininter=
         {
           # ternary language without intercept (maximize the auc)
           if(clf$params$verbose){print("Setting environment for the language 'bininter'")}
           if(clf$params$evalToFit != "auc_" & clf$params$evalToFit != "accuracy_")
           {
             clf$params$evalToFit <- "auc_"
             warning("terga1_fit: not a valid evalToFit. Changing to auc_.")
           }
         },
         ratio=
         {
           # ternary language without intercept (maximize the auc)
           if(clf$params$verbose){print("Setting environment for the language 'ratio'")}
           if(clf$params$evalToFit != "auc_" & clf$params$evalToFit != "accuracy_")
           {
             clf$params$evalToFit <- "auc_"
             warning("terga1_fit: not a valid evalToFit. Changing to auc_.")
           }
         },
         {
           stop(paste("The language",clf$params$language, "is not implemented !"))
         }
  )

  # Print the experiment configuration
  if(clf$params$verbose) printClassifier(obj = clf)
  
  # if(clf$params$max.nb.features < nrow(X))
  # {
  #   if(clf$params$verbose) print(paste("... ... restricting X to the",clf$params$max.nb.features,"most significant variables"))  
  #   selected.features <- rownames(clf$feature.cor[order(clf$feature.cor$p),][1:min(clf$params$max.nb.features, nrow(X)),])
  #   X <- X[selected.features,]
  #   clf$coeffs_ <- clf$coeffs_[selected.features]
  #   clf$feature.cor <- clf$feature.cor[selected.features,]
  # }

  if(clf$params$verbose) print("... ... X needs to be transposed")
  tX = as.data.frame(t(X))
  if(clf$params$verbose) print("... ... X transposed")

  # setting kBest
  if(clf$params$kBest=="NULL")
  {
    clf$params$kBest <- round(nrow(X)/2)
  }
  # setting kStep
  if (clf$params$kStep=="NULL")
  {
    clf$params$kStep = round(ncol(tX)/10)
  } 
  
  if(clf$params$method == "custom") # algo vraiment basique, je résous le LP et je round
  { 
    if(clf$params$verbose) print("... ... custom method (not stable)")
    l <- buildlp(tX, y, 
                 k.sparse = clf$params$sparsity, 
                 vartype = "real", 
                 gamma = clf$params$gamma, 
                 lb = clf$params$lb, 
                 ub = clf$params$ub)
    
    wb <- runsolver(l, tX)
    # remove intercept, which will be computed again latter for better precision
    wb <- wb[-length(wb)]
    res <- multipleRR(clf, X, y, wb, clf$params$nRR)    
    
  } else if(clf$params$method == "myRRSkbest") 
  {
    if(clf$params$verbose) print("... ... myRRSkbest method (not stable)")
    res <- myRRSkbest(clf, X, y)
  } else if (clf$params$method == "glmnetRR") 
  {
    if(clf$params$verbose) print("... ... glmnetRR method")
    res <- glmnetRR(clf, X, y)
  }

  return(res)
}

