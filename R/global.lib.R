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
# @script: global.lib.R                                          
# @author: Edi Prifti
# @author: Lucas Robin
# @author: Yann Chevaleyre
# @author: Jean-Daniel Zucker
# @date: August 2016    
# @date: November 2023                                         

################################################################
# CONTENTS
# ========= COMPUTING MODEL OBJECTS, SCORES, ERRORS...
# ========= GETTERS, SETTERS ...
# ========= CHECKERS
# ========= CONVERTERS
# ========= ADDERS, REMOVERS
# ========= COMPUTING PREVALENCE, ABUNDANCE & other POP ATTRIBUTES
# ========= CROSS VALIDATION
# ========= FILTERING PROCEDURES
# ========= NULL DISTRIBUTION, FUNCTIONS
# ========= DIGESTING RESULTS
# ========= EXPERIMENT HANDLERS
# ========= FILE MANIPULATION FUNCTIONS
################################################################



################################################################
# COMPUTING MODEL OBJECTS, SCORES, ERRORS...
################################################################


#' Evaluates the confusion Matrix of the predicted class and the class to predict
#'
#' @description This function evaluates the accuracy of a model
#' @param mod: a model object to be evaluated
#' @param X: dataset to classify
#' @param y: variable to predict
#' @param clf: an object containing the different parameters of the classifier
#' @return a confusion matrix
computeConfusionMatrix <- function(mod, X, y, clf)
{
  yhat <- evaluateYhat(mod = mod, X = X, y = y, clf = clf)
  if(is.null(yhat))
  {
    return(NULL)
  }
  # yhat <- factor(yhat, levels = names(table(clf$data$y)))
  
  if(length(yhat) != length(y))
  {
    return(NULL)
  }
  cm <- table(y, yhat, dnn = c("y", "yhat"))
  
  return(cm)
}


#' Compute other prediction scores such as precision, recall and f-score
#'
#' @description This function computes prediction scores based on the confusion matrix such as accuracy, precision, recall and f-score
#' @param mod: a model object to be evaluated
#' @param X: dataset to classify
#' @param y: variable to predict
#' @param clf: an object containing the different parameters of the classifier
#' @param mode: training or testing mode
#' @return a model whose evaluation parameters are updated
evaluateAdditionnalMetrics <- function(mod, X, y, clf, mode = "train")
{
  
  # if the following attributes are selected, then we need to fix it since they are derivates of a score
  if(clf$params$objective == "auc" & (clf$params$evalToFit != "fit_" | clf$params$evalToFit != "unpenalized_fit_"))
  {
    clf$params$evalToFit <- "accuracy_"
  }
  
  if(clf$params$evalToFit == "accuracy_") # additional metrics is auc_
  {
    # compute the auc
    aucg                 <- evaluateAUC(score = mod$score, y = y, sign = ">")
    aucl                 <- evaluateAUC(score = mod$score, y = y, sign = "<")
    mod$auc_             <- max(aucg, aucl)
  }
  
  if(clf$params$evalToFit == "auc_") # additional metrics is auc_
  {
    # evalute the accuracy that is not measured
    mod <- evaluateAccuracy(mod = mod, X = X, y = y, clf = clf, force.re.evaluation = TRUE, mode = mode)
  }
  
  if(!myAssertNotNullNorNa(mod$confusionMatrix_)) 
  {
    # visit this website for more information on the measures https://en.wikipedia.org/wiki/Precision_and_recall
    mod$confusionMatrix_ <- computeConfusionMatrix(mod, X, y, clf)
    
    if(is.null(mod$confusionMatrix_))
    {
      return(NULL)
    }
  }
  
  # precision = tp/(tp+fp)
  # cm = confusion matrix (2,2 is the positive class)
  mod$precision_ <- mod$confusionMatrix[2, 2] / (mod$confusionMatrix[2, 2] + mod$confusionMatrix[2, 1])
  #mod$precision_ <- mod$confusionMatrix[1, 1] / (mod$confusionMatrix[1, 1] + mod$confusionMatrix[2, 1])
  
  # recall = tp/(tp+fn), aka sensitivity
  #mod$recall_    <- mod$confusionMatrix[1, 1] / (mod$confusionMatrix[1, 1] + mod$confusionMatrix[1, 2])
  mod$recall_    <- mod$confusionMatrix[2, 2] / (mod$confusionMatrix[2, 2] + mod$confusionMatrix[1, 2])
  
  mod$f1_       <- 2 * (mod$precision_ * mod$recall_) / (mod$precision_ + mod$recall_)
  
  return(mod)
}


#' Computes the predected classification using a given model
#'
#' @description This function evaluates the predicted classification either using (1) a model object that contains intercept and sign or (2) directly the attributes score, intercept, sign
#' @param mod: a model object to be used in the class prediction
#' @param X: dataset to classify
#' @param y: variable to predict
#' @param clf: an object containing the different parameters of the classifier
#' @param score: the score passed directly
#' @param intercept: the intercept passed directly
#' @param sign: the sign passed directly
#' @return a vector with the predicted classification of the samples
evaluateYhat <- function(mod = NULL, X, y, clf, score=NULL, intercept=NULL, sign=NULL)
{
  if(!isModel(mod))
  {
    stop("evaluateYhat: please provide a valid model object BTR or SOTA.")
  }
  
  if(isModelSotaRF(mod))
  {
    yhat <- predict(object = mod$obj, t(X[mod$indices_,]))
  }
  
  if(isModelSotaSVM(mod))
  {
    yhat <- predict(object = mod$obj, t(X[mod$indices_,]))
  }
  
  if(isModelSotaGLMNET(mod))
  {
    yhat <- predict(object = mod$obj, s = mod$lambda, newx = t(X), type="class")[,1]
  }
  
  if(!isModelSota(mod)) # if BTR
  {
    # score_
    if(!myAssertNotNullNorNa(mod$score_))
    {
      scorelist <- getModelScore(mod = mod, X = X, clf = clf, force.re.evaluation = FALSE) # compute the score of the model
      score <- scorelist$score_
      
      if(is.null(score)) 
      {
        return(NULL)
      }
      
    } else
    {
      score       <- mod$score_
    }
    
    if(!myAssertNotNullNorNa(mod$intercept_) | !myAssertNotNullNorNa(mod$sign_)) 
    {
      mod <- evaluateIntercept(mod = mod, X = X, y = y, clf = clf)
      intercept <- mod$intercept_
      sign <- mod$sign_
    }else{
      intercept <- mod$intercept_
      sign <- mod$sign_
    }
    
    if(!myAssertNotNullNorNa(score, "missing score from evaluateYhat", stop = FALSE)) return(NULL)
    if(!myAssertNotNullNorNa(intercept, "missing intercept from evaluateYhat", stop = FALSE)) return(NULL)
    if(!myAssertNotNullNorNa(sign, "missing sign from evaluateYhat", stop = FALSE)) return(NULL)
    
    lev   <- levels(as.factor(clf$data$y))
    
    # NOTE: in the score we may have infite values that come for instance from the ratio language
    # This means that whatever the intercept these examples are positive ones. As such we can omit
    # them when computing the intercept.
    ind.infinite  <- is.infinite(mod$score_)
    ind.nan       <- is.nan(mod$score_)
    
    # For the ratio language
    # CASE 1  
    # a/b > teta
    # a > teta * b
    # if b is infinite than whatever the teta the inequation is true
    # CASE 2
    # a/b < teta
    # a < teta * b
    # if b is infinite than whatever the teta the inequation is false
    
    if(mod$sign==">")
    {
      yhat.bool <- (score - intercept > 0)
      yhat.bool[ind.infinite] <- TRUE
      yhat.bool[ind.nan] <- FALSE # since 0 is not > than 0
      # compute the class
      yhat  <- factor(lev[as.factor(yhat.bool)], levels=lev)
      
    }else
    {
      yhat.bool <- (score - intercept < 0)
      yhat.bool[ind.infinite] <- FALSE
      yhat.bool[ind.nan] <- FALSE # since 0 is not < than 0
      # compute the class
      yhat  <- factor(lev[as.factor(yhat.bool)], levels=lev)
    }
  }
  
  return(yhat)
}


#' Compute other prediction scores such as precision, recall and f-score
#'
#' @description This function computes prediction scores based on the confusion matrix such as accuracy, precision, recall and f-score
#' @param X: dataset to classify
#' @param y: variable to predict
#' @param mod: a predomics object to be updated
#' @param clf: an object containing the different parameters of the classifier
#' @return a model whose evaluation parameters are updated or a list containing coefficients and intercept if mod is not set.
#' @export
computeCoeffSVMLin <- function(X, y, clf=NULL, mod=NULL)
{
  if(is.null(mod))
  {
    print("Compute coeffs from a dataset directly not a predomics object")
    p = nrow(X)
    svm <- ksvm(as.matrix(t(X)), y, type="C-svc", kernel='vanilladot',C=1)
  }else
  {
    if(!isModel(obj = mod))
    {
      stop("computeCoeffSVMLin: is not a valid predomics model.")
    }
    if(!isModelSotaSVM(obj = mod))
    {
      stop("computeCoeffSVMLin: is a valid predomics model, but not a SotaSVM.")
    }
    if(is.null(clf))
    {
      stop("computeCoeffSVMLin: the clf is also needed.")
    }
    if(clf$params$kernel!="vanilladot")
    {
      # if not linear we don't compute the coeffs, does not mean anything
      return(mod)
    }
    p = length(mod$indices_)
    svm <- mod$obj
  }
  
  # compute the coefficients
  M <- diag(p) # create an empty diagonal matrix (sparse Matrix object)
  M <- rbind(M,rep(0,p)) # add intercept
  D <- predict(svm, M, type='decision') # compute probabilities
  intercept = -D[length(D)] # Compute the intercept
  
  D <- D + intercept # add the them of the intercept
  w <- D[1:length(D)-1] # compute the wi
  
  #print(c('coefs = ',D))
  #print(c('intercept=',intercpt))
  
  if(is.null(mod))
  {
    return(list(coeffs=w, intercept=intercept))
  }else
  {
    names(w) <- mod$names_
    mod$coeffs_ <- w
    return(mod)
  }
}


#' Evaluates wether an object is a model
#'
#' @description Evaluates wether an object is a model
#' @export
#' @param obj: an object to test
#' @return TRUE if the object is a model
isModel <- function(obj)
{
  # test if the model exists
  if(is.null(obj))
  {
    return(FALSE)
  }
  # test if the model is a list
  if(!is.list(obj))
  {
    return(FALSE)
  }
  
  # test weather it contains two main attributes
  res <- !is.null(obj$indices_) & !is.null(obj$names_)
  return(res)
}


#' Evaluates wether an object is a population of models
#'
#' @description Evaluates wether an object is a population of models
#' @param obj: an object to test
#' @return TRUE if the object is a population
#' @export
isPopulation <- function(obj)
{
  # test if the object exists
  if(is.null(obj))
  {
    return(FALSE)
  }
  
  # test if the list exists
  if(!is.list(obj))
  {
    return(FALSE)
  }
  
  # test if the list is not empty
  if(length(obj) == 0)
  {
    return(FALSE)
  }
  
  # test weather all the elements of the list are models
  if(!all(unlist(lapply(obj, isModel))))
  {
    return(FALSE)
  }
  
  return(TRUE)
}


#' Evaluates wether an object is a model collection objecct
#'
#' @description Evaluates wether an object is a model collection objecct
#' @param obj: an object to test
#' @return TRUE if the object is a model collection objecct
#' @export
isModelCollection <- function(obj)
{
  # test if the object exists
  if(is.null(obj))
  {
    return(FALSE)
  }
  
  # test if the list exists
  if(!is.list(obj))
  {
    return(FALSE)
  }
  
  # test if the list is not empty
  if(length(obj) == 0)
  {
    return(FALSE)
  }
  
  # test weather all the elements of the list are models
  if(!all(unlist(lapply(obj, isPopulation))))
  {
    return(FALSE)
  }
  return(TRUE)
}


#' Evaluates wether an object is a classifier
#'
#' @description Evaluates wether an object is a classifier
#' @export
#' @param obj: an object to test
#' @return TRUE if the object is a classifier
isClf <- function(obj)
{
  # test if the model exists
  if(is.null(obj))
  {
    return(FALSE)
  }
  # test if the model is a list
  if(!is.list(obj))
  {
    return(FALSE)
  }
  
  # test weather it contains two main attributes
  res <- !is.null(obj$params) & !is.null(obj$learner)
  return(res)
}

#' Evaluates wether an object is a model SOTA SVM
#'
#' @description Evaluates wether a learner is SOTA or not
#' @export
#' @param obj: a model to test
#' @return TRUE if the object is a SOTA learner
isLearnerSota <- function(obj)
{
  res <- isClf(obj) & ((length(grep("sota",obj$learner))==1) | (length(grep("logreg",obj$params$language))==1))
  return(res)
}

#' Evaluates wether an object is an experiment
#'
#' @description Evaluates wether an object is an experiment
#' @export
#' @param obj: an object to test
#' @return TRUE if the object is an experiment
isExperiment <- function(obj)
{
  # test if the model exists
  if(is.null(obj))
  {
    return(FALSE)
  }
  # test if the model is a list
  if(!is.list(obj))
  {
    return(FALSE)
  }
  
  # test weather it contains two main attributes
  res <- !is.null(obj$classifier) & isClf(obj$classifier)
  return(res)
}


#' Evaluates wether an object is a model SOTA 
#'
#' @description Evaluates wether an object is a model of type sota
#' @export
#' @param obj: a model to test
#' @return TRUE if the object is a model sota 
isModelSota <- function(obj)
{
  res <- isModelSotaSVM(obj) | isModelSotaRF(obj) | isModelSotaGLMNET(obj)
  return(res)
}

#' Evaluates wether an object is a model SOTA SVM
#'
#' @description Evaluates wether an object is a model of type sota
#' @export
#' @param obj: a model to test
#' @return TRUE if the object is a model sota SVM
isModelSotaSVM <- function(obj)
{
  res <- isModel(obj)
  if(!res)
  {
    return(res)
  }
  res <- (length(grep("sota.svm",obj$learner))==1)
  return(res)
}


#' Evaluates wether an object is a model SOTA RF
#'
#' @description Evaluates wether an object is a model of type sota
#' @export
#' @param mod: a model to test
#' @return TRUE if the object is a model sota RF
isModelSotaRF <- function(obj)
{
  res <- isModel(obj)
  if(!res)
  {
    return(res)
  }
  res <- (length(grep("sota.rf",obj$learner))==1) 
  return(res)
}

#' Evaluates wether an object is a model SOTA GLMNET
#'
#' @description Evaluates wether an object is a model of type sota
#' @param mod: a model to test
#' @return TRUE if the object is a model sota GLMNET
#' @export
isModelSotaGLMNET <- function(obj)
{
  res <- isModel(obj)
  if(!res)
  {
    return(res)
  }
  res <- (obj$learner == "terda") 
  if(is.null(obj$obj))
  {
    res <- FALSE
  }else
  {
    res <- res & !("character" %in% class(obj$obj))  
  }
  return(res)
}

#' Evaluates wether an object is a model BTR Terda
#'
#' @description Evaluates wether an object is a model of type BTR
#' @param mod: a model to test
#' @return TRUE if the object is a model BTR Terda
#' @export
isModelTerda <- function(obj)
{
  res <- isModel(obj)
  if(!res)
  {
    return(res)
  }
  res <- (obj$learner == "terda") 
  if(!obj$language %in% c("ter","terinter","bin","bininter","ratio"))
  {
    res <- FALSE
  }
  return(res)
}

#' Evaluates wether an object is a model BTR 
#'
#' @description Evaluates wether an object is a model of type BTR
#' @export
#' @param obj: a model to test
#' @return TRUE if the object is a model BTR 
isModelBTR <- function(obj)
{
  res <- isModel(obj) & !isModelSotaSVM(obj) & !isModelSotaRF(obj) & !isModelSotaGLMNET(obj)
  return(res)
}


#' Evaluates the accuracy of a model
#'
#' @description This function evaluates the accuracy of either (1) a model object that contains intercept and sign or (2) directly the attributes score, intercept, sign
#' @param mod: a model object to be used in the class prediction
#' @param X: dataset to classify
#' @param y: variable to predict
#' @param clf: an object containing the different parameters of the classifier
#' @param force.re.evaluation: evaluate again all the elements needed for accuracy (default:FALSE)
#' @param mode: training or test mode. If training, the funciton maximizes accuracy.
#' @return either (1) a model whose evaluation parameters are updated or (2) the accuracy
#' @export
evaluateAccuracy <- function(mod = NULL, X, y, clf, force.re.evaluation = FALSE, mode = "train")
{
  
  #If mod is not a valid model
  if(!isModel(obj = mod)) 
  {
    stop("evaluateAccuracy: please provide a valid model object BTR or SOTA")
  }else
  {
    # test if the confusion matrix exists
    if(!myAssertNotNullNorNa(mod$confusionMatrix_) | force.re.evaluation)
    {
      # NOTE: we consider that evaluateFit is the main function where we would have computed the score and intercept if force.re.evaluation.
      # here we need to update the confusionMatrix_
      
      # compute the score if it does not exist
      if(!myAssertNotNullNorNa(mod$score_))
      {
        scorelist       <- getModelScore(mod = mod, X = X, clf, force.re.evaluation = force.re.evaluation) # compute the score of the model
        mod$score_      <- scorelist$score_
        mod$pos_score_  <- scorelist$pos_score_
        mod$neg_score_  <- scorelist$neg_score_
        
        if(is.null(mod$score_))
        {
          return(NULL)
        }
      } 
      
      # compute the intercept and/or the sign if they do not exist
      if((!myAssertNotNullNorNa(mod$intercept_) | !myAssertNotNullNorNa(mod$sign_)) & !isModelSota(mod)) 
      {
        mod <- evaluateIntercept(mod = mod, X = X, y = y, clf = clf)
      }
      
      if(!isModelSota(mod))
      {
        if(!myAssertNotNullNorNa(mod$intercept_)) return(NULL)
        if(!myAssertNotNullNorNa(mod$sign_)) return(NULL)
      }
      
      # compute the confusion matrix
      mod$confusionMatrix_ <- computeConfusionMatrix(mod, X, y, clf)
      
      if(is.null(mod$confusionMatrix_))
      {
        return(NULL)
      }
    } # end re-evaluation or missing confusion matrix
    
    # EXPERIMENTAL TODO: add other accuraciy (weighted) 
    # compute a confusionMatrix that is weighted for a better AUC
    #props <- prop.table(table(y))
    #confusionMatrix_.prop <- apply(mod$confusionMatrix_, 1, fun <- function(x){x*rev(props)}) 
    # accuracy = (tp+tn)/(tp+fp+tn+fn)
    #a1 <- sum(diag(confusionMatrix_.prop)) / sum(confusionMatrix_.prop)
    # EXPERIMENTAL  
    
    if(mode == "train")
    {
      a1 <- sum(diag(mod$confusionMatrix_)) / sum(mod$confusionMatrix_)
      a2 <- sum(diag(apply(mod$confusionMatrix_, 2, rev))) / sum(mod$confusionMatrix_)
      
      if(a1 < a2  & !isModelSota(mod))
      {
        # inverse the sign
        if(mod$sign_ == ">") 
        {
          mod$sign_ <- "<"
        }else 
        {
          mod$sign_ <- ">"
        }
        
        # and recompute everything. This works with a recursive calling but due to R limits recursive calling in some deep cases
        # throws errors. This is why we are recoding another solution
        # mod <- evaluateAccuracy(mod = mod, X = X, y = y, force.re.evaluation = force.re.evaluation, mode = mode)
        
        # recompute the confusion matrix
        mod$confusionMatrix_ <- computeConfusionMatrix(mod, X, y, clf)
        if(is.null(mod$confusionMatrix_))
        {
          return(NULL)
        }
        
        # and accuracy
        mod$accuracy_ <- sum(diag(mod$confusionMatrix_)) / sum(mod$confusionMatrix_)
      }else
      {
        mod$accuracy_ <- a1  
      }
    }else
    {
      mod$accuracy_ <- sum(diag(mod$confusionMatrix_)) / sum(mod$confusionMatrix_)
    }
  } # end train/test
  
  return(mod)
}


#' Computes the ^y score of the model
#'
#' @description Returns the ^y score of the model
#' @importFrom kernlab predict
#' @param mod: a model object where the score will be computed
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param force.re.evaluation: we recompute the score (default:TRUE)
#' @return a vector containing the predicted ^y score for each observation
#' @export
getModelScore <- function(mod, X, clf, force.re.evaluation = TRUE)
{

  if(!isModel(obj = mod))
  {
    stop("getModelScore: please provide a valid model object BTR or SOTA!")    
  }
  
  if(any(is.na(mod$indices_)))
  {
    printModel(mod); print(mod$indices_)
    if(clf$params$warnings) warning("getModelScore: indices contain NAs")
    return(NULL)
  }
  
  if(max(mod$indices_) > nrow(X))
  {
    printModel(mod); print(mod$indices_)
    if(clf$params$warnings) warning(paste("getModelScore: indices should smaller then", nrow(X)))
    return(NULL)
  }

  # If we already have the score no need to recompute it
  if(myAssertNotNullNorNa(mod$score_) & !force.re.evaluation)
  {
    res <- list(score_ = as.numeric(mod$score_),
                pos_score_ = as.numeric(mod$pos_score_),
                neg_score_ = as.numeric(mod$neg_score_)
    )
    return(res)
  } # end score exists
  
  # if model is sota we need to recompute score
  if(isModelSotaSVM(obj = mod))
  {
    if(is.null(mod$obj))
    {
      if(clf$params$warnings) warning("getModelScore: please provide a valid SVM model in mod$obj!")
      return(NULL)
    }
    # else we compute it again
    res <- list(score_ = as.numeric(predict(mod$obj, t(X[mod$indices_,]), type="decision")),
                pos_score_ = as.numeric(mod$pos_score_),
                neg_score_ = as.numeric(mod$neg_score_)
    )
    return(res)
  } # end SVM model
  
  if(isModelSotaRF(obj = mod))
  {
    if(is.null(mod$obj))
    {
      if(clf$params$warnings) warning("getModelScore: please provide a valid random forest model in mod$obj!")
      return(NULL)
    }
    
    prob <- NA # flag it to check wether it will work
    # else we compute it again
    try(prob  <-  predict(mod$obj, t(X[mod$indices_,]), type="prob"), silent = TRUE)
    if(length(prob) == 1)
    {
      if(is.na(prob))
      {
        if(clf$params$warnings) warning("getModelScore: the RF probability score could not be computed!")
        return(NULL)
      }
    }
    
    if(mean(prob[,1]) > mean(prob[,2]))
    {
      yhat <- prob[,1]
    }else
    {
      yhat <- prob[,2]
    }
    
    # else we compute it again
    res <- list(score_ = as.numeric(yhat),
                pos_score_ = as.numeric(mod$pos_score_),
                neg_score_ = as.numeric(mod$neg_score_)
    )
    return(res)
  } # end RF model
  
  
  # if model is sota we need to recompute score
  if(isModelSotaGLMNET(obj = mod))
  {
    if(is.null(mod$obj))
    {
      if(clf$params$warnings) warning("getModelScore: please provide a valid random forest model in mod$obj!")
      return(NULL)
    }
    
    prob <- NA # flag it to check wether it will work
    # else we compute it again
    #try(prob  <-  predict(mod$obj, data.matrix(t(X[mod$indices_,]), s = mod$lambda), type="prob"), silent = TRUE)
    try(prob  <-  predict(object = mod$obj, s = mod$lambda, newx = data.matrix(t(X)))[,1], silent = TRUE)
    if(length(prob) == 1)
    {
      if(is.na(prob))
      {
        if(clf$params$warnings) warning("getModelScore: the GLMNET probability score could not be computed!")
        return(NULL)
      }
    }
    
    # else we compute it again
    res <- list(score_ = as.numeric(prob),
                pos_score_ = as.numeric(mod$pos_score_),
                neg_score_ = as.numeric(mod$neg_score_)
    )
    return(res)
  } # end GLMNET model
  
  
  # for NAs set coefficients randomly
  if(any(is.na(mod$coeffs_)))
  {
    ind                 <- which(is.na(mod$coeffs_)) # find which coefs are NA
    set.seed(clf$params$seed) # fix seed
    replacement.coeffs  <- sample(c(-1,1),length(ind),replace = TRUE) # sample randomly
    mod$coeffs_[ind]    <- replacement.coeffs # replace
    # and give a warning
    if(clf$params$warnings) warning("getModelScore: there are NAs in the model's coefficents. This have been set randomly to -1 or 1.") 
  }
  
  # if this is classification problem (TODO test other AIC etc ...)
  if(length(mod$indices_) == 1) 
  { 
    if(mod$language == "ratio") # if custom language (ratio, ...)
    {
      score <- rep(0, length(ncol(X)))
      pos_score <- rep(0, ncol(X))
      neg_score <- rep(0, ncol(X))
    }else
    {
      # if k_sparsity 1 no need to use a complex formula
      data <- t(as.matrix(X[mod$indices_,]))
      score <- mod$coeffs_ * data # compute score ^y  

      if(sign(mod$coeffs_) > 0)
      {
        pos_score <- score
        neg_score <- rep(0, ncol(X))
      }else
      {
        pos_score <- rep(0, ncol(X))
        neg_score <- score
      }
    }
  }else 
  { 
    
    # if size > 1 we spit the positive and negative parts
    pos <- which(mod$coeffs_ > 0)
    neg <- which(mod$coeffs_ < 0)
    pos <- list(indices_ = mod$indices_[pos], coeffs_ = mod$coeffs_[pos])
    neg <- list(indices_ = mod$indices_[neg], coeffs_ = mod$coeffs_[neg])
    
    if(length(pos$indices_) == 0) 
    {
      # send an zero filled vector (another idea to send out a NULL object)
      pos_score <- rep(0, ncol(X))
      names(pos_score) <- colnames(X)
      pos_score <- t(as.matrix(pos_score))
    } else if(length(pos$indices_) == 1) # if the positive part is of length one
    {
      pos_data <- as.matrix(X[pos$indices_,])
      pos_score <- pos_data
    } else
    { 
      pos_data <- X[pos$indices_,]
      pos_score <- pos$coeffs_ %*% as.matrix(pos_data) # compute score ^y
    }
    
    if(length(neg$indices_) == 0)
    {
      # send an zero filled vector (another idea to send out a NULL object)
      neg_score <- rep(0, ncol(X)); names(neg_score) <- colnames(X)
      neg_score <- t(as.matrix(neg_score))
    } else if(length(neg$indices_) == 1) 
    {
      # if the positive part is of length one
      neg_data <- as.matrix(X[neg$indices_,])
      neg_score <- neg_data
    } else # If greater than one
    {
      neg_data <- X[neg$indices_,]
      neg_score <- (-neg$coeffs_) %*% as.matrix(neg_data) # compute score ^y
    }
    
    # compute the score using contextual information from the model not from the clf
    if(mod$language=="ratio") # if custom language (ratio, ...)
    {
      score <- clf$params$scoreFormula(class_1_score = pos_score, class_2_score = neg_score, epsilon = clf$params$epsilon)
    } else # else default language ter + teta
    {
      
      if(nrow(pos_score) < ncol(pos_score)) 
      {
        pos_score <- t(pos_score)
      }
      if(nrow(neg_score) < ncol(neg_score)) 
      {
        neg_score <- t(neg_score)
      }
      score <- pos_score - neg_score
    }
  }
  
  res <- list(score_ = as.numeric(score),
              pos_score_ = as.numeric(pos_score),
              neg_score_ = as.numeric(neg_score)
  )
  
  return(res)
}


#' Computes the AUC of a model
#'
#' @description Computes the AUC of a model
#' @param score: the ^y score of the model
#' @param y: the response vector
#' @param sign: in which direction to make the comparison? "auto" (default): automatically define in which group 
#' the median is higher and take the direction accordingly. ">": if the predictor values for the control group 
#' are higher than the values of the case group (controls > t >= cases). "<": if the predictor values for the 
#' control group are lower or equal than the values of the case group (controls < t <= cases).
#' @return an auc value
#' @importFrom pROC roc
evaluateAUC <- function(score, y, sign = '>')
{
  # NOTE: in the score we may have infite values that come for instance from the ratio language
  # This means that whatever the intercept these examples are positive ones. As such we can omit
  # them when computing the intercept.
  ind.infinite  <- is.infinite(score)
  ind.nan       <- is.nan(score)
  ind.filter    <- ind.infinite | ind.nan
  
  # if all the observations are not numbers return NA
  if(sum(ind.filter) == length(y))
  {
    return(NA)
  }
  
  score <- score[!ind.filter]
  y     <- y[!ind.filter]
  
  # if y does not contain exactly 2 levels, and if these don't have at least 
  # 1 count then we don't compute the AUC
  if(length(table(y)) != 2 | any(table(y)==0))
  {
    auc <- NA
    
  }else # otherwise we compute it
  {
    rocobj <- suppressMessages(roc(response = y, predictor = score, direction = sign))
    auc <- as.numeric(rocobj$auc)
  }
  
  if(sign == "auto")
  {
    resg <- evaluateAUC(score = score, y = y, sign = ">")
    resl <- evaluateAUC(score = score, y = y, sign = "<")
    auc <- max(resg, resl)
  }
  
  return(auc)
}


#' Evaluates the fitting score of a model object
#' 
#' @description Evaluates the fitting score of a model object.
#' @param mod : a model object
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the response vector
#' @param clf: the classifier parameter object
#' @param force.re.evaluation: re-evaluate all the scores even if they exist (default:FALSE)
#' @param mode: A choice from c("train", "test") indicates wether we wish to learn the threthold 
#' of the model (default:"train") or not "test" for the c("terinter","bininter","ratio") languages
#' @return a model object with the fitting score
evaluateFit <- function(mod, X, y, clf, force.re.evaluation = FALSE, mode = "train")
{
  # if the models is not a valid object
  if(!isModel(mod))
  {
    if(!is.character(mod) | !is.numeric(mod))
    {
      stop("evaluateFit: please provide a valid model object or a feature index vector.")
    }
    
    if(is.character(mod))
    { 
      # if model is of the form of variable names
      mod <- names2index(X = X, var.names = mod)
    }
    mod <- individual(X, y, clf = clf, ind = mod)
  }
  
  # compute the score in the case it is asked to recompute
  if(force.re.evaluation)
  {
    scorelist <- getModelScore(mod = mod, X = X, clf = clf, force.re.evaluation = force.re.evaluation)
    mod$score_ <- scorelist$score_
    mod$pos_score_ <- scorelist$pos_score_
    mod$neg_score_ <- scorelist$neg_score_
  }else
  {
    # compute the score if it does not exist
    if(!myAssertNotNullNorNa(mod$score_))
    {
      scorelist <- getModelScore(mod = mod, X = X, clf = clf, force.re.evaluation = force.re.evaluation)
      # if(!any(is.na(scorelist)))
      # {
      #   mod$score_ <- scorelist$score_
      #   mod$pos_score_ <- scorelist$pos_score_
      #   mod$neg_score_ <- scorelist$neg_score_  
      # }
      if(!any(is.na(scorelist)) | isModelSota(mod))
      {
        mod$score_ <- scorelist$score_
        mod$pos_score_ <- scorelist$pos_score_
        mod$neg_score_ <- scorelist$neg_score_  
      }
    }else
    {
      # in the case the score has been computed before but for an other X, we recompute
      if(length(mod$score_) != ncol(X))
      {
        scorelist <- getModelScore(mod = mod, X = X, clf = clf, force.re.evaluation = force.re.evaluation)
        # if(!any(is.na(scorelist)))
        # {
        #   mod$score_ <- scorelist$score_
        #   mod$pos_score_ <- scorelist$pos_score_
        #   mod$neg_score_ <- scorelist$neg_score_  
        # }
        if(!any(is.na(scorelist)) | isModelSota(mod))
        {
          mod$score_ <- scorelist$score_
          mod$pos_score_ <- scorelist$pos_score_
          mod$neg_score_ <- scorelist$neg_score_  
        }
      }
    }
  }
  
  # if after all the above steps we still don't have a score than we kill the model.
  if(!myAssertNotNullNorNa(mod$score_))
  {
    return(NULL)
  }
  
  if(is.null(mod$eval.sparsity)) # if sparsity is not set
  {
    mod$eval.sparsity <- length(mod$indices_)
  }
  
  switch(clf$params$objective, 
         # THE AUC objective
         auc = {
           # compute the intercept and sign
           if(mode == 'train')
           {
             if((mod$language != "ter" & mod$language != "bin") & !isModelSota(mod))
             {
               mod                <- evaluateIntercept(X = X, y = y, clf = clf, mod = mod)  
             }else
             {
               # force intercept to NA for sota
               if(isModelSota(mod))
               {
                 mod$intercept_   <- NA
                 mod$sign_        <- NA
               }else
               {
                 mod$intercept_   <- 0  
               }
             }
             
             # sanity check
             if(!clf$params$evalToFit %in% names(mod))
             {
               stop("evaluateFit: the evalToFit parameter seems not to be a valid one. Please make sure it is among the available ones")
             }
             
             # if the following attributes are selected, then we need to fix it since they are derivates of a score
             if(clf$params$evalToFit != "fit_" | clf$params$evalToFit != "unpenalized_fit_")
             {
               clf$params$evalToFit <- "accuracy_"
             }
             
             # in case it is auc
             if(clf$params$evalToFit == "auc_") # in this case the auc will be computed in evaluate other metrics
             {
               # compute the auc
               aucg                 <- evaluateAUC(score = mod$score, y = y, sign = ">")
               aucl                 <- evaluateAUC(score = mod$score, y = y, sign = "<")
               mod$auc_             <- max(aucg, aucl)
               mod$unpenalized_fit_ <- mod$auc_             
             }
             
             # in case it is accuracy
             if(clf$params$evalToFit == "accuracy_") # in this case the auc will be computed in evaluate other metrics
             {
               mod <- evaluateAccuracy(mod, X, y, clf, force.re.evaluation = force.re.evaluation, mode = mode)
               mod$unpenalized_fit_ <- mod$accuracy_             
             }
             
             # otherwise compute the rest 
             if(clf$params$evalToFit != "auc_" & clf$params$evalToFit != "accuracy_")
             {
               mod <- evaluateAdditionnalMetrics(mod = mod, X = X, y = y, clf = clf, mode = mode)
               mod$unpenalized_fit_ <- mod[[clf$params$evalToFit]]
               # compte accuracy also
               mod <- evaluateAccuracy(mod, X, y, clf, force.re.evaluation = force.re.evaluation, mode = mode)
               # and auc, since these are helpful information
               aucg                 <- evaluateAUC(score = mod$score, y = y, sign = ">")
               aucl                 <- evaluateAUC(score = mod$score, y = y, sign = "<")
               mod$auc_             <- max(aucg, aucl)
             }
           } # if test mode, we don't recompute the intercept
           else
           {
             # sanity check
             if(!clf$params$evalToFit %in% names(mod))
             {
               stop("evaluateFit: the evalToFit parameter seems not to be a valid one. Please make sure it is among the available ones")
             }
             
             # if the following attributes are selected, then we need to fix it since they are derivates of a score
             if(clf$params$evalToFit != "fit_" | clf$params$evalToFit != "unpenalized_fit_")
             {
               clf$params$evalToFit <- "accuracy_"
             }
             
             # in case it is auc
             if(clf$params$evalToFit == "auc_") # in this case the auc will be computed in evaluate other metrics
             {
               # compute the auc
               aucg                 <- evaluateAUC(score = mod$score, y = y, sign = ">")
               aucl                 <- evaluateAUC(score = mod$score, y = y, sign = "<")
               mod$auc_             <- max(aucg, aucl)
               mod$unpenalized_fit_ <- mod$auc_      
               
               # compute the accuracy as it won't be computed elsewhere
               mod <- evaluateAccuracy(mod, X, y, clf, force.re.evaluation = force.re.evaluation, mode = mode)
               
             }
             
             # in case it is accuracy
             if(clf$params$evalToFit == "accuracy_") # in this case the auc will be computed in evaluate other metrics
             {
               mod <- evaluateAccuracy(mod, X, y, clf, force.re.evaluation = force.re.evaluation, mode = mode)
               mod$unpenalized_fit_ <- mod$accuracy_
               
               # compute also the AUC as otherwise it won't be computed anywhere
               aucg                 <- evaluateAUC(score = mod$score, y = y, sign = ">")
               aucl                 <- evaluateAUC(score = mod$score, y = y, sign = "<")
               mod$auc_             <- max(aucg, aucl)
             }
             
             # otherwise compute the rest when evalToFit is not auc_ nor accuracy_
             if(clf$params$evalToFit != "auc_" & clf$params$evalToFit != "accuracy_")
             {
               mod <- evaluateAdditionnalMetrics(mod = mod, X = X, y = y, clf = clf, mode = mode)
               mod$unpenalized_fit_ <- mod[[clf$params$evalToFit]]
               # compute accuracy also
               mod <- evaluateAccuracy(mod, X, y, clf, force.re.evaluation = force.re.evaluation, mode = mode)
               # and auc, since these are helpful information
               aucg                 <- evaluateAUC(score = mod$score, y = y, sign = ">")
               aucl                 <- evaluateAUC(score = mod$score, y = y, sign = "<")
               mod$auc_             <- max(aucg, aucl)
             }
           } # end mode = test
           
         },
         # THE REGRESSION correlation method
         cor = { 
           #added may 8th 2016 fix the mix of population negative & positive
           
           tryCatch({
             # as of 2018/07/09 the regression process changes, We will search 
             # for a preason correlation and will maximise the r2. We also 
             # implemented standard error of the mean (which needs to be minimized)
             ina      <- is.na(mod$score_) | is.na(y) | is.infinite(mod$score) | is.infinite(y)
             mod$cor_ <- abs(cor(mod$score[!ina], y[!ina], method = "pearson"))
             
             # # for optimization reasons the cor will be a pearson implementation, 
             # y is already ranked for objective being "cor" performed in the fit()
             # # We just need to tank the score to obtain the same result as a spearman correlation.
             #abs(cor.test(mod$score, y, objective = "spearman")$estimate)
             #lmfit    <- lm(mod$score_~y)
             #mod$ser_ <- summary(lmfit)$coefficients[,2][2]
             #mod$rsq_ <- summary(lmfit)$r.squared
             
             # use the r2 instead
             score.scaled <- as.vector(scale(mod$score_[!ina], center = TRUE, scale = TRUE))
             y.scaled <- as.vector(scale(y[!ina], center = TRUE, scale = TRUE))
             mod$rsq_ <- abs(cor(mod$score[!ina], y[!ina], method = "pearson"))^2
             mod$ser_ <- sqrt(sum((score.scaled - y.scaled)^2, na.rm = TRUE)/(length(score.scaled) - 2))
             
           })
           # , warning = function(war) 
           # {
           #   # warning handler picks up where error was generated
           #   if(clf$params$debug) print(paste("Warning",err, "Setting result to NA"))
           #   mod <- NULL
           # }, error = function(err) 
           # {
           #   # error handler picks up where error was generated
           #   if(clf$params$debug) print(paste("ERROR",err, "Setting result to NA"))
           #   mod <- NULL
           # }, finally = {
           #   # NOTE:  Finally is evaluated in the context of the inital tryCatch block 
           #   # and 'err' will not exist if a warning or error occurred.
           # })
           
           if(is.null(mod)) return(mod)
           if(myAssertNotNullNorNa(mod$cor_)) 
           {
             mod$cor_           <- as.numeric(mod$cor_)
           }
           if(myAssertNotNullNorNa(mod$rsq_)) 
           {
             mod$rsq_           <- as.numeric(mod$rsq_)
           }
           if(myAssertNotNullNorNa(mod$ser_)) 
           {
             mod$ser_           <- as.numeric(mod$ser_)
           }
           
           # get the value to maximize in the general optimization variable
           #mod$unpenalized_fit_ <- mod$rsq_
           mod$unpenalized_fit_ <- mod$rsq_
         },
         aic={ # THE AIC objective
           # TODO test it out
           mod$aic_             <- estimateCoefficientsIndividual(X=X, y=y, ind = mod$indices_)$aic
           mod$unpenalized_fit_ <- mod$aic_
         },
         { # else
           if(clf$params$warnings) warning('This objective method does not exist !')
         }
  )
  
  # apply the penalty based on model size
  mod$fit_ <- max(mod$unpenalized_fit_ - clf$params$k_penalty * mod$eval.sparsity, 0)
  #mod$fit_ <- mod$fit_ - clf$params$k_penalty * sqrt(mod$eval.sparsity) # Square root when to make it softer
  
  return(mod)
}


#' Computes the best intercept for the model while minimizing error
#' 
#' @description Computes the best intercept for the model
#' @param score: the ^y score of the model
#' @param y: the response vector
#' @param verbose: print running information when set to TRUE
#' @param sign: weather the score should be greater or smaller than the intercept (default:"auto")
#' @param return.all: if TRUE, the function will return the intercept as well as the table used to compute it.
#' @param plot: if TRUE, the score will be visialized (default:FALSE)
#' @return the intercept, the sign and the accuracy
#' @export
computeIntercept <- function(score, y, verbose = FALSE, sign = "auto", plot = FALSE) {
  # make vectors of 0/1 to identify positive labels and negative labels
  
  if(!is.factor(y))
  {
    y <- as.factor(y)
  }
  
  lev <- levels(y)
  if(length(lev) != 2)
  {
    stop("computeIntercept: please make sure the y class has only two levels.")
  }
  
  # if(!all(levels(y) %in% c("-1",1)))
  # {
  #   stop("computeIntercept: please make sure the y class levels are -1 and 1.")
  # }
  
  if(sign == "auto") {
    sup_res = computeIntercept(score, y, verbose, sign = ">")
    inf_res = computeIntercept(score, y, verbose, sign = "<")
    
    if(sup_res$accuracy > inf_res$accuracy)
    {
      return(sup_res)
    }else
    {
      return(inf_res)
    }
  }
  
  # PosClass <- (y== 1)+0
  # NegClass <- (y==-1)+0
  
  # replaced on 21/02/18
  PosClass <- (y == lev[1])+0
  NegClass <- (y == lev[2])+0
  
  CountNeg <- sum(NegClass)
  CountPos <- sum(PosClass)
  
  stopifnot(sign=='>' | sign =='<')
  
  # put all that in a dataframe 'm', sorted by increasing score values
  # in case identical scores appear, merge them with "aggregate"
  
  #m.raw2 <- m.raw
  if(length(score) != length(PosClass) | length(score) != length(NegClass))
  {
    if(clf$params$warnings) warning("ComputeIntercept: score should be the same length as the other variables, returning NULL")
    return(NULL)
  }
  
  m.raw <- data.frame(score, PosClass, NegClass)   # combine them in a frame
  m.raw <- m.raw[order(m.raw$score),]  # sort the frame by increasing scores, and rearrange
  m <- m.raw

  # if(length(unique(m.raw$score)) != length(m.raw$score)) # if multiple values to be aggregated
  # {
  #   m <- aggregate(x = m.raw, by = list(m.raw$score), FUN = "sum")
  #   m <- m[,-1]
  # }else #otherwise no need to aggregate this is much faster
  # {
  #   m <- m.raw
  # }
  
  #plot(m[,1], m[,2])
  
  # the following line adds a fake score s, below the minimum one, such that "sum(wi.xi) > s" includes all observations
  # if(nrow(m)>1)
  # {
  #   #m <- rbind(c(2*m[1,1] - m[2,1],0,0),m) # optional
  #   # a faster implementatio than rbind
  #   nextrow = nrow(m)+1
  #   m[nextrow:(nextrow+1),] = c(2*m[1,1] - m[2,1],0,0)
  #   # we need to assure unique row names
  #   row.names(m) = 1:nrow(m)
  # }
  
  FN <- cumsum(m["PosClass"])
  TN <- cumsum(m["NegClass"])
  # add the true negative counts and false negative counts to dataframe 'm' for the sign '>'
  if (sign == '>') 
  {
    ERR <- FN + CountNeg - TN
  } else 
  {
    ERR <- TN + CountPos - FN
  }
  colnames(FN) <- "FN"
  colnames(TN) <- "TN"
  colnames(ERR) <- "ERR"
  m <- cbind.data.frame(m, FN, TN, ERR)
  #print(m)
  
  min.err.ind <- which.min(m[,"ERR"])
  intercept <- m[min.err.ind,]$score
  #abline(h=intercept, lty=2, col = "blue")
  #abline(v=min.err.ind, lty = 2, col = "blue")
  res <- list(intercept = intercept, 
              sign = sign, 
              accuracy = (1-(min(m[,"ERR"])/length(y)))
  )
  
  if(plot)
  {
    par(mfrow=c(2,1))
    # plot score
    plot(m$score, cex = 1, col = 'white', ylab = "model score", xlab = "observations ordered by score")
    points(which(m$PosClass==1),m$score[m$PosClass==1], col="firebrick", pch="*", cex = 2)
    points(which(m$NegClass==1),m$score[m$NegClass==1], col="darkblue", pch="*", cex = 2)
    abline(h = intercept, col = "red", lty = 2)
    # plot error
    plot(m$ERR, cex = 1, col = 'black', pch = 19, ylab = "cumulative error", xlab = "observations ordered by score")
    abline(v = min.err.ind, col = "red", lty = 2)
    par(mfrow=c(1,1))
  }

  return(res)
}


#' Evaluates the fitting score of a model object
#'
#' @description Evaluates the fitting score of a model object.
#' @param mod : a model object
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the response vector
#' @param clf: the classifier parameter object
#' @return a model object with the fitting score
evaluateIntercept <- function(mod, X, y, clf)
{
  if(!isModel(mod))
  { # if model is not an object
    if(is.character(mod))
    { # if model is of the form of variable names
      mod <- index2names(X, mod)
    }
    mod <- individual(X, y, clf = clf, ind = mod)
  }
  
  if(isModelSota(mod))
  {
    if(clf$params$warnings) warning("evaluateIntercept: no intercept for sota models. Returning unchanged.")
    return(mod)
  }
  
  if((clf$params$intercept=="NULL"))
  {
    
    # compute the fitting scores depending on the method
    if(!myAssertNotNullNorNa(mod$score_))
    {
      scorelist <- getModelScore(mod = mod, X = X, clf = clf)
      mod$score_ <- scorelist$score_
      mod$pos_score_ <- scorelist$pos_score_
      mod$neg_score_ <- scorelist$neg_score_
      
      if(is.null(mod$score_))
      {
        return(NULL)
      }
    }
    
    # NOTE: in the score we may have infite values that come for instance from the ratio language
    # This means that whatever the intercept these examples are positive ones. As such we can omit
    # them when computing the intercept.
    ind.infinite  <- is.infinite(mod$score_)
    ind.nan       <- is.nan(mod$score_)
    
    score         <- mod$score_
    # the ind.infinite scores should be replaced by a big value let say max(X)+1
    try(score[ind.infinite] <- clf$data$X.max + 1, silent = TRUE)
    # the ind.nan should be replaced by a small value for instance min(X)-1
    try(score[ind.nan] <- clf$data$X.min - 1, silent = TRUE)
    
    
    switch(clf$params$objective, 
           auc={ # THE AUC objective
             interc             <- computeIntercept(score = score, y = y, clf)
             mod$intercept_     <- interc$intercept
             mod$sign_          <- interc$sign
           },
           cor={
             mod$intercept_     <- "NULL"
             mod$sign_          <- "NULL"
           },
           {
             if(clf$params$warnings) warning("evaluateIntercept: This objective does not exist !")
           }
    )
  } else
  {
    mod$intercept_ <- clf$params$intercept
    mod$sign_ <- ">"
  }
  
  return(mod)
}


#' Evaluates the fitting coefficents of a model object
#'
#' @description Evaluates the fitting coefficients of a model object.
#' @param mod: a model object
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the response vector
#' @param clf: the classifier parameter object
#' @param eval.all: should the function evaluate all the scores (default:FALSE)
#' @param force.re.evaluation: re-evaluate all the scores even if they exist (default:FALSE)
#' @return a model object with the fitting scores evaluated
#' @export 
evaluateModelRegression <- function(mod, X, y, clf, eval.all = FALSE, force.re.evaluation = FALSE)
{
  if(!isModel(mod))
  {
    stop("evaluateModelRegression: the model to be evaluated does not exist.")
  }
  
  if(isModelSota(mod))
  {
    if(clf$params$warnings) warning("evaluateModelRegression: no intercept for sota models. Returning unchanged.")
    return(mod)
  }
  
  mod.res <- mod
  
  # If sparsity is not the same
  if(mod.res$eval.sparsity > length(unique(mod.res$indices_)))
  {
    if(clf$params$warnings) warning("evaluateModelRegression: An individual has at least one indice replicated")
    values2keep             <- which(mod.res$indices_ == unique(mod.res$indices_))
    mod.res$indices_        <- mod.res$indices_[values2keep]
    mod.res$names_          <- mod.res$names_ [values2keep]
    mod.res$coeffs_         <- mod.res$coeffs[values2keep]
    mod.res$eval.sparsity   <- length(unique(mod.res$indices_))
  }
  
  # Compute score
  if(is.null(mod.res$score_) | force.re.evaluation)
  {
    scorelist <- getModelScore(mod = mod.res, X = X, clf = clf)
    mod.res$score_ <- scorelist$score_
    mod.res$pos_score_ <- scorelist$pos_score_
    mod.res$neg_score_ <- scorelist$neg_score_
    
    if(is.null(mod.res$score_))
    {
      return(NULL)
    }
  }
  mod.res                   <- evaluateFit(mod = mod.res, X=X, y=y, clf=clf, force.re.evaluation = force.re.evaluation)
  return(mod.res)
}


#' Evaluates the fitting score of a model object
#'
#' @description Evaluates the fitting score of a model object.
#' @param mod: a model object
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the response vector
#' @param clf: the classifier parameter object
#' @param eval.all: should the function evaluate all the scores (default:FALSE)
#' @param force.re.evaluation: re-evaluate all the scores even if they exist (default:FALSE)
#' @param estim.feat.importance: evaluate the importance in the model object (default:FALSE)
#' @param mode: A choice from c("train", "test") indicates wether we wish to learn the threthold 
#' of the model (default:"train") or not "test" for the c("terinter","bininter","ratio") languages
#' @return a model object with the fitting scores evaluated
#' @export 
evaluateModel <- function(mod, X, y, clf, eval.all = FALSE, force.re.evaluation = FALSE, estim.feat.importance = FALSE, mode = 'train')
{
  
  if(mode != "train" & mode != "test")
  {
    stop("evaluateModel: mode should be one of c('train','test')")
  }
  
  if(!isModel(mod))
  {
    # if not a model object but a valid index, create a model
    if(!is.list(mod))
    {
      mod <- individual(X = X, y = y, clf = clf, ind = mod) # transform into a model
    }else
    {
      if(clf$params$warnings) warning("evaluateModel: the model to be evaluated does not exist, returning NULL.")
      return(NULL)  
    }
  }
  
  # DON'T EVALUATE RATIO, TER and TERINTER MODELS WITHOUT NEGATIVE AND POSITIVE TERMS
  if(!isModelSota(mod))
  {
    # at this stage the model should be a valid one. If not return NULL
    if(!isModel(mod))
    {
      if(clf$params$warnings) warning("evaluateModel: the model to be evaluated does not exist and at this stage it should be one, returning NULL.")
      return(NULL)  
    }
    
    if(mod$language == "ratio" | mod$language == "ter" | mod$language == "terinter")
    {
      if(length(table(sign(mod$coeffs_))) != 2)
      {
        return(NULL)
      }
    }
  }
  
  # make a copy of the model object
  mod.res <- mod
  
  if(mode == "train")
  {
    # reset the attributes that need to be recomputed
    # general attributes
    mod.res$fit_                <- NA
    mod.res$unpenalized_fit_    <- NA
    # classification
    mod.res$auc_                <- NA
    mod.res$accuracy_           <- NA
    mod.res$precision_          <- NA
    mod.res$recall_             <- NA
    mod.res$f1_                 <- NA
    mod.res$intercept_          <- NA # the intercept
    mod.res$sign_               <- NA # the sign of the model
    # regression
    mod.res$cor_                <- NA
    mod.res$rsq_                <- NA # r2
    mod.res$ser_                <- NA # standard error of the mean
    mod.res$aic_                <- NA
    mod.res$score_              <- NA # model score
    mod.res$pos_score_          <- NA # the positive model score
    mod.res$neg_score_          <- NA # the negative model score
  }
  
  
  # If this is a regression problem no need to find intercepts etc...
  if(clf$params$objective == "cor")
  {
    if(!isModelSota(mod.res))
    {
      mod.res <- evaluateModelRegression(mod = mod.res, X = X, y = y, clf = clf, eval.all = eval.all, force.re.evaluation = force.re.evaluation)
      return(mod.res)
    }else
    {
      if(clf$params$warnings) warning("evaluateModel: evaluating a sota model in correlation objective")
    }
  }
  
  if(isModelSota(mod.res))
  {
    # feature importance estimation will be switched off for the sotas, since the internal model structure is very different
    if(estim.feat.importance) 
    {
      estim.feat.importance = FALSE
    }
  }
  
  # if(mod.res$eval.sparsity != length(mod.res$indices_))
  # {
  #   stop("evaluateModel: The sparsity of the model does not march the number of indices. This is not normal... killing the execution")
  # }
  
  # If sparsity is not the same
  if(mod.res$eval.sparsity != length(unique(mod.res$indices_)))
  {
    if(clf$params$warnings) warning("evaluateModel: An individual has at least one indice replicated")
    values2keep           <- which(mod.res$indices_ == unique(mod.res$indices_))
    mod.res$indices_      <- mod.res$indices_[values2keep]
    mod.res$names_        <- mod.res$names_[values2keep]
    mod.res$coeffs_       <- mod.res$coeffs[values2keep]
    mod.res$eval.sparsity <- length(unique(mod.res$indices_))
  }
  
  # first evaluate the fit
  mod.res <- evaluateFit(mod = mod.res, X=X, y=y, clf=clf, force.re.evaluation = force.re.evaluation, mode = mode)
  
  # compute all other evaluation metrics
  if(eval.all)
  {
    # At this stage this should not happen but for debug stop it
    if((!myAssertNotNullNorNa(mod.res$intercept_) | !myAssertNotNullNorNa(mod.res$sign_)) & !isModelSota(mod.res)) 
    {
      if(clf$params$warnings) warning("evaluateModel: model without intercept at this stage is not normal.")
      return(NULL)
    }
    
    mod.res         <- evaluateAdditionnalMetrics(mod = mod.res, X = X, y = y, clf = clf, mode = mode)
    if(!isModel(mod.res))
    {
      # if(clf$params$warnings) warning("evaluateModel: returning an empty model.")
      return(NULL)
    }
  }
  
  # if BTR
  if(!isModelSota(mod.res))
  {
    if(estim.feat.importance)
    {
      mod.res         <- estimateFeatureImportance(mod = mod.res, X = X, y = y, clf = clf, plot.importance = FALSE)
    }
  }
  
  # At this stage this should not happen but for debug stop it
  if(!myAssertNotNullNorNa(mod.res$unpenalized_fit_)) 
  {
    if(clf$params$warnings) warning("evaluateModel: model does not have a valid unpenalized_fit_ attribute.")
    return(NULL)
  }
  
  return(mod.res)
}


#' Estimates the importance of each feature in the model object
#'
#' @description Estimates the importance of each feature in the model object
#' @param mod: a model object
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the response vector
#' @param clf: the classifier parameter object
#' @param attribute: which attribute should be used to compute the importance (default:unpenalized_fit_)
#' @param plot.importance: should the function plot the improtance of the features (default:FALSE)
#' @return a model object with the importance of each feature computed. Negative importance of a feature means that the feature is not beneficial.
#' @export 
estimateFeatureImportance <- function(mod, X, y, clf, attribute = "unpenalized_fit_", plot.importance = FALSE)
{
  if(!isModel(obj = mod))
  {
    stop("estimateFeatureImportance: please provide a valid predomics model")
  }
  
  if(isModelSota(mod))
  {
    if(clf$params$warnings) warning("estimateFeatureImportance: estimating feature importance is active only for BTR models")
    return(mod)
  }
  
  # recompute the model if the attribute is NA
  if(is.na(mod[[attribute]]))
  {
    mod                       <- evaluateModel(mod = mod.tmp, X = X, y = y, clf = clf, 
                                               eval.all = TRUE, force.re.evaluation = TRUE, mode = "train")
    
  }
  
  importance <- rep(NA, length(mod$indices_))
  names(importance) <- mod$names_
  
  if(length(mod$indices_) == 1)
  {
    importance[1]             <- mod[[attribute]]
  }else
  {
    for(i in 1:length(mod$indices_))
    {
      mod.tmp <- mod
      
      # omit one feature
      mod.tmp$indices_        <- mod.tmp$indices_[-i]
      mod.tmp$names_          <- mod.tmp$names_[-i]
      mod.tmp$coeffs_         <- mod.tmp$coeffs_[-i]
      mod.tmp$eval.sparsity   <- mod.tmp$eval.sparsity -1
      
      # recompute the model with one feature less
      mod.tmp                 <- evaluateModel(mod = mod.tmp, X = X, y = y, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE, mode = "train")
      
      # compute the importance
      available.attributes <- c("fit_", "unpenalized_fit_", "auc_", "accuracy_", "cor_", "aic_", 
                                "intercept_", "eval.sparsity", "precision_", "recall_", "f1_")
      if(!attribute %in% available.attributes)
      {
        stop("estimateFeatureImportance: attribute does not exist.")
      }
      
      if(!attribute %in% names(mod) | !attribute %in% names(mod.tmp))
      {
        stop("estimateFeatureImportance: attribute does not exist in the models.")
      }
      
      importance[i]           <- mod[[attribute]] - mod.tmp[[attribute]]
    }
  }
  
  
  if(plot.importance)
  {
    barplot(sort(importance), col="darkgreen", horiz=TRUE, las=2)
  }
  mod$importance_           <- importance
  return(mod)
}


#' Creates an object individual
#'
#' @description Creates an object individual
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the class vector
#' @param clf: the object containing the classifier information
#' @param ind: the indexes of the variables forming the individual could be null if we give the function a dense vector (via the coeff parameter) or if we also want to generate the individual
#' @param coeffs: the coefficients of the model, it could be a dense vector (in this case, ind need to be null), or it could be only the non zero values, or if it's null a new individual will be genrated
#' @param obj: an object to be incorporated in the model (default:NULL). We use this usually for SOTA.
#' @param res_clf: if provided information on mda etc can be found and transmitted to the model object
#' @return an individual (model) object 
#' @export
individual <- function(X, y, clf, coeffs = NULL, ind = NULL, eval.all= FALSE, signs = NULL, obj = NULL, res_clf = NULL)
{
  if(clf$params$objective != "cor") # force y as a factor if it is not the case
  {
    if(!is.factor(y)) y <- as.factor(y)
  }
  
  # if we don't create a model from a vector of variable indices
  if(is.null(ind))
  {
    if(is.null(coeffs)) # if there's no coeffs we generate them
    {
      # if the function is not specified (this is probably because the function is used in a particular learner)
      if(!is.null(clf$functions$individual_vec))
      {
        coeffs <- clf$functions$individual_vec(clf, signs = signs)  
      }else
      {
        # use the default version by computing the features
        if(clf$params$warnings) warning("individual: the coefficients can not be computed since no function is specified and no indexes either")
        return(NULL)
      }
    }
    # for BTR languages
    ind <- which(coeffs != 0)
    coeffs <- coeffs[ind]
    
  } else # if indices exist
  {
    # sanity check
    if(length(ind) == 0) 
    { 
      stop(print("individual: Model does not contain any features"))
    }
    
    if(any(is.na(ind))) 
    { 
      stop(print("individual: Model indices contain NAs!")) 
    }
    
    if(is.character(ind))
    {
      ind <- names2index(X, ind)
      if(clf$params$warnings) warning("individual: This model contains characters and not indexes. Converting ...")
    }
    ind <- as.numeric(ind)
    
    # reorder ind for easier comparison betwween models
    ind <- sort(ind)
  }
  
  individual                  <- list()
  # overall experiment information
  individual$learner          <- clf$learner
  individual$language         <- clf$params$language
  individual$objective        <- clf$params$objective
  individual$evalToFit        <- clf$params$evalToFit
  
  # information on the model itself
  individual$indices_         <- ind # vector of indices of non-zero coefficients (e.g. [19,22])
  individual$names_           <- rownames(X)[ind] # vector of names of features of non-zero coefficient.
  
  if(!is.null(obj))
  {
    individual$obj            <- obj
  }
  
  if(is.null(individual$names_))
  {
    stop("individual: no names for the model... looks like a highjack of the code")
  }
  
  if(is.null(coeffs)) # if no coefficients are set
  {
    if(is.null(clf$params$estimate_coefs)) 
    {
      clf$params$estimate_coefs = FALSE # by default if not specified
    }
    
    if(clf$params$estimate_coefs) # if we need to estimate the coefficients using a linear model
    {
      # If a multiple logistic regression model is selected to estimate the coefficents
      individual$coeffs_ <- estimateCoefficientsIndividual(X=X, y=y, ind = individual$indices_)$coefs
      if (any(is.na(individual$coeffs_)))
      {
        # This has happened when two vectors of the data are the same and the coefficient can't be computed.
        # In that case we don't take it into account and put it to zero
        individual$coeffs_[is.na(individual$coeffs_)] <- 0
        if(clf$params$warnings) warning("individual: coeff problem to be checked !") # We should associate an error code to this
      }
      
    }else
    { 
      # if we don't need to estimate coefficents using a linear model
      # no sign needed all are positive for the binary languages
      if(clf$params$language=="bin" | clf$params$language=="bininter") 
      {
        individual$coeffs_  <- rep(1, length(individual$indices_))
      }else # for the oather languages get the sign
      {
        if(!is.null(clf$coeffs_))
        {
          if(nrow(X) != length(clf$coeffs_))
          {
            stop("individual: there is a problem with the coeffs in the clf. Does not match with the number of features in X")
          }
          individual$coeffs_  <- clf$coeffs_[ind] # get coeffs if comuted at the begining
        }else
        {
          individual$coeffs_  <- getSign(X = X[ind,], y = y, clf = clf)  # +1, -1 for the ternary learning  
        }
      }
    }
  } else # else if coefficients are set add them to the individual
  { 
    individual$coeffs_        <- coeffs
  }
  
  # Setting empty values
  individual$fit_             <- NA
  individual$unpenalized_fit_ <- NA
  individual$auc_             <- NA
  individual$accuracy_        <- NA
  individual$cor_             <- NA
  individual$aic_             <- NA
  
  # if intercept is not set in the global parameters
  if(clf$params$intercept == "NULL")
  {
    individual$intercept_     <- NA
  } else 
  {
    individual$intercept_     <- clf$params$intercept
    individual$sign_          <- ">" # for a given intercept set the sign by default
  }
  
  individual$eval.sparsity    <- length(ind)
  
  if(eval.all)
  {
    individual                <- evaluateModel(mod = individual, 
                                               X = X, 
                                               y = y, 
                                               clf = clf, 
                                               eval.all = TRUE, 
                                               force.re.evaluation = TRUE)
  }
  
  # add more information to the model (mda etc.) 
  if(!is.null(res_clf))
  {
    # check existance of fip data
    if(!is.null(res_clf$crossVal$fip))
    {
      feature.importance.cv <- rowMeans(res_clf$crossVal$fip$mda, na.rm = TRUE)
      feature.prevalence.cv <- res_clf$crossVal$fip$fpf / ncol(res_clf$crossVal$fip$mda) # normalize %
      if(!is.null(res_clf$classifier$fip$mda))
      {
        feature.importance  <- res_clf$classifier$fip$mda  
      }else
      {
        feature.importance  <- rep(NA, length(clf$data$features))
        names(feature.importance) <- names(clf$data$features)
      }
      
      
      # MDA generalization
      individual$mda.cv_ <- rep(0, length(individual$names_))
      names(individual$mda.cv_) <- individual$names_
      ind.features <- intersect(individual$names_, names(feature.importance.cv))
      individual$mda.cv_[ind.features] <- feature.importance.cv[ind.features]
      
      # prevalence in top models in the folds
      individual$prev.cv_ <- rep(0, length(individual$names_))
      names(individual$prev.cv_) <- individual$names_
      ind.features <- intersect(individual$names_, names(feature.prevalence.cv))
      individual$prev.cv_[ind.features] <- feature.prevalence.cv[ind.features]
      
      # MDA empirical
      individual$mda_ <- rep(0, length(individual$names_))
      names(individual$mda_) <- individual$names_
      ind.features <- intersect(individual$names_, names(feature.importance))
      individual$mda_[ind.features] <- feature.importance[ind.features]
      
    } # end fip existance testing
  } # end res_clf existance testing
  
  return(individual)
}


#' @description Evaluates an entire population of models, that be predomics objects or individuals
#' 
#' @import foreach
#' @title evaluatePopulation
#' @name evaluatePopulation
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the class vector
#' @param clf: the object containing the classifier information
#' @param pop: the population of models to be evaluated
#' @param eval.all: should the function evaluate all the scores for each of the models (default:FALSE)
#' @param force.re.evaluation: re-evaluate all the scores even if they exist for each of the models (default:FALSE)
#' @param estim.feat.importance: evaluate the importance in the model object for each of the models (default:FALSE)
#' @param mode: A choice from c("train", "test") indicates wether we wish to learn the threthold 
#' of each of the models (default:"train") or not "test" for the c("terinter","bininter","ratio") languages
#' @param delete.null.models: should null indivuals be deleted (default:TRUE)
#' @param lfolds: compute evaluation in crossval (default:NULL)
#' @return an individual object 
#' @export
evaluatePopulation <- function(X, y, clf, pop, eval.all = FALSE, 
                               force.re.evaluation = FALSE, 
                               estim.feat.importance = FALSE, 
                               mode = "train", 
                               delete.null.models = TRUE, 
                               lfolds = NULL) 
{
  # test the classifier object
  if(!isClf(clf))
  {
    stop("fit: please provide a valid classifier object!")
  }
  
  check.X_y_w(X, y, w=NULL)
  
  # clean the null models
  if(delete.null.models)
  {
    pop <- pop[!sapply(pop, is.null)]
  }
  
  # if crossval call this in recursive
  if(!is.null(lfolds))
  {
    pop.lfolds <- list()
    for(f in 1:length(lfolds))
    {
      if(mode == "train")
      {
        pop.lfolds[[f]] <- evaluatePopulation(X[,-lfolds[[f]]],
                                              y[-lfolds[[f]]],
                                              clf,
                                              pop,
                                              eval.all = eval.all,
                                              force.re.evaluation = force.re.evaluation,
                                              estim.feat.importance = estim.feat.importance,
                                              mode = mode,
                                              delete.null.models = delete.null.models,
                                              lfolds = NULL)
      }else # test
      {
        pop.lfolds[[f]] <- evaluatePopulation(X[,lfolds[[f]]],
                                              y[lfolds[[f]]],
                                              clf,
                                              pop,
                                              eval.all = eval.all,
                                              force.re.evaluation = force.re.evaluation,
                                              estim.feat.importance = estim.feat.importance,
                                              mode = mode,
                                              delete.null.models = delete.null.models,
                                              lfolds = NULL)
      }
    }
    names(pop.lfolds) <- names(lfolds)
    return(pop.lfolds)
  }
  
  # otherwise we continue to evaluate normally the population
  
  # # Parallel computing
  # res <- list()
  # if(clf$params$parallel.local)
  # {
  #   res <- foreach(i = 1:length(pop))  %dorng% 
  #   { 
  #     mod <- pop[[i]]
  #     if(!is.null(mod))
  #     {
  #       evaluateModel(mod = mod,
  #                     X = X, 
  #                     y = y, 
  #                     clf = clf, 
  #                     eval.all = eval.all, 
  #                     force.re.evaluation = force.re.evaluation, 
  #                     estim.feat.importance = estim.feat.importance, 
  #                     mode = mode)
  #     } # end else existance pop
  #   } # end foreach loop
  # } else
  # {
  #   res <- list()
  #   for (i in 1:length(pop)) # for all the individuals in the population
  #   {
  #     mod <- pop[[i]]
  #     if(!is.null(mod))
  #     {
  #       res[[i]] <- evaluateModel(mod = mod, 
  #                                 X = X, 
  #                                 y = y, 
  #                                 clf = clf, 
  #                                 eval.all = eval.all, 
  #                                 force.re.evaluation = force.re.evaluation, 
  #                                 estim.feat.importance = estim.feat.importance, 
  #                                 mode = mode)
  #     } # end else existance pop
  #   } # end for loop
  # } # end else //
  
  res <- list()
  for (i in 1:length(pop)) # for all the individuals in the population
  {
    mod <- pop[[i]]
    if(!is.null(mod))
    {
      res[[i]] <- evaluateModel(mod = mod, 
                                X = X, 
                                y = y, 
                                clf = clf, 
                                eval.all = eval.all, 
                                force.re.evaluation = force.re.evaluation, 
                                estim.feat.importance = estim.feat.importance, 
                                mode = mode)
      #print(i)
    } # end else existance pop
  } # end for loop
  
  # clean population after evaluation as well
  if(delete.null.models)
  {
    res <- cleanPopulation(pop = res, clf = clf)
  }
  
  return(res)
}


#' sortPopulation
#' 
#' @description Sort a population according to a given attribute (evalToOrder)
#' @param pop: a population (list) of evaluated predomics objects
#' @param evalToOrder: the attribute to be used in the sorting (default:fit_)
#' @param decreasing: whether the sorting should be be decreasing or not (default:decreasing)
#' @return a sorted population of predomics objects
#' @export
sortPopulation <- function(pop, evalToOrder = "fit_", decreasing = TRUE)
{
  
  if(!isPopulation(pop))
  {
    warning("sortPopulation: the population object is not a a valid one, returning NULL")
    return(NULL)
  }
  
  # # OVERRIDE  
  # # for correlation use the standard error of the mean and reverse
  # if(pop[[1]]$objective == "cor")
  # {
  #   #
  #   evalToOrder   <- "ser_"
  #   decreasing    <- TRUE
  # }
  
  evaluation <- populationGet_X(element2get = evalToOrder, toVec = TRUE, na.rm = TRUE)(pop)
  
  if(length(evaluation) > 0)
  {
    ord <- order(evaluation, decreasing = decreasing)  
  }else
  {
    return(NULL)
  }
  
  return(pop[ord])
}



#' Computes the ^y score of the model as a ratio
#'
#' @description Computes the ^y score of the model as a ratio
#' @param class_1_score: the sum score for the features of class 1
#' @param class_2_score: the sum score for the features of class 2
#' @param epsilon: is a very small value that will would avoid Inf values in the ratio. This can be either specified in the when setting the classifier and if not specified will be set as the minimum number of the machine (e.g. 2.23e-308). Caution this should be adapted when working with other types of data.
#' @return a vector containing the predicted ^y score for each observation
#' @export
scoreRatio <- function(class_1_score, class_2_score, epsilon = NULL)
{
  class_1_score <- as.numeric(class_1_score)
  class_2_score <- as.numeric(class_2_score)
  
  if(epsilon=="NULL" | is.null(epsilon))
  {
    warning("Please provide epsilon. Inf values can be produced if not provided. Setting to default 2.23e-308")
    epsilon <- .Machine$double.xmin
  } else
  {
    # add the epsilon to the denominator
    class_2_score <- class_2_score + epsilon
    # add the epsilon to the nominator for equal treatment
    class_1_score <- class_1_score + epsilon
  }
  
  # initialize everything as NA. It is easier to spot issues when debugging
  res <- rep(NA, length(class_1_score))
  for(i in 1:length(res))
  {
    res[i] <- class_1_score[i] / class_2_score[i]
  }
  return(res)
}


################################################################
# GETTERS, SETTERS ...
################################################################

#' Get the fitting score of an individual object
#'
#' @description Get the fitting score of an individual object.
#' @param individual : an individual object
#' @return a fitting score
getFitIndividual <- function(individual)
{
  #return(ifelse(!is.na(individual), individual$fit_, {print("TEST"); NA}))
  return(individual$fit_)
}

#' Get the index of the features in a given individual
#'
#' @description Get the indices of the features used in the individuals
#' @param individual : an individual object
#' @return the indices of the features
getIndicesIndividual <- function(individual)
{
  return(individual$indices_)
}

#' Get the indices of the features used in a population of individuals
#'
#' @description Get the indices of the features used in a population of individuals
#' @param pop : a list of individuals
#' @return a matrix of indices (rows), and individuals (cols)
getIndicesPopulation <- function(pop)
{
  return(sapply(pop, getIndicesIndividual))
}


getAccuracyIndividual <- function(mod){ return(mod$accuracy_) }

getAccuracyPopulation <- function(pop)
{
  return(sapply(pop, getAccuracyIndividual))
}



#' Get the models from a classifier result for each k-sparsity
#'
#' @description Get the N best models from a classifier result for each k-sparsity.
#' @param obj: the classifier result output from the function fit. This can also be a ModelCollection or Population object
#' @param significance: if TRUE, (default:FALSE) a statistical test will be applied to find the lowest threshold that will delimit the window
#' of the best models. If FALSE, the models will be selected according to the rest of the criteria.
#' @param by.k.sparsity: if TRUE (default:TRUE), the filtering will be performed for each sparsity level
#' @param k.penalty: (default:0), it will penalize the models with large sparsity if different, when by.k.sparsity is set to TRUE
#' @param n.best: the number of best models to be returned for each sparsity if by.k.sparsity is set to TRUE or for the whole population 
#' otherwise (default:5).
#' @param nbest: the number of best models we wish to get from the population, per each sparsity or not. If there are less best models then this
#' number, less will be returned
#' @param single.best: if TRUE, this will return the best model of all (default:FALSE) and the n.best will be set to 1. 
#' @param single.best.cv: if single.best is TRUE, we could chose the best model based on data from cross validation (default:TRUE) and in this 
#' case obj should be an experiment or from empirical results not in CV.
#' @param single.best.k: if single.best is TRUE, we could chose the best model of a given sparsity that is specified by a number here. 
#' If this value is specified (default:NULL), then this will de-actvate single.best.cv.
#' @param max.min.prevalence: if TRUE (default:FALSE), the best models will be selected based on their performance but also on the prevalence of 
#' the features that compose it.
#' @param X: the dataset to be learned (default:NULL). This is neeeded when max.min.prevalence is set to TRUE.
#' @param verbose: provide more information about execution (default = FALSE)
#' @param evalToOrder: which attribute of the model object should we use to order the models and select them (default:fit_)
#' @param return.population: if set to TRUE (default:FALSE), the result will be send as a population of models
#' @param unique.control: if set to TRUE (default:TRUZ), we correct the population so that no dupplication of models takes place
#' @return a list of model objects or a model when it is a single one or a model collection
#' @export
getNBestModels <- function(obj, 
                           significance = FALSE, 
                           by.k.sparsity = TRUE,
                           k.penalty = 0,
                           n.best = 5,
                           single.best = FALSE,
                           single.best.cv = TRUE,
                           single.best.k = NULL,
                           max.min.prevalence = FALSE,
                           X = NULL,
                           verbose = FALSE, 
                           evalToOrder = "fit_",
                           return.population = FALSE,
                           unique.control = TRUE
)
{
  
  # sanity check
  if(!isExperiment(obj) & !isModelCollection(obj) & !isPopulation(obj))
  {
    warning("getNBestModels: please provide a valid experiment, modelCollection or population object ... returning NULL")
    return(NULL)
  }
  
  # if an experiment
  if(isExperiment(obj = obj)) 
  {
    if(verbose) print(paste0("getNBestModels: the object is an experiment"))
    mc              <- obj$classifier$models
  }
  
  # convert to a model collection if it is not
  if(isPopulation(obj = obj)) # if a population
  {
    if(verbose) print(paste0("getNBestModels: the object is an population of models"))
    mc              <- listOfModels2ModelCollection(pop = obj)
  }
  # if a modelCollection
  if(isModelCollection(obj = obj)) 
  {
    if(verbose) print(paste0("getNBestModels: the object is a model collection"))
    mc              <- obj
  }
  
  if(!isModelCollection(mc))
  {
    warning("getNBestModels: the object is not a valid model collection. Returning empty handed.")
    return(NULL)
  }
  
  if(single.best) 
  {
    n.best          <- 1
    if(verbose) print(paste0("getNBestModels: single best"))
  }
  
  # set up switch variables that are no needed to parameterize but that we could in the future
  penalty_by_kfold <- FALSE
  
  if(by.k.sparsity | single.best)
  {
    ####################################################################
    # # if by.k.sparsity
    ####################################################################
    
    if(verbose) print(paste0("getNBestModels: by k sparsity"))
    res <- list()
    pop.valids <- c()
    # for each k_sparsity
    for(i in 1:length(mc))
    {
      if(unique.control)
      {
        pop           <- unique(mc[[i]])
      }else
      {
        pop           <- mc[[i]]  
      }
      
      if(verbose) print(paste0("getNBestModels: the population of sparsity ", i," contains ", length(pop), " models"))
      
      if(significance)
      {
        # restrict the population to the confidence interval
        pop         <- selectBestPopulation(pop = pop, score = evalToOrder, p = 0.05, k_penalty = k.penalty)
        if(verbose) print(paste0("getNBestModels: after significance selection with k.penalty ", k.penalty," it contains ", length(pop), " models"))
      }
      
      # if we wish to select best models with max min prevalence
      if(max.min.prevalence)
      {
        if(!is.null(X))
        {
          eval      <- populationGet_X(element2get = evalToOrder, toVec = TRUE, na.rm = FALSE)(pop)
          best.eval <- max(eval, na.rm = T)
          # epsilon is used to be able to select a window of best models
          epsilon   <- sqrt(best.eval * (1 - best.eval) / ncol(X))
          pop       <- pop[eval >= (best.eval - epsilon) & !is.na(eval)]
          pop       <- getMaxMinPrevalenceModel(pop = pop, X = X, selected = 0)
          if(verbose) print(paste0("getNBestModels: after max.min.prevalence it contains ", length(pop), " models"))
        }
      }
      
      pop           <- sortPopulation(pop = pop, evalToOrder = evalToOrder)
      pop           <- pop[1:min(n.best, length(pop))]
      if(verbose) print(paste0("getNBestModels: the final population contains ", length(pop), " models"))
      res[[i]]      <- pop
      # mark valididity
      pop.valids <- c(pop.valids, isPopulation(pop))
      
    } # end for loop
    names(pop.valids) <- names(mc)
    names(res)      <- names(mc)[pop.valids]
    
    mc <- mc[pop.valids]
    
    if(!isModelCollection(mc))
    {
      warning("digestModelCollection: after treating the mc object no result is available. Returning NULL")
      return(NULL)
    }
    
    if(single.best)
    {
      single.best.cv <- FALSE
      # set best model type switch
      if(isExperiment(obj) & !myAssertNotNullNorNa(single.best.k))
      {
        if(!is.null(obj$crossVal))
        {
          single.best.cv <- TRUE  
        }
      }
      
      k_catalogue <- paste0("k_",obj$classifier$params$sparsity)
      spar        <- populationGet_X("eval.sparsity", toVec = TRUE, na.rm = FALSE)(modelCollectionToPopulation(res))
      
      # Here we are left with two options
      if(single.best.cv)
      {
        # get the best model from crossval information
        if(verbose) print(paste0("getNBestModels: single.best.cv mode ... returning the best model"))
        if(obj$classifier$params$objective == "auc" & !(evalToOrder == "accuracy_" | evalToOrder == "auc_"))
        {
          evalToOrder <- "accuracy_"
        }
        if(obj$classifier$params$objective == "cor" & !(evalToOrder == "cor_"))
        {
          evalToOrder <- "cor_"
        }
        
        # Abreviations for the results
        key <- data.frame(real=c("auc_","accuracy_","cor_"), abbrv=c("auc","acc","cor")); rownames(key) <- key$real
        emp.name <- paste("empirical", as.character(key[evalToOrder,]$abbrv), sep=".")
        gen.name <- paste("generalization", as.character(key[evalToOrder,]$abbrv), sep=".")
        # for each classifier
        emp.data <- obj$crossVal$scores[[emp.name]][k_catalogue, ]
        gen.data <- obj$crossVal$scores[[gen.name]][k_catalogue, ]
        # plot for debug
        # par(mfrow=c(2,1)); image(as.matrix(t(emp.data))); image(as.matrix(t(gen.data)))
        
        if(!is.null(emp.data) & !is.null(gen.data))
        {
          # compute the penalty data
          emp.data.penalty          <- emp.data
          k                         <- as.numeric(gsub("k_","",k_catalogue))
          
          
          # if we want to compute the penalty by k_fold
          if(penalty_by_kfold)
          {
            for(j in 1:nrow(emp.data))
            {
              emp.data.penalty[j,]  <- emp.data[j,] - penalty[i] * k[j]
            }
            # select the k_sparse for each k_fold
            ind <- apply(emp.data.penalty, 2, which.max)
            k_sparse <- rownames(emp.data.penalty)[ind]
            
            best_empirical <- diag(as.matrix(emp.data[ind,]))
            best_generalization <- diag(as.matrix(gen.data[ind,]))
            
          }else # otherwise we compute a meaned penalty
          {
            mean_score <- rowMeans(emp.data.penalty, na.rm = TRUE)
            mean_score_penalty  <- mean_score - k.penalty * k
            #plot(mean_score, mean_score_penalty, ylim=c(0,1),xlim=c(0.5,1))
            
            # make sure to be in the space of available sparsity models in the emperical models
            ind <- which.max(mean_score_penalty[names(mc)])
            k_sparse <- rep(rownames(emp.data.penalty[names(mc),])[ind], ncol(emp.data))
            best_empirical <- as.numeric(emp.data[names(mc),][ind,])
            best_generalization <- as.numeric(gen.data[names(mc),][ind,])
            
            # if no values are found in empirical
            if(length(ind) == 0)
            {
              best_empirical <- logical(0)
              best_generalization <- logical(0)
            }
            
            # => TEST if(all(is.na(best_empirical))) best_empirical <- logical(0)
            # => TEST if(all(is.na(best_generalization))) best_generalization <- logical(0)
            # plot(best_empirical, best_generalization, ylim=c(0.5,1),xlim=c(0.5,1))
            # boxplot(list(best_empirical,best_generalization), ylim=ylim)
          }
          res.k.cv <- data.frame(best_empirical, best_generalization, k_sparse)
        }
        else
        {
          res.k.cv <- data.frame(best_empirical = NA, best_generalization = NA)[-1,]
        }
        best.k <- as.numeric(gsub("k_","",as.character(unique(res.k.cv$k_sparse))))
        
      }else
      {
        # get the best model from empirical information
        if(verbose) print(paste0("getNBestModels: single.best mode ... returning the best model"))
        
        eval <- NULL
        for(i in 1:length(res))
        {
          eval.i <- populationGet_X(evalToOrder)(res[[i]])
          if(length(eval.i) > 0)
          {
            eval <- c(eval, eval.i)
          }else
          {
            eval <- c(eval, NA)
          }
        }
        # eval <- unlist(lapply(res, function(x){populationGet_X(evalToOrder)(x)[[1]]}))
        # apply the k_penalty
        eval <- eval - (k.penalty * spar)
        
        best.k <- as.numeric(gsub("k_","",names(eval[which.max(eval)])))
      }
      
      # select the best model for a given k
      if(!myAssertNotNullNorNa(single.best.k))
      {
        # set from above if it does not exist
        single.best.k <- best.k
      }else
      {
        if(single.best.k == 0)
        {
          # this is a special case, and means that there will not be a selection but the maximum number of variables will be taken into account
          #k   <- as.numeric(gsub("k_","",k_catalogue))
          single.best.k <- max(spar)
        }
      }
      
      k_spar <- paste0("k_", single.best.k)
      
      # otherwise when single.best.k is specified this will be the first choice
      if(length(single.best.k) == 1 & is.numeric(single.best.k))
      {
        if(k_spar %in% names(mc)) # check if we have it in the model collection
        {
          if(verbose) print(paste0("getNBestModels: single.best.k mode with k_spar ", k_spar, " returning the best model"))
          pop     <- mc[[k_spar]]
          eval    <- populationGet_X(element2get = evalToOrder, toVec = TRUE, na.rm = FALSE)(pop)
          mod     <- pop[[which.max(eval)]]
          if(return.population)
          {
            return(list(mod))
          }else
          {
            return(mod)  
          }
        }else # not found
        {
          if(verbose) print(paste0("getNBestModels: single.best.k mode with k_spar ", k_spar, " not found in the results"))
          if(return.population)
          {
            return(list(NA))
          }else
          {
            return(NA)  
          }
        }
      }else
      {
        print(paste0("getNBestModels: single.best.k mode but ",k_spar, " is not found in the model collection. Executing the default settings."))
      }
    } # end of single.best.k 
    
    if(return.population)
    {
      if(verbose) print(paste0("getNBestModels: returning a population of models"))
      # Transform the model collection onto a population
      return(modelCollectionToPopulation(mod.collection = res))
    }else
    {
      if(verbose) print(paste0("getNBestModels: returning a model collection"))
      # a model collection
      return(res)
    }
    
  }else 
  {
    ####################################################################
    # # else not by.k.sparsity
    ####################################################################
    if(verbose) print(paste0("getNBestModels: regardless of k sparsity"))
    # first convert the pop
    if(unique.control)
    {
      pop             <- unique(modelCollectionToPopulation(mc))
    }else
    {
      pop             <- modelCollectionToPopulation(mc)  
    }
    
    if(verbose) print(paste0("getNBestModels: the population with all sparsities contains ", length(pop), " models"))
    
    if(significance)
    {
      # restrict the population to the confidence interval
      pop           <- selectBestPopulation(pop = pop, score = evalToOrder, p = 0.05, k_penalty = k.penalty)
      if(verbose) print(paste0("getNBestModels: after significance selection with k.penalty ", k.penalty," it contains ", length(pop), " models"))
    }
    
    # if we wish to select best models with max min prevalence
    if(max.min.prevalence)
    {
      if(!is.null(X))
      {
        eval        <- populationGet_X(element2get = evalToOrder, toVec = TRUE, na.rm = FALSE)(pop)
        k           <- populationGet_X(element2get = "eval.sparsity", toVec = TRUE, na.rm = TRUE)(pop)
        eval.penalty<- eval - (k*k.penalty)
        best.eval   <- max(eval.penalty)
        epsilon     <- sqrt(best.eval * (1 - best.eval) / ncol(X))
        pop         <- pop[eval.penalty >= (best.eval - epsilon)]
        pop         <- getMaxMinPrevalenceModel(pop = pop, X = X, selected = 0)
        if(verbose) print(paste0("getNBestModels: after max.min.prevalence it contains ", length(pop), " models"))
      }
    } # end max.min.prevalence
    
  } # end by.k.sparsity ifelse
  
  
  if(return.population)
  {
    if(verbose) print(paste0("getNBestModels: returning a population of models"))
    return(pop)
  }else
  {
    if(verbose) print(paste0("getNBestModels: returning a model collection"))
    return(listOfModels2ModelCollection(pop))
  }
} 



getTheBestIndividual <- function(pop, evalToFit = "fit_")
{
  if(length(pop) > 0)
  {
    pop <- sortPopulation(pop,evalToOrder = evalToFit)
    
    return(pop[[1]])
    
  } else
  {
    return(NULL)
  }
}



#' Get the best model from a classifier result
#'
#' @description Gets a given attribute from a population of predomics objects
#' @param element2get: the name of the attribute to get
#' @param toVec: should the results be unlisted (default:TRUE)
#' @param na.rm: delete the elements that are NA (default) when returning tovec
#' @return a vector of attributes 
#' @export
populationGet_X <- function(element2get, toVec = TRUE, na.rm = TRUE)
{
  # custom function
  func <- function(pop)
  {
    # sanity check
    if(length(pop) > 0)
    {
      res <- lapply(pop, function(indiv) 
        if(!is.list(indiv)) 
        {
          NA
        }else
        {
          indiv[[element2get]]
        }
      )
    } else 
    {
      res <- NA
    }
    
    if(toVec) # return vector
    {
      if(na.rm) 
      {
        ind.na <- which(!is.na(res))
        return(unlist(res[ind.na]))  
      }else
      {
        return(unlist(res))
      }
      
    } else # return list 
    {
      if(na.rm)
      {
        ind.na <- which(!is.na(res))
        return(res[ind.na])
      }else
      {
        return(res)
      }
    }
  }
  
  return(func)
}


#' Set models with a given liist of objects
#'
#' @description Sets a given attribute to the objects of the a given population
#' @param element2set: the name of the attribute to set
#' @param listwithelements: the list containing the elements to add
#' @return an updated population
#' @export
populationSet_X <- function(pop, element2set = NULL, listwithelements = NULL)
{
  # sanity check
  if(is.null(pop))
  {
    stop(paste("populationSet_X: please provide a non empty population"))
  }
  
  if(!is.list(pop))
  {
    stop(paste("populationSet_X: please provide a list object as a population"))
  }
  
  # check the type of elements to set
  if(is.list(listwithelements))
  {
    existing <- populationGet_X(element2get = element2set, toVec = FALSE, na.rm=TRUE)(pop)  
  }else
  {
    existing <- populationGet_X(element2get = element2set, toVec = TRUE, na.rm=TRUE)(pop)
  }
  
  if(length(existing) > 0)
  {
    print(paste("This attribute exists and",length(existing),"are found"))
  }
  
  if(length(listwithelements)==1)
  {
    if(length(existing) > 0)
    {
      if(class(existing) != class(listwithelements)) 
      {
        stop(paste("populationSet_X: please make sure you provide the same type as existing for",element2set))
      }
    }
  }else
  {
    if(length(listwithelements) != length(pop))
    {
      stop(paste("populationSet_X: please make sure the elements to add are the same number as elements in the population"))
    }
  }
  
  # for each element of the population
  for(i in 1:length(pop))
  {
    if(length(listwithelements)==1)
    {
      if(!is.null(pop[[i]]))
      {
        pop[[i]][element2set] <- listwithelements
      }else
      {
        warning(paste("populationSet_X: element",i,"of the list is null ... skiping"))
      }
    }else
    {
      if(is.list(listwithelements))
      {
        pop[[i]][[element2set]] <- listwithelements[[i]]
      }else
      {
        pop[[i]][[element2set]] <- listwithelements[i]
      }
    }
  } # end for loop pop
  
  return(pop)
}

#' Get the fitting score of a list a models
#'
#' @description Get the fitting score of a list a models.
#' @param pop : a list of models
#' @return a vector of fitting scores
getFitModels <- function(pop)
{
  res <- rep(NA,length(pop))
  for(i in 1: length(res))
  {
    res[i] <- getFitModel(pop[[i]])
  }
  return(res)
}

#' Get the fitting score of a model object
#'
#' @description Get the fitting score of a model object.
#' @param mod : a model object
#' @return a fitting score
getFitModel <- function(mod)
{
  return(mod$fit_)
}

#' Get the fitting score of a list of individuals
#'
#' @description Get the fitting score of a list of individuals.
#' @param pop : a list of individuals
#' @return a vector of fitting scores
getFitPopulation <- function(pop)
{
  return(sapply(pop, getFitIndividual))
}



################################################################
# CHECKERS
################################################################

check.X_y_w <- function(X, y, w=NULL) {
  if (dim(X)[2] != length(y)) 
  { 
    stop("check.X_y_w: dimension of X and y is not coherent") 
  }
  if (!is.null(w)) 
  {
    if (dim(X)[1] != length(w)) 
    { 
      stop("check.X_y_w: dimension of X and w is not coherent") 
    }
  }
}

check.tX_y_w <- function(X,y,w=NULL) 
{
  if (dim(X)[1] != length(y)) 
  { 
    stop("check.txyw: dimension of X and y is not coherent") 
  }
  if (!is.null(w)) 
  {
    if (dim(X)[2] != length(w)) 
    { 
      stop("check.txyw: dimension of X and w is not coherent") 
    }
  }
}


################################################################
# CONVERTORS
################################################################
# natts is the number of attributes,

#' Transform the model object onto dense format (long) one
#'
#' @description Builds a model object based on model that is in the dense (long) format.
#' @param natts: the number of attributes
#' @param mod: a predomics model object
#' @return a dense (long) format model
#' @export
modelToDenseVec <- function(natts, mod) 
{
  # test model validity
  if(!isModel(obj = mod))
  {
    stop("modelToDenseVec: the model object is not valid.")
  }
  
  # test the number of attributes
  if(class(natts)!="integer" | length(natts)!=1)
  {
    stop("modelToDenseVec: please make sure the natts attribute is an integer number of features.")
  }
  
  # initialize the v
  v <- rep(0, natts)
  
  if(isModelSota(mod))
  {
    mod$coeffs_ <- rep("*", length(mod$indices_))
    names(mod$coeffs_) <- names(mod$indices_)
  }
  
  for(i in 1:length(mod$indices_)) 
  {
    v[[ mod$indices_[[i]] ]] <- mod$coeffs_[[i]]
  }
  return(v)
}


#' denseVecToModel
#'
#' @description Builds a model object based on model that is in the dense (long) format.
#' @param X: dataset
#' @param y: labels
#' @param v: A vector of coeffs (example v=c(0.0,1.0,0.0,-1.0))
#' @param clf: classifier information
#' @param eval.all: If TRUE the fitting of the function and intercept will be computed
#' @param obj: an object model to add to the model (default:NULL)
#' @return an model object
denseVecToModel <- function(X, y, v, clf, eval.all=FALSE, obj = NULL) 
{
  if(is.null(v))
  {
    return(NULL)
  }
  res <- individual(X, y, clf, coeffs = v[which(v != 0)], ind = which(v != 0), eval.all = eval.all, obj = obj)
  return(res)
}

#' sparseVecToModel
#'
#' @description Builds a model object based on model that is in the sparse (short) format.
#' @param X: dataset
#' @param y: labels
#' @param v: A vector of indexes (example v=c(1,11))
#' @param clf: classifier information
#' @param eval.all: Should the model be evaluated (default:FALSE)
#' @param obj: an object model to add to the model (default:NULL)
#' @return an model object
sparseVecToModel <- function(X, y, v, clf, eval.all=FALSE, obj = NULL) 
{
  if(is.null(v))
  {
    return(NULL)
  }
  res <- individual(X, y, clf, ind = v, eval.all = eval.all, obj = obj)
  return(res)
}

#' Builds a model object from a list of vector coefficients
#'
#' @description Builds a model object from a list of vector coefficients.
#' @import snow
#' @param X: dataset
#' @param y: labels
#' @param clf: classifier
#' @param v: list of vectors of coeffs. For example, v=list( c(0.0,1.0,0.0,-1.0) , c(1.0,1.0,0.0,0.0) , c(0.0,1.0,1.0,-1.0) )
#' @param lobj: a list of objects to add as elements in the model objects if not null (default:NULL)
#' @return an model object
#' @export
listOfDenseVecToListOfModels <- function(X, y, clf, v, lobj = NULL) 
{
  if(is.null(v))
  {
    return(NULL)
  }
  
  if(!is.null(lobj))
  {
    if(!is.list(lobj))
    {
      stop("listOfDenseVecToListOfModels: lobj should be a list of objects.")
    }
    
    if(length(lobj) != length(v))
    {
      stop("listOfDenseVecToListOfModels: lobj should be a list the same length as v")
    }
  }
  
  pop <- list()
  for(i in 1:length(v))
  {
    model <- v[[i]]
    check.X_y_w(X, y, model)
    pop[[i]] <- denseVecToModel(X, y, model, clf, eval.all=TRUE, obj = lobj[[i]])
  }
  
  return(pop)
}


#' listOfSparseVecToListOfModels
#'
#' @description Converts an list of "SparseVec" objects onto a list of predomics objects
#' @import snow
#' @param X: dataset
#' @param y: labels
#' @param clf: classifier
#' @param v: list of vectors of coeffs. For example, v=list( c(0.0,1.0,0.0,-1.0) , c(1.0,1.0,0.0,0.0) , c(0.0,1.0,1.0,-1.0) )
#' @param lobj: a list of objects to add as elements in the model objects if not null (default:NULL)
#' @param eval.all: evaluate population (default:FALSE)
#' @return an model object
#' @export
listOfSparseVecToListOfModels <- function(X, y, clf, v, lobj = NULL, eval.all = FALSE) 
{
  
  if(!is.null(lobj))
  {
    if(!is.list(lobj))
    {
      stop("listOfDenseVecToListOfModels: lobj should be a list of objects.")
    }
    
    if(length(lobj) != length(v))
    {
      stop("listOfDenseVecToListOfModels: lobj should be a list the same length as v")
    }
  }
  
  if(length(v) == 0)
  {
    if(clf$params$warnings) warning("listOfSparseVecToListOfModels: empty list returning NULL")
    return(NULL)
  }
  
  pop <- list()
  for(i in 1:length(v))
  {
    model <- v[[i]]
    pop[[i]] <- sparseVecToModel(X, y, model, clf, eval.all = eval.all, obj = lobj[[i]])
  }
  
  return(pop)
}


#' Builds a list of dense vector coefficients from a list of models
#' 
#' @param clf: classifier
#' @param X: dataset
#' @param y: labels
#' @param list.models: list of models
#' @return a list of dense vectors of coefficient
#' @export
listOfModelsToListOfDenseVec <- function(clf, X, y, list.models) 
{
  if(!isPopulation(obj = list.models))
  {
    stop("listOfModelsToListOfDenseVec: Please specify a population of model objects")
  }
  
  res <- list()
  for(i in 1:length(list.models))
  {
    res[[i]] <- modelToDenseVec(nrow(X), list.models[[i]])
  }
  return(res)
}


#' Builds a list of sparse vector coefficients from a list of models
#' 
#' @param list.models: list of models
#' @return a list of dense vectors of coefficient
#' @export
listOfModelsToListOfSparseVec <- function(list.models) 
{
  if(!isPopulation(obj = list.models))
  {
    stop("listOfModelsToListOfSparseVec: Please specify a population of model objects")
  }
  
  res <- populationGet_X(element2get = "indices_", toVec = FALSE, na.rm = FALSE)(list.models)
  return(res)
}


#' Builds a list of dense vector coefficients from a list of models
#' 
#' @param clf: classifier
#' @param X: dataset
#' @param y: labels
#' @param v: list of dense vectors
#' @return a model collection
listOfDenseVecToModelCollection <- function(clf, X, y, v) 
{
  # convert the list of dense vectors to a population of models
  pop <- listOfDenseVecToListOfModels(X = X, y = y, clf = clf, v = v)
  # convert it to a model collection
  mc <- listOfModels2ModelCollection(pop = pop)
  return(mc)
}


#' names2index
#'
#' @description Transforms feature names feature indexes
#' @param X: the dataset
#' @param var.names: the feature names vector
#' @return the index of the features
names2index <- function(X, var.names)
{ 
  if(class(var.ind)!="character" & class(var.ind)!="factor")
  {
    stop("index2names: feature names should be character")
  }
  index <- match(var.names, rownames(X))
  if (any(is.na(index))) 
  {
    stop("index2names: Features not found!")
  }
  return(index)
}

#' index2names
#'
#' @description Transforms feature indexes into feature names
#' @param X: the dataset
#' @param var.ind: the feature index vector
#' @return the names of the features
index2names <- function(X, var.ind)
{
  if(class(var.ind) != "numeric")
  {
    stop("index2names: indexes should be numeric")
  }
  
  if(max(var.ind) > nrow(X))
  {
    stop("index2names: index is out of the bounds")
  }
  
  return(rownames(X)[var.ind]) 
}


#' updateModelIndex
#'
#' @description Update the index of a model objectn.
#' @param obj: the object is a model
#' @param features: the list of features which overrides the clf$data$features if this exists.
#' @return the same object type as input, but updated
updateModelIndex <- function(obj, features = NULL)
{
  if(!(isModel(obj)))
  {
    stop("updateIndexes: please provide a valid model object")
  }
  
  if(is.null(features))
  {
    stop("updateIndexes: the features object is missing, which is necessary for this function.")
  }
  
  obj$indices_ <- match(obj$names_, features)
  
  if(any(is.na(obj$indices_)))
  {
    warning("Some indices are not found.")
  }
  
  return(obj)
}


#' updateObjectIndex
#'
#' @description Update the index of a model, population, or modelCollection.
#' @param obj: the object can be a model, population, or modelCollection
#' @param features: the list of features which overrides the clf$data$features if this exists.
#' @return an the same object type as input, but updated
#' @export
updateObjectIndex <- function(obj, features = NULL)
{
  if(!(isModelCollection(obj) | isModel(obj) | isPopulation(obj)))
  {
    warning("updateIndexes: please provide a model collection or a population or a single model object")
    return(NULL)
  }
  
  if(is.null(features))
  {
    stop("updateIndexes: the features object is missing, which is necessary for this function.")
  }
  
  # Model
  if(isModel(obj))
  {
    res <- updateModelIndex(obj = obj, features = features)
  }
  
  # Population
  if(isPopulation(obj))
  {
    res <- list()
    for(i in 1:length(obj))
    {
      res[[i]] <- updateModelIndex(obj = obj[[i]], features = features)
    }
  }
  
  # Model collection
  if(isModelCollection(obj))
  {
    pop <- modelCollectionToPopulation(obj)
    pop.new <- list()
    for(i in 1:length(pop))
    {
      pop.new[[i]] <- updateModelIndex(obj = pop[[i]], features = features)
    }
    res <- listOfModels2ModelCollection(pop.new)
  }
  
  return(res)
}


#' listOfModels2ModelCollection
#'
#' @description Structures a list of predomics objects into a structured collection by k_sparsity. 
#' @param pop: is population (a list) of predomics objects
#' @param nBest: number of elements to return for each sparsity (default:NA)
#' @return an model collection object
#' @export
listOfModels2ModelCollection <- function(pop, nBest = NA)
{
  # this is the old select_nBest_BySparsity while the first listOfModels2ModelCollection is deactivated
  spar <- populationGet_X(element2get = "eval.sparsity", toVec = TRUE, na.rm = TRUE)(pop)
  
  # get for every sparsity, the indices of the individuals that have this sparsity
  real.sparsity <- as.numeric(names(table(spar)))
  real.sparsity <- real.sparsity[real.sparsity != 0] # delete sparsity 0 if any
  
  # get the index of samples with that sparsity
  indivBySparsity <- lapply(real.sparsity, 
                            function(x, spar) 
                            {
                              which(spar == x)
                            }, spar)
  
  res <- lapply(seq_along(indivBySparsity), 
                function(x, indivBySparsity)
                { 
                  pop[indivBySparsity[[x]]] 
                }, 
                indivBySparsity)
  names(res) <- lapply(real.sparsity, function(x) paste("k", x, sep = "_"))
  
  
  # selection
  if(!is.na(nBest))
  {
    res <- lapply(res, function(x) # for each k_sparsity
    {
      x <- sortPopulation(pop = x, evalToOrder = "fit_", decreasing = TRUE)
      
      # filter by number of elements we would like to keep
      if(length(x) > nBest)
        x <- x[1:nBest]
      return(x)
    })
  }
  
  return(res)
}


#' Transform a model collection to a population (or list of model objects)
#' @param mod.collection: a modelCollection object organized by k_sparsity
#' @export
modelCollectionToPopulation <- function(mod.collection)
{
  
  if(!isModelCollection(mod.collection))
  {
    warning("modelCollectionToPopulation: unvalid model collection ... returning NULL")
    return(NULL)
  }
  
  k_sparsity <- names(mod.collection)
  
  res <- list()
  pop.names <- c()
  for(i in 1:length(mod.collection))
  {
    pop <- mod.collection[[i]]
    # if pop exists
    if(!is.null(pop))
    {
      # check it out
      if(!isPopulation(obj = pop))
      {
        stop("modelCollectionToPopulation: unvalid population")
      }
      pop.names <- c(pop.names, rep(k_sparsity[i],length(pop)))
      res[(length(res)+1):(length(res) + length(pop))] <- pop
    }
  }
  names(res) <- pop.names
  
  return(res)
}


#' listOfModelsToDenseCoefMatrix
#'
#' @description For each model in the list of models it will convert to dense format and convert to a data.frame
#' @param clf: the classifier object
#' @param X: the dataset
#' @param y: the class vector
#' @param list.model: a list of model objects
#' @param rm.empty: remove null models in the list if any (default:TRUE)
#' @param order.row: order rows by occurence (default:TRUE)
#' @return an data frame with model coefficients in rows
#' @export
listOfModelsToDenseCoefMatrix <- function(clf, X, y, list.models, rm.empty = TRUE, order.row = TRUE)
{
  if(!isPopulation(list.models))
  {
    stop("listOfModelsToDenseCoefMatrix: please provide a valid population of models")
  }
  check.X_y_w(X = X, y = y)
  
  # Transform the model objects into dense vectors
  pop.dense <- listOfModelsToListOfDenseVec(clf = clf, X=X, y=y, list.models = list.models)
  
  # concatenate them into a matrix
  pop.dense <- as.matrix(do.call(cbind.data.frame, pop.dense))
  rownames(pop.dense) <- rownames(X); colnames(pop.dense) <- paste(clf$learner, c(1:ncol(pop.dense)), sep="_")
  mod.names <- colnames(pop.dense)
  
  if(any(is.na(rownames(pop.dense))))
  {
    print("listOfModelsToDenseCoefMatrix: some features are NA in pop.noz ... omitting them")
    pop.dense <- pop.dense[!is.na(rownames(pop.dense)),]
  }
  
  # if order the rows by participation occurance
  if(order.row)
  { 
    pop.dense.noz <- pop.dense[order(rowSums(pop.dense !=0 , na.rm = TRUE), decreasing = TRUE),]
    if(any(class(pop.dense.noz) == "numeric")) # if a single model
    {
      pop.dense.noz <- as.matrix(as.data.frame(pop.dense.noz))
    }
  } else
  {
    pop.dense.noz <- pop.dense
  }
  
  if(any(is.na(rownames(pop.dense.noz))))
  {
    print("listOfModelsToDenseCoefMatrix: some features are NA in pop.noz ... omitting them")
    pop.dense.noz <- pop.dense.noz[!is.na(rownames(pop.dense.noz)),]
  }
  
  if(rm.empty)
  {
    # filter rows that are zero for all the models
    features.ord <- rownames(pop.dense.noz)
    ind.tokeep <- rowSums(pop.dense.noz != 0, na.rm = TRUE) != 0
    pop.dense.noz <- pop.dense.noz[ind.tokeep,]
    if(any(class(pop.dense.noz) == "numeric")) # if a single model
    {
      pop.dense.noz <- as.matrix(as.data.frame(pop.dense.noz))
      rownames(pop.dense.noz) <- features.ord[ind.tokeep]
      colnames(pop.dense.noz) <- mod.names
    }
  }
  return(pop.dense.noz)
}


# Create a function that transforms a population of model objects onto a dataframe to be plotted
#' populationToDataFrame
#'
#' @description For each model in the list of models it will extract each attribute and create a dataframe needed for further exploration
#' @param pop: a list of model objects, (i.e a population of models)
#' @param attributes: the list of attributes that we wish to have in the data.frame (default:"learner","language","fit_", "unpenalized_fit_", "auc_", "accuracy_", "cor_", "aic_", "intercept_", "eval.sparsity", "sign_","precision_", "recall_","f1_")
#' @return an data frame with attributes for each model
#' @export
populationToDataFrame <- function(pop, attributes = c("learner","language","fit_", "unpenalized_fit_",
                                                      "auc_", "accuracy_", "cor_", "aic_", "intercept_",
                                                      "eval.sparsity", "sign_","precision_", "recall_","f1_"))
{
  
  if(!isPopulation(pop))
  {
    stop("populationToDataFrame: Please provide a valid population object")
  }
  mod <- pop[[1]]
  # attributes <- names(mod)
  # ind.toomit <- match(c("indices_", "names_", "coeffs_", "score_", "confusionMatrix_", 
  #                       "mate", "selected", "toBeMutated"),attributes)
  # attributes <- attributes[-ind.toomit]
  
  ind.match <- match(attributes, names(mod))
  if(any(is.na(ind.match)))
  {
    stop(paste("populationToDataFrame: unknown attributes",attributes[is.na(ind.match)]))
  }
  
  df <- paste("mod",1:length(pop),sep="_")
  ind.null <- c()
  # for each attribute get the values
  for(i in 1:length(attributes))
  {
    x <- populationGet_X(element2get = attributes[i], toVec = TRUE, na.rm = FALSE)(pop)
    #df <- data.frame(df,x)
    if(is.null(x))
    {
      ind.null <- c(ind.null,i)
    }else
    {
      if(attributes[i]=="eval.sparsity")
      {
        x <- factor(as.character(x),levels = as.character(unique(x)))
      }
      df <- data.frame(df,x)
    }
  }
  df <- df[,-1] # get out the first one that was used to start the df
  if(length(ind.null)>1)
  {
    colnames(df) <- attributes[-ind.null]  
  }else
  {
    colnames(df) <- attributes
  }
  rownames(df) <- paste("mod",1:length(pop),sep="_")
  return(df)
}

################################################################
# ADDERS, REMOVERS
################################################################

# Function used to add a model to a model collection
addModelToModelCollection <- function(mod.collection, model)
{
  spar <- length(model$indices_)
  name <- paste("k", spar, sep = "_")
  if(length(mod.collection) > 0)
  {
    mod.collection[[name]][[length(mod.collection[[name]]) + 1]] <- model
  } else 
  {
    mod.collection[[name]] <- list(model)
  }
  return(mod.collection)
}

# Function used to add a list of model to a model collection
addListOfModelsToModelCollection <- function(mod.collection, model.list)
{
  for(i in 1:length(model.list))
  {
    mod.collection <- addModelToModelCollection(mod.collection, model.list[[i]])
  }
  return(mod.collection)
}



################################################################
# COMPUTING PREVALENCE, ABUNDANCE & other POP ATTRIBUTES
################################################################

#' Evaluate the prevalence of a given model
#'
#' @description Evaluate the prevalence of a given model
#' @param mod: a model object
#' @param X: dataset where to compute the prevalence
#' @return A vector containing the prevalence of each feature
#' @export
evaluatePrevalence <- function(mod, X)
{
  # sanity checks todo
  data <- X[mod$indices_,]
  data <- (data > 0)+0.0
  if(length(mod$indices_) > 1)
  {
    prev1 <- apply(data, 1, sum)
  }else
  {
    prev1 <- sum(data)
  }
  return(prev1)
}


#' Evaluates the prevalence of a list of features in the whole dataset and per each class
#'
#' @description Evaluate the prevalence of a given model
#' @param features: a list of features or features indexes for which we wish to compute prevalence
#' @param X: dataset where to compute the prevalence
#' @param y: if provided it will also compute hte prevalence per each class (default:NULL)
#' @param prop: weather to compute the prevalence in number or as a proportion (default:TRUE)
#' @param zero.value: the value that specifies what is zero. This can be a different than 0 in log transformed data for instance (default = 0)
#' @return A list containing the prevalence in the whole dataset as well as classes (if provided)
#' @export
getFeaturePrevalence <- function(features, X, y = NULL, prop = TRUE, zero.value = 0)
{
  
  if(all(is.character(features)))
  {
    # check weather they exist in the dataset's rownames
    ind <- match(features,rownames(X))
    if(any(is.na(ind)))
    {
      stop(paste("getFeaturePrevalence: ",sum(is.na(ind)),"features are not found as rownames of X"))
    }
  }
  
  # if features are indices than check weather they are in the range
  if(all(is.numeric(features)))
  {
    if(max(features > nrow(X)))
    {
      stop("getFeaturePrevalence: feature indexes are out of range. Please check again.")
    }
    else
    {
      ind <- features
    }
  }
  
  # get the corresponding data
  if(ncol(X) == 1)
  {
    data <- as.matrix(X[ind,])
  }else{
    data <- X[ind,]
  }
  
  if(length(features) == 1)
  {
    # when only 1 feature
    data <- t(as.matrix(data))
    rownames(data) <- features
  }
  
  res.all <- list()
  # compute the prevalence for the whole vector
  res <- rowSums(data!= zero.value, na.rm = TRUE) # they can be negative if logged
  if(prop) res <- res / ncol(data)
  names(res) <- features
  res.all[[1]] <- res
  
  # also for each class
  if(!is.null(y))
  {
    lev <- names(table(y))
    for(i in 1:length(lev))
    {
      if(length(features) == 1)
      {
        res <- sum(data[,y==lev[i]] != zero.value, na.rm = TRUE) # they can be negative if logged
      }else
      {
        if(ncol(X) == 1)
        {
          if(table(y)[lev[i]] == 0)
          {
            res <- rep(0, nrow(data))
            names(res) <- features
          }else
          {
            res <- (data[,y==lev[i]] != zero.value)+0.0 # they can be negative if logged  
          }
        }else
        {
          res <- rowSums(data[,y==lev[i]] != zero.value, na.rm = TRUE) # they can be negative if logged
          names(res) <- features
        }
      }
      if(table(y)[lev[i]] != 0)
      {
        if(prop) 
        {
          res <- res / table(y)[lev[i]]  
        }
      }
      res.all[[i+1]] <- res
    }
    names(res.all) <- c("all",lev)
  }else
  {
    names(res.all) <- c("all")
  }
  return(res.all)
}


#' Get the model that has the highest minimal prevalence in its features
#'
#' @description Get the model that has the highest minimal prevalence in its features
#' @param pop: a population of model objects
#' @param X: dataset where to compute the prevalence
#' @param evalToOrder: which score should we use to order the models and select them (default:fit_)
#' @param selected: the number of selected models (default:0). If 0, everything is returned.
#' @return a model or a list of model objects
#' @export
getMaxMinPrevalenceModel <- function(pop, X = NULL, evalToOrder = "fit_", selected = 0)
{
  if(!isPopulation(pop))
  {
    stop("getMaxMinPrevalenceModel: please provide a valid population")
  }
  
  if(is.null(X))
  {
    stop("getMaxMinPrevalenceModel: please provide a valid X object")
  }
  
  # sort the population
  pop             <- sortPopulation(pop, evalToOrder = evalToOrder)
  
  # compute the prevalence for each feature of each model
  prev            <- lapply(pop, evaluatePrevalence, X) 
  # minimum of the prevalence of each feature
  prev.min        <- unlist(lapply(prev, min)) 
  names(prev.min) <- c(1:length(prev.min))
  prev.sort       <- sort(prev.min, decreasing = TRUE) 
  
  # we return everything
  if(selected == 0)
  {
    res           <- pop
    return(res)
  }
  
  # if we need to select one or some
  if(selected == 1)
  {
    res           <- pop[[which.max(prev.min)]] #max of minimum
  }else
  {
    res           <- lapply(as.numeric(names(prev.sort[1:min(selected, length(prev.min))])),function(x) pop[[x]])
  }
  
  return(res)
}


# this function gets the abundance of a dataset and produces an object to be used for plots
# for the whole dataset and by class
getAbundance <- function(data, y=NULL, prop=TRUE){
  res.all <- list()
  # for the whole vector
  res.all[[1]] <- t(data)
  
  # also for each class
  if(!is.null(y)){
    lev <- names(table(y))
    for(i in 1:length(lev)){
      res <- t(data[,y==lev[i]])
      res.all[[i+1]] <- res
    }
    names(res.all) <- c("all",lev)
  }else{
    names(res.all) <- c("all")
  }
  return(res.all)
}


#' computeCardEnrichment
#'
#' @description Computes statistic for enrichment of the cardinality of a score for a two class vector
#' @param v.card.mat: a dataframe with the cardinality of each feature (columns) and each group in the y vector (rows)
#' @param y: the vector containing the class specification for each sample
#' @return a data.frame with the statistics computed
#' @export
computeCardEnrichment <- function(v.card.mat, y)
{
  # build the res object
  res <- list()
  res$v.card.mat <- v.card.mat
  res$y <- y
  
  if(is.null(y))
  {
    return(res)
  }
  
  card.all <- table(y)
  if(length(card.all)!=2) stop("The number of classes should be 2!")
  dat <- v.card.mat[names(card.all),]
  if(ncol(v.card.mat) == 1)
  {
    dat <- as.matrix(dat)
    colnames(dat) <- colnames(v.card.mat)
  }
  # commpute the negative
  dat.negative <- dat
  for(i in 1:length(card.all))
  {
    dat.negative[i,] <- card.all[i] - dat.negative[i,]  
  }
  
  # compute the chisq.test() p-values
  chisq.p <- c()
  chisq.mat.list <- list()
  for(i in 1:ncol(v.card.mat))
  {
    chisq.mat.list[[i]] <- chisq.mat <- rbind(dat[,i],dat.negative[,i]); rownames(chisq.mat) <- c("present","missing")
    chisq.p <- c(chisq.p, signif(chisq.test(chisq.mat)$p.value,3))
    #barplot(chisq.mat,col=c("orangered", "darkolivegreen2"))
  }
  names(chisq.mat.list) <- names(chisq.p) <- colnames(v.card.mat)
  chisq.p[is.nan(chisq.p)] <- 1 # put maximal p-val when warning
  chisq.q <- p.adjust(chisq.p, method = "bonferroni")
  #plot(chisq.q<0.05)
  
  # add statistics to object
  res$card.all <- card.all
  res$chisq.p <- chisq.p
  res$chisq.q <- chisq.q
  return(res)
}

# this function takes a list of vectors with proportions and transforms it in a melted version
# class -1 prevalence is multiplied by -1 for the plot
meltScoreList <- function(v=v.prop, prepare.for.graph=TRUE, topdown=TRUE)
{
  if(!is.list(v)) stop("v should be a list of vectors of the same size.")
  v.prop.melt <- data.frame(matrix(NA, nrow=length(v), ncol=0))
  for(i in 1:length(v))
  {
    v.prop.melt <- data.frame(v.prop.melt,
                              t(data.frame(names(v[[i]]), 
                                           v[[i]], 
                                           rep(names(v)[i], length(v[[i]]))
                              )
                              )
    )
  }
  v.prop.melt <- data.frame(t(v.prop.melt)); 
  colnames(v.prop.melt) <- c("name","score","group")
  # transform score into a number
  v.prop.melt$score <- as.numeric(as.character(v.prop.melt$score))
  # fix factor level order
  if(topdown) 
  {
    v.prop.melt$name <- factor(v.prop.melt$name, levels = rev(names(v$all)))
  }else
  {
    v.prop.melt$name <- factor(v.prop.melt$name, levels = names(v$all))
  }
  # put prevalence score 1 in negative for the plot
  if(prepare.for.graph) v.prop.melt[v.prop.melt$group=="-1",]$score <- v.prop.melt[v.prop.melt$group=="-1",]$score *-1
  
  return(v.prop.melt)
}


#' Plot performance scores for multiple learners.
#'
#' @description summarySE gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95\%).
#' @param data: a data frame
#' @param groupvars: a vector containing names of columns that contain grouping variables
#' @param na.rm: a boolean that indicates whether to ignore NA's
#' @param conf.interval: the percent range of the confidence interval (default is 95\%)
#' @return A transformed data frame with information on the different errors and confidence.
## @import plyr
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) 
{
  # This was taken from 'http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) 
  {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop, .fun = function(xx, col) 
  {
    c(N    = length2(xx[[col]], na.rm=na.rm),
      mean = mean   (xx[[col]], na.rm=na.rm),
      sd   = sd     (xx[[col]], na.rm=na.rm)
    )
  },
  measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


#' Computes different metrics for a given distributions
#'
#' @description This function computes to compute a certain number of metrics on a dataset for each variable 
#' (rows, such as prevalence, quartile distribution, etc.)
#' @param data: a data frame containing the data to be treated.
#' @return a data frame containing different metrics: variance_to_mean, signal_to_noise, variation_coefficient, efficiency and quartile_dispertion
#' @export
computeFeatureMetrics <- function (data) 
{
  # la prevalence
  prev <- rowSums(data>0)/ncol(data)
  # l'information mutuelle
  quart_disp <- apply(data, 1, quartileDispersion)
  moy <- apply(data, 1, mean)
  std <- apply(data, 2, sd)
  var <- apply(data, 2, var)
  q1 <- apply(data, 2, quantile, 0.25)
  q2 <- apply(data, 2, quantile, 0.5)
  q3 <- apply(data, 2, quantile, 0.75)
  variance_to_mean <- (var/moy)
  variance_to_mean[is.nan(variance_to_mean)] <- 0
  signal_to_noise <- (moy/std)
  signal_to_noise[is.nan(signal_to_noise)] <- 0
  variation_coefficient <- (std/moy)
  variation_coefficient[is.nan(variation_coefficient)] <- 0
  efficiency <- (std^2/moy^2)
  efficiency[is.nan(efficiency)] <- 0
  quartile_dispertion <- (q3 - q1)/q2
  quartile_dispertion[is.nan(quartile_dispertion) | is.infinite(quartile_dispertion)] <- 0
  res <- data.frame(variance_to_mean, 
                    signal_to_noise, 
                    variation_coefficient, 
                    efficiency, 
                    quartile_dispertion)
  return(res)
}




#' evaluates the feature importance in a population of models
#'
#' @description This function perturbes the dataset by shuffling one at a time a subset of features that appear in a population of models
#' and recomputes the evaluation of those models. The mean deltas of the score to consider will give a measure of importance. Two methods 
#' are implemented: the first (extensive), will shuffle feature by feature multiple times and will compute the evaluation for the whole 
#' population of models, which can be very time consuming. The second (optimized) and the default approach consists on using a different 
#' seed when shuffling a given feature and computing the population. In this setting it is not needed to run multiple seeds on the whole 
#' dataset. This procedure is designed to be applied in cross validation.
#' @param pop: a population of models to be considered. This population will be filtered if filter.ci = TRUE (default) using the interval 
#' confidence computed around the best model using a binomial distribution.
#' @param X: dataset used to classify
#' @param y: variable to predict
#' @param clf: an object containing the different parameters of the classifier
#' @param score: the attribute of the model to be considered in the evaluation (default:fit_)
#' @param filter.ci: filter the population based on the best model confidence interval (default:TRUE)
#' @param method: Two methods are implemented: the first (extensive), will shuffle feature by feature multiple times and will compute the 
#' evaluation for the whole population of models, which can be very time consuming. The second (optimized) and the default approach consists 
#' on using a different seed when shuffling a given feature and computing the population.
#' @param seed: one or more seeds to be used in the extensive method shuffling (default:c(1:10). For the optimized method only the first seed will be used 
#' and the rest of the seeds that are needed for each model will be incremented from there.
#' @param aggregation: the method to be used to aggregate the evaluation for a the whole population (default: mean), but can be either mean or median.
#' @param verbose: wether to print out information during the execution process.
#' @return a data.frame with features in rows and the population mean/median score for each model*seed of the population
#' @export
evaluateFeatureImportanceInPopulation <- function(pop, X, y, clf, score = "fit_", filter.ci = TRUE, method = "optimized", 
                                                  seed = c(1:10), aggregation = "mean", verbose = TRUE)
{
  if(!isPopulation(pop))
  {
    if(clf$params$warnings) warning("evaluateFeatureImportanceInPopulation: please masure to provide a valid population.")
    return(NULL)
  }
  
  # Sanity checks
  check.X_y_w(X = X, y = y)
  
  if(!method %in% c("optimized","extensive"))
  {
    stop("evaluateFeatureImportanceInPopulation: please provide a valid method either optimized or extensive.")
  }else
  {
    if(method == "optimized") 
    {
      optimized <- TRUE
    }else 
    {
      optimized <- FALSE
    }
  }
  
  if(!aggregation %in% c("mean","median"))
  {
    stop("evaluateFeatureImportanceInPopulation: please provide a valid aggregation method either mean or median.")
  }else
  {
    if(aggregation == "mean")
    {
      aggr.mean <- TRUE
    }else 
    {
      aggr.mean <- FALSE
    }
  }
  
  if(verbose) print(paste("There are",length(pop), "models in this population"))
  if(filter.ci)
  {
    # select the best models in the population after penalizing for sparsity
    pop         <- selectBestPopulation(pop, p = 0.05, k_penalty = 0.75/100)
    if(verbose) print(paste("There are",length(pop), "models in this population after filtering"))
  }
  
  if(!isPopulation(pop))
  {
    if(clf$params$warnings) warning("evaluateFeatureImportanceInPopulation: no models are left after best model selection.")
    return(NULL)
  }
  
  # Reevaluate the population in X (which is the x_test) in generalization
  pop <- evaluatePopulation(X = X, y = y, clf = clf, pop = pop, eval.all = TRUE, force.re.evaluation = TRUE, mode = "test")
  
  # compute the presence of features in the models
  fa            <- makeFeatureAnnot(pop = pop, X = X, y = y, clf = clf)
  feat.prev     <- rowSums(fa$pop.noz != 0, na.rm = TRUE)
  # order by prevalence in the population of models
  feat.prev     <- feat.prev[order(feat.prev, decreasing = TRUE)]
  
  # create a matrix presence mask
  feat.prez.mask     <- fa$pop.noz != 0 & !is.na(fa$pop.noz)
  
  # the pool of features
  feat.pool     <- names(feat.prev)
  if(verbose) print(paste("Feature prevalence is computed. There are ", length(feat.pool), "features to consider."))
  
  eval.orig     <- populationGet_X(element2get = score, toVec = TRUE, na.rm = FALSE)(pop)
  
  if(optimized)
  {
    # for each feature perturb the data to test for its importance
    res.all <- list()
    
    for(i in 1:length(feat.pool))
    {
      if(verbose) cat(paste(feat.pool[i], "\t"))
      
      # the eval shuffle mask
      eval.shuf.mask <- rep(NA, length(pop))
      names(eval.shuf.mask) <- names(pop)
      # the feature mask
      feat.prez.ind <- feat.prez.mask[feat.pool[i],]
      
      # put data in the new population
      pop.new <- list()
      
      # For each model in the population
      for(j in 1:length(pop[feat.prez.ind]))
      {
        # store a copy to change and reset
        X.shuf <- X
        if(seed[1] == 0)
        {
          # for testing purposes if seed is 0 than no shuffling occurs and the DA should be 0
          ind                   <- c(1:ncol(X))
        }else
        {
          # shuffle the feature j
          set.seed(seed[1] + j)
          ind <- sample(x = c(1:ncol(X)), size = ncol(X), replace = FALSE)
        }
        X.shuf[feat.pool[i],] <- X.shuf[feat.pool[i],ind]
        pop.new[[j]] <- evaluateModel(mod = pop[[i]], X = X.shuf, y = y, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE, mode = "test")
        if(verbose) cat("*")
      } # end loop models
      
      # get the evaluation after perturbation
      eval.shuf <- populationGet_X(element2get = score, toVec = TRUE, na.rm = FALSE)(pop.new)
      # compute the delta of the evaluation before and after perturbation and store it. We call this DA (decreased accuracy)
      eval.shuf.mask[feat.prez.ind] <- c(eval.orig[feat.prez.ind] - eval.shuf)
      
      res.all[[i]] <- eval.shuf.mask
      if(verbose) cat(paste("\n"))
      
    } # end loop each feature
    names(res.all) <- feat.pool
    
  }else # extensive
  {
    # for each feature perturb the data to test for its importance
    res.all <- list()
    
    for(i in 1:length(feat.pool))
    {
      if(verbose) cat(paste(feat.pool[i], "\t")) 
      
      # the eval shuffle mask
      eval.shuf.mask <- rep(NA, length(pop))
      names(eval.shuf.mask) <- names(pop)
      # the feature mask
      feat.prez.ind <- feat.prez.mask[feat.pool[i],]
      
      # we can do this multiple times
      res.f <- c()
      for(j in 1:length(seed))
      {
        # store a copy to change and reset
        X.shuf                <- X
        if(seed[j] == 0)
        {
          # for testing purposes if seed is 0 than no shuffling occurs and the DA should be 0
          ind                   <- c(1:ncol(X))
        }else
        {
          # shuffle the feature j
          set.seed(seed[j])
          ind                   <- sample(x = c(1:ncol(X)), size = ncol(X), replace = FALSE)  
        }
        
        X.shuf[feat.pool[i],] <- X.shuf[feat.pool[i],ind]
        pop.eval              <- evaluatePopulation(X = X.shuf, y = y, clf = clf, pop = pop[feat.prez.ind], 
                                                    eval.all = TRUE, force.re.evaluation = TRUE, mode = "test")
        # get the evaluation after perturbation
        eval.shuf             <- populationGet_X(element2get = score, toVec = TRUE, na.rm = FALSE)(pop.eval)
        # compute the delta of the evaluation before and after perturbation and store it. We call this DA (decreased accuracy)
        eval.shuf.mask[feat.prez.ind] <- c(eval.orig[feat.prez.ind] - eval.shuf)
        # concatenate if multiple perturbations
        res.f                 <- c(res.f, eval.shuf.mask)
        if(verbose) cat("*")
      } # end loop seeds
      
      if(verbose) cat(paste("\n")) 
      res.all[[i]] <- res.f
    } # end loop each feature
    names(res.all) <- feat.pool
  }
  
  # put all the DAs in a data frame
  res.all.df                <- t(data.frame(res.all))
  # standard deviation
  SDA                       <- apply(res.all.df, 1, sd, na.rm = TRUE)
  # Prevalence of DA
  PDA                       <- rowSums(!is.na(res.all.df))
  
  # transform the data
  if(aggr.mean)
  {
    # compute the MDA (mean decreased accuracy)
    eval.aggr                 <- apply(res.all.df, 1, mean, na.rm = TRUE)
  }else
  {
    # compute the MDA (median decreased accuracy)
    eval.aggr                 <- apply(res.all.df, 1, median, na.rm = TRUE)
  }
  
  names(eval.aggr)            <- feat.pool
  res                         <- list(feat.catalogue = feat.pool,
                                      feat.catalogue.annot =  fa$feature.df[feat.pool,],
                                      feat.pop.prev = feat.prev,
                                      feat.pop.da.list = res.all,
                                      feat.pop.da.df = res.all.df,
                                      mda = eval.aggr,
                                      sda = SDA,
                                      pda = PDA/ncol(res.all.df) # as a percentage
  )
  
  return(res)
}


# Coefficient of Quartile Deviation
# note this does not work well for sparse data, thus adding an option for non zero
#' @export
quartileDispersion <- function(v, nonzero=TRUE){
  if(nonzero) v <- v[v!=0] 
  # https://en.wikipedia.org/wiki/Quartile_coefficient_of_dispersion
  # http://www.emathzone.com/tutorials/basic-statistics/quartile-deviation-and-its-coefficient.html
  q1 <- quantile(v, 0.25)
  q3 <- quantile(v, 0.75)
  return((q3 - q1)/(q3 + q1))
}


# #' Get the fitting score of a model object
# #'
# #' @description Get the fitting score of a model object.
# #' @param mod : a model object
# #' @return a fitting score
# getFitModel <- function(mod){return(mod$fit_)}
# 
# 
# #' Get the fitting score of a list a models
# #'
# #' @description Get the fitting score of a list a models.
# #' @param pop : a list of models
# #' @return a vector of fitting scores
# getFitModels <- function(pop){
#   res <- rep(NA,length(pop))
#   for(i in 1: length(res)){
#     res[i] <- getFitModel(pop[[i]])
#   }
#   return(res)
# }



#' Evaluates the sign for a given feature this is the old getMgsVsTraitSignDiscr function
#'
#' @description Evaluates the sign for a given feature this is the old getMgsVsTraitSignDiscr function.
#' @import foreach
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the response vector
#' @param clf: the classifier parameter object
#' @param parallel.local: weather or not to run in //
#' @return a vector of +1 & -1 for each variable
#' @export
getSign <- function(X, y, clf = NULL, parallel.local = FALSE)  
{
  
  if(!is.matrix(X)) X <- as.matrix(X)
  
  # check dimensions
  if(ncol(X) != length(y))
  {
    if(nrow(X) == length(y))
    {
      # this happens when only one element is in X it will transpose
      X <- t(X)
    }else
    {
      stop("getSign: wrong dimensions, this should not happen ... stopping.")
    }
  }
  
  if(is.numeric(y)) # if regression (correlation)
  {
    # for optimization, in pearson, convert first in rank
    y <- rank(y)
    # if dimension 1
    if(is.null(dim(X)) & length(X)==length(y))
    {
      X <- rank(X)
      # res <-  as.numeric(sign(cor.test(X,y, method="spearman")$estimate))
      # for optimization reasons the cor will be a pearson implementation, y is already ranked for objective being "cor" performed in the fit()
      # We just need to take the score to obtain the same result as a spearman correlation.
      res <- sign(cor(na.omit(data.frame(X[i,], y)), method = "pearson")[1,2])
    } else
    { 
      # if more than one
      # transform X to rank for spearman correlation (using pearson faster one). This is a change that is active only for this function
      X.rank <- X
      for(i in 1: nrow(X)){ X.rank[i,] <- rank(X[i,]) }
      X <- X.rank # so that we don't change the code hereafter
      
      if(parallel.local) # test if //
      {
        res <- foreach(i = 1:nrow(X))  %dopar%
        { 
          #as.numeric(sign(cor.test(X[i,],y, method="spearman")$estimate))
          # for optimization reasons the cor will be a pearson implementation, y is already ranked for objective being "cor" performed in the fit()
          # We just need to tank the score to obtain the same result as a spearman correlation.
          sign(cor(na.omit(data.frame(X[i,], y)), method = "pearson")[1,2])
        }
        res <- unlist(res); names(res) <- rownames(X)
      } else
      {
        res <- rep(NA, nrow(X))
        names(res) <- rownames(X)
        for (i in 1:nrow(X))
        {
          #res[i]  <- as.numeric(sign(cor.test(X[i,],y, method="spearman")$estimate))
          # for optimization reasons the cor will be a pearson implementation, y is already ranked for objective being "cor" performed in the fit()
          # We just need to tank the score to obtain the same result as a spearman correlation.
          
          res[i]  <- sign(cor(na.omit(data.frame(X[i,], y)), method = "pearson")[1,2])
        }
      }
    }
    
  }else # testing classes (classification)
  {
    # if dimension 1
    if(is.null(dim(X)) & length(X)==length(y))
    {
      elements <- table(as.vector(y))
      cat1 <- names(elements)[1]
      cat1.cl <- -1
      cat2 <- names(elements)[2]
      cat2.cl <- 1
      
      m1 <- mean(X[y == cat1], na.rm = TRUE) 
      m2 <- mean(X[y == cat2], na.rm = TRUE)
      
      if(all(is.na(c(m1,m2))))
      {
        res <- NA
      }else
      {
        if(is.na(m1))
        {
          res <- cat2.cl
        }else
        {
          if (m1 < m2)
          {
            res <- cat2.cl
          } else
          {
            res <- cat1.cl
          }
        }
        
        if(is.na(m2))
        {
          res <- cat1.cl
        }else
        {
          if (m1 < m2)
          {
            res <- cat2.cl
          } else
          {
            res <- cat1.cl
          }
        }
      }
    } else
    { # if more than one
      elements    <- table(as.vector(y))
      if(length(elements) != 2)
      {
        stop("individual: classification setting, only two classes should be provided.")
      }
      cat1 <- names(elements)[1]
      cat1.cl <- -1
      cat2 <- names(elements)[2]
      cat2.cl <- 1
      
      # create a mask for the results
      res         <- rep(NA, nrow(X))
      names(res)  <- rownames(X)
      
      if(parallel.local) # test if //
      {
        res <- foreach(i = 1:nrow(X))  %dopar%
        {
          # These two oparations are quite slow
          m1 <- mean(X[i,y == cat1], na.rm = TRUE) 
          m2 <- mean(X[i,y == cat2], na.rm = TRUE)
          
          if(all(is.na(c(m1,m2))))
          {
            NA
          }else
          {
            if(is.na(m1))
            {
              cat2.cl
            }else
            {
              if (m1 < m2)
              {
                cat2.cl
              } else
              {
                cat1.cl
              }
            }
            
            if(is.na(m2))
            {
              cat1.cl
            }else
            {
              if (m1 < m2)
              {
                cat2.cl
              } else
              {
                cat1.cl
              }
            }
          }
        }
        res <- unlist(res); names(res) <- rownames(X)
      } else
      {
        for (i in 1:nrow(X))
        {
          # These two oparations are quite slow
          m1 <- mean(X[i,y == cat1], na.rm = TRUE) 
          m2 <- mean(X[i,y == cat2], na.rm = TRUE)
          
          if(all(is.na(c(m1,m2))))
          {
            res[i] <- NA
          }else
          {
            if(is.na(m1))
            {
              res[i] <- cat2.cl
            }else
            {
              if (m1 < m2)
              {
                res[i] <- cat2.cl
              } else
              {
                res[i] <- cat1.cl
              }
            }
            
            if(is.na(m2))
            {
              res[i] <- cat1.cl
            }else
            {
              if (m1 < m2)
              {
                res[i] <- cat2.cl
              } else
              {
                res[i] <- cat1.cl
              }
            }
          }
        } # end for
      } # end else //
    } # end else dimension
  } # end else numerical
  
  if(any(is.na(res)))
  {
    cat("... getSign: some signs were not found (NAs)... setting them to default 1\n")
    res[is.na(res)] <- 1
  }
  
  if(any(res == 0))
  {
    cat("... getSign: some signs were not found (NAs)... setting them to default 1\n")
    res[res == 0] <- 1
  }
  
  return(res)
}


################################################################
# CROSS VALIDATION FUNCTIONS
################################################################

# This function is of the R package "caret" . Not mine.
# https://github.com/topepo/caret/blob/master/pkg/caret/R/createFolds.R
# Splits the data into k groups to be used with k cross-validation
# update: from Edi (2017 nov 13)
#' @export
create.folds <- function (y, k = 10, list = TRUE, returnTrain = FALSE, seed = NULL)
{
  if (k < 0) {
    AssertTwoLevels <- function(fold, y) {
      if (length(fold) != 1) {
        while (length(table(y[fold])) == 1 | length(table(y[-fold])) ==
               1) {
          fold <- sample(c(1:length(y)), length(fold))
        }
      }
      return(fold)
    }
    out <- lapply(1:length(y), function(x) AssertTwoLevels(sample(c(1:length(y)),
                                                                  -k), y))
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
                        sep = "")
    return(out)
  }
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2)
      cuts <- 2
    if (cuts > 5)
      cuts <- 5
    y <- cut(y, unique(quantile(y, probs = seq(0, 1, length = cuts))),
             include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      seqVector <- rep(1:k, numInClass[i]%/%k)
      if (is.null(seed)) {
        if (numInClass[i]%%k > 0)
          seqVector <- c(seqVector, sample(1:k, numInClass[i]%%k))
      }
      else {
        set.seed(seed)
        if (numInClass[i]%%k > 0)
          seqVector <- c(seqVector, sample(1:k, numInClass[i]%%k))
      }
      if (is.null(seed)) {
        foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
      }
      else {
        set.seed(seed)
        foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
      }
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
                        sep = "")
    if (returnTrain)
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}

################################################################
# FILTERING PROCEDURES
################################################################

#' Selects a the top k features that are significantly associated with the class to predict
#' 
#' @description Runs statistics on the data and selects a subset of k features that are the most significant. 
#' An accelerated version is implemented based on the BioQC package for the Mann-Whittney tests. Besides filtering 
#' this function can be used in a more larger statistical context.
#' @import BioQC
#' @param data: the dataset X
#' @param trait: is the equivalent of y (class, or numerical)
#' @param k: the number of features (default:10)
#' @param type: the statistics to be run (default:wilcoxon)
#' @param restrict: Run the statistics in a subset of the dataset (default: a vector of all TRUE)
#' @param multiple.adjust: the multiple testing adjustment method (default:BH)
#' @param paired: wether paired statistics should be run (default:FALSE)
#' @param sort: return variables sorted by p-value significance (default:TRUE)
#' @param verbose: print out information indicating progress (default:FALSE)
#' @param verbose.step: Showing a 1 percent progress.
#' @param return.data: if (default:FALSE) this returns the statistics of X, otherwise the restricted data subset
#' @param accelarate: use a turbo method developped by bioQC (default:FALSE). There is an issue when executing in batch.
#' @export
filterfeaturesK <- function(data, 
                            trait, 
                            k = 10, 
                            type = "wilcoxon", 
                            restrict = rep(TRUE, ncol(data)),  
                            multiple.adjust = "BH", 
                            paired = FALSE, 
                            sort = TRUE,
                            verbose = FALSE,
                            verbose.step = NULL,
                            return.data = FALSE, 
                            accelerate = FALSE
) 
{
  if(!type %in% c("spearman","pearson","wilcoxon","t.test"))
  {
    stop("filterfeaturesK: unknown type! Please provide one of the following: spearman, pearson, wilcoxon, t.test")
  }
  
  # sanity checks
  if(length(trait) != ncol(data))
  {
    stop("filterfeaturesK: incompatible dimensions between data and trait")
  }
  
  # if trait is not a vector get it
  if(!is.null(nrow(trait)))
  {
    if(nrow(trait)==1)
    {
      cl <- class(trait[1,1])
      if(cl=="numeric")
      {
        trait <- as.numeric(trait)
      }
      if(cl=="character")
      {
        trait <- as.character(trait)
      }
    }else
    {
      stop("filterfeaturesK: trait should be a vector")
    }
  }
  
  # fix to correlation if not specified
  if(any(class(trait) == "numeric"))
  {
    if(!(type == "spearman" | type == "pearson"))
    {
      if(length(table(trait)) == 2)
      {
        trait <- as.factor(trait)
        trait.val <- names(table(trait))
      }else
      {
        if(verbose) cat("... ... trait seems to be a class with two levels, changing type to Spearman\n")
        type <- "spearman" 
        trait <- rank(trait)
        trait.val <- NA
      }
    }
  }else
  {
    # if this is a classification problem
    trait.val <- names(table(trait))
  }
  
  # convert data to matrix
  if(!is.matrix(data))
  {
    data <- as.matrix(data)
  }
  
  
  # create reults mask
  res <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = 5))
  rownames(res) <- rownames(data)
  colnames(res) <- c("rho", "rho2", "p", "q", "status")
  
  # Handle switches
  mean.test <- TRUE
  cor.spearman <- FALSE
  mean.test.t <- TRUE
  
  # correlation
  if(type == "spearman" | type == "pearson") 
  {
    mean.test <- FALSE
    if(type == "spearman")
    {
      cor.spearman <- TRUE
      trait <- rank(trait) # for speedup
    }
  }
  
  # classification problem
  if(mean.test)
  {
    if(type == "wilcoxon")
    {
      mean.test.t <- FALSE
    }
    
    if (length(table(trait)) != 2) 
    {
      stop("filterfeaturesK: You can't use this test. The trait should contain only two categories!")
    }
  }
  
  # printing out logs in verbose mode
  if(mean.test)
  {
    if(mean.test.t) 
    {
      if(verbose) cat("... ... filterfeaturesK: Executing T test\n")
    }else
    {
      if(verbose) cat("... ... filterfeaturesK: Executing Wilcoxon test\n")
    }
  }else{
    if(cor.spearman) 
    {
      if(paired)
      {
        if (verbose) cat("... ... filterfeaturesK: Spearman correlation mode between two classes\n")  
      }else
      {
        if (verbose) cat("... ... filterfeaturesK: Spearman correlation\n")
      }
    }else
    {
      if(paired)
      {
        if (verbose) cat("... ... filterfeaturesK: Pearson correlation mode between two classes\n")
      }else
      {
        if (verbose) cat("... ... filterfeaturesK: Pearson correlation\n")
      }
    }
  }
  
  validity <- rep(FALSE, nrow(data))
  
  if(verbose) # start progress bar 
  {
    comp <- ""
    cat("|") 
  }
  
  if(is.null(verbose.step)) verbose.step <- round(nrow(data)/100)
  
  if(verbose.step == 0) verbose.step <- 1
  
  # for each variable
  for (i in 1:nrow(data)) 
  {
    # slim completion progress indicator 100 points (%)
    if (verbose & i %% verbose.step == 0) 
    {
      cat(paste(comp, "*", sep = ""))
    }
    # or more detail
    #if (verbose & i %% verbose.step == 0) print(paste("... ... ... ",round(i/nrow(data)*100)," % - (",i," features)", sep=""))
    
    vt <- trait[restrict] # trait restricted
    vd <- data[i, restrict] # data restricted
    valid <- FALSE # default is not valid
    
    # MEAN Tests
    if(mean.test)
    {
      # check validity
      val.tab <- table(vt, vd)
      if(all(rowSums(val.tab) > 2) & 
         nrow(val.tab) == 2 & 
         ncol(val.tab) > 1 # when there are only zeros for instance
      )
      {
        valid <- TRUE
        validity[i] <- TRUE
      }
      
      if(valid)
      {
        # T-TEST
        if(mean.test.t) 
        {
          # apply test
          try(tmp <- stats::t.test(vd ~ vt, paired = paired), silent = TRUE)
          try(res[i, "p"] <- tmp$p.value, silent = TRUE) # store results
          # determine the status
          if (mean(vd[trait == trait.val[1]], na.rm = TRUE) > mean(vd[trait == trait.val[2]], na.rm = TRUE)) 
          {
            res[i, "status"] <- trait.val[1]
          }
          else {
            res[i, "status"] <- trait.val[2]
          }
        }
        else { # WILCOXON
          # we will use here an optimization that greatly speeds things several times from the BioQC package
          # https://accio.github.io/BioQC/bioqc-efficiency.html
          # we will compute this outside the for loop
          
          if (!accelerate)
          {
            # apply test
            try(tmp <- stats::wilcox.test(vd ~ vt, paired = paired), silent = TRUE)
            try(res[i, "p"] <- tmp$p.value, silent = TRUE)
          }
          
          # determine the status
          if (mean(vd[trait == trait.val[1]], na.rm = TRUE) > mean(vd[trait == trait.val[2]], na.rm = TRUE)) 
          {
            res[i, "status"] <- trait.val[1]
          }
          else 
          {
            res[i, "status"] <- trait.val[2]
          }
        } # end else
      } # end validity
      
    }else # CORRELATION
    {
      # PAIRED
      if (paired) 
      {
        # validity check
        if (length(table(trait)) != 2 | (table(trait)[1] != table(trait)[2])) 
        {
          stop("filterfeaturesK: trait does not seem to be a 2-level categorical variable or with the same prevalence")
        }
        cl1 <- names(table(vt)[1])
        cl1.ind <- (vt == cl1)
        cl2 <- names(table(vt)[2])
        cl2.ind <- (vt == cl2)
        
        # SPEARMAN
        if(cor.spearman) 
        {
          # run test and store results
          try(tmp <- stats::cor.test(rank(vd[cl1.ind]), rank(vd[cl2.ind]), method = "pearson"), silent = TRUE)
          try(res[i, 1] <- tmp$estimate, silent = TRUE)
          if (!is.na(res[i, 1])) 
          {
            if (res[i, 1] > 0) 
            {
              res[i, 5] <- "POS"
            }
            else {
              res[i, 5] <- "NEG"
            }
          } # end if
          res[i, 2] <- res[i, 1]^2
          res[i, 3] <- tmp$p.value
          
          if(!is.na(res[i, 2])) 
          {
            valid <- TRUE
          }
        }else # PEARSON
        {
          # run test and store results
          try(tmp <- stats::cor.test(vd[cl1.ind], vd[cl2.ind], method = "pearson"), silent = TRUE)
          try(res[i, 1] <- tmp$estimate, silent = TRUE)
          if (!is.na(res[i, 1])) 
          {
            if (res[i, 1] > 0) 
            {
              res[i, 5] <- "POS"
            }
            else {
              res[i, 5] <- "NEG"
            }
          } # end if
          res[i, 2] <- res[i, 1]^2
          res[i, 3] <- tmp$p.value
          
          if(!is.na(res[i, 2])) 
          {
            valid <- TRUE
          }
        } # end spearman/pearson
        
      }else # UNPAIRED
      {
        # SPEARMAN
        if(cor.spearman) 
        {
          try(tmp <- stats::cor.test(rank(vd), vt, method = "pearson"), silent = TRUE)
          try(res[i, 1] <- tmp$estimate, silent = TRUE)
          if (!is.na(res[i, 1])) 
          {
            if (res[i, 1] > 0) 
            {
              res[i, 5] <- "POS"
            }
            else {
              res[i, 5] <- "NEG"
            }
          } # end if
          res[i, 2] <- res[i, 1]^2
          res[i, 3] <- tmp$p.value
          
          if(!is.na(res[i, 2])) 
          {
            valid <- TRUE
          }
          
        }else # PEARSON
        {
          try(tmp <- stats::cor.test(vd, vt, method = "pearson"), silent = TRUE)
          try(res[i, 1] <- tmp$estimate, silent = TRUE)
          if (!is.na(res[i, 1])) 
          {
            if (res[i, 1] > 0) 
            {
              res[i, 5] <- "POS"
            }
            else {
              res[i, 5] <- "NEG"
            }
          } # end if
          res[i, 2] <- res[i, 1]^2
          res[i, 3] <- tmp$p.value
          
          if(!is.na(res[i, 2])) 
          {
            valid <- TRUE
          }
        } # end spearman/pearson
      } # end paired/unpaired
      
      if(!valid)
      {
        res[i, "p"] <- NA
        res[i, "status"] <- NA
      } # end validity
      
    } # end mean test / correlation
  } # end for
  
  if(verbose) cat("|\n") # end progress bar
  
  if(accelerate)
  {
    # OPTIMIZATION particular case of the optimized wilcoxon
    if(mean.test & type == "wilcoxon")
    {
      res[, 3] <- wmwTest(x = t(data), indexList = (trait==names(table(trait))[1]), valType="p.two.sided", simplify = TRUE)
      res[!validity, 3] <- NA # to make sure we don't keep non valid data that we won't select as feature selection
    }
  }
  
  # multiple adjustment
  res[, 4] <- stats::p.adjust(res[, "p"], method = multiple.adjust)
  
  # if the result is sorted 
  if(sort)
  {
    # reorder features
    res <- res[order(res$p),]
  }
  
  if(return.data)
  {
    features <- rownames(res)
    res <- as.data.frame(data)[features,]
  }
  
  # send a subset of k rows if k is smaller than the number of features
  if(k < nrow(res))
  {
    res <- res[1:k,]
  } 
  
  return(res)
}

# library(microbenchmark)
# microbenchmark(filterfeaturesK(data = X[1:100,], trait = y, type = "wilcoxon"),
#                filterfeaturesK(data = X[1:100,], trait = y, type = "t.test"),
#                filterfeaturesK(data = X[1:100,], trait = c(1:ncol(X)), type = "spearman"),
#                filterfeaturesK(data = X[1:100,], trait = c(1:ncol(X)), type = "pearson"),
#                times = 10)
# Unit: milliseconds
#                                                                          expr       min        lq      mean    median        uq        max neval cld
# filterfeaturesK(data = X[1:100, ], trait = y, type = "wilcoxon")               75.77135  77.41353 333.64460  77.76377  78.61290 2632.23538    10   a
# filterfeaturesK(data = X[1:100, ], trait = y, type = "t.test")                430.92199 437.45751 455.39362 446.81926 462.45985  531.81801    10   a
# filterfeaturesK(data = X[1:100, ], trait = c(1:ncol(X)), type = "spearman")    45.30576  47.42271  49.26134  47.53064  50.66467   57.84277    10   a
# filterfeaturesK(data = X[1:100, ], trait = c(1:ncol(X)), type = "pearson")     37.95488  40.29546 298.54324  42.15932  50.14797 2579.22371    10   a


#' Selects the most prevalent features in the dataset baset on the provided thresholds.
#' 
#' @description Filters out all features that display a prevalence below a given threshold provided as a number of 
#' observations or percentage. This for the total dataset or by class.
#' @param X: the dataset X
#' @param y: the class vector (default:NULL)
#' @param nb.prevalence: the minimum number of non zero observations (default: 10)
#' @param perc.prevalence: the percentage of non zero observations (default: NULL)
#' @param by.class: wether the filter should be applied by class (default: TRUE)
#' @return the filtered dataset, without the features that do not pass the filter.
#' @export
filterFeaturesByPrevalence <- function(X, y = NULL, nb.prevalence = NULL, perc.prevalence = NULL, by.class = TRUE)
{
  
  if(is.null(y))
  {
    # this is a case when we consider y to be only one class
    y <- as.factor(rep("all",ncol(X)))
    
  }else
  {
    if(any(class(y) == "numeric"))
    {
      y <- as.factor(as.character(y))
    }
    
    if(any(class(y) == "character"))
    {
      y <- as.factor(y)
    }
  }
  
  # sanity check
  check.X_y_w(X, y)
  # Get the prevalence for all the features
  prev <- getFeaturePrevalence(features = rownames(X), X = X, y = y, prop = FALSE)
  prev <- prev[!duplicated(prev)]
  
  card                            <- length(y)
  cards                           <- table(y)
  cards.min                       <- cards;
  cards.min[1:length(cards.min)]  <- rep(NA,length(cards.min))
  card.min                        <- NA
  
  if(is.null(nb.prevalence))
  {
    if(is.null(perc.prevalence))
    {
      stop("filterFeaturesByPrevalence: please provide a value for nb.prevalence or perc.prevalence!")
    }else
    {
      if(perc.prevalence < 0 | perc.prevalence > 100)
      {
        stop("filterFeaturesByPrevalence: please provide a perc.prevalence value between 0 and 100!")
      }
      
      if(by.class)
      {
        if(length(cards.min) > 1)
        {
          # if we have more than one class
          for(i in 1:length(cards))
          {
            cards.min[i] <- round(perc.prevalence * cards[i] /100)
          }
        }else{
          cards.min <- round(perc.prevalence * card /100)  
        }
      }else 
      {
        card.min <- round(perc.prevalence * card /100)
      }
    }
  }else # if nb prevalence is given
  {
    if(!is.null(perc.prevalence))
    {
      stop("filterFeaturesByPrevalence: please chose either nb.prevalence or perc.prevalence!")
    }
    
    if(by.class)
    {
      if(length(cards.min) > 1)
      {
        # if we have more than one class
        for(i in 1:length(cards))
        {
          cards.min[i] <- min(nb.prevalence, cards[i]) # to set up an upper limit
        }
      }else
      {
        cards.min <- nb.prevalence
      }
    }else 
    {
      card.min <- nb.prevalence
    }
  }
  
  if(by.class)
  {
    # create a mask
    ind <- rep(TRUE,nrow(X))
    for(i in 1:length(cards))
    {
      # adjust the mask by taking into account all the filters for each level
      ind <- ind & prev[[names(cards)[i]]] >= cards.min[i]
    }
    tokeep <- names(which(ind))
  }else
  {
    tokeep <- names(which(prev$all >= card.min))
  }
  res <- X[tokeep,]
  
  return(res)
}



#' filterNoSignal: Omits the variables with no information
#' @description This function will clean a dataset from the variables that have no or little information.
#' @param X: the dataset to clean
#' @param side: side=1 means that variables are in the rows. Other than 1 it will transpose the dataset
#' @param threshold: auto, will compute the first derivate of the median(sd)/x and will find an automatic threshold. When threshold is a numerical it will be used as a threshold and when is something else, will automatically be 0.
#' @param verbose: print out information when TRUE (default:FALSE)
#' @export
filterNoSignal <- function(X, side = 1, threshold = "auto", verbose = FALSE)
{
  
  # Side = 1 this is for the rows, if something else we will transpose
  if(side == 1)
  {
    res <- X
  }else
  {
    res <- t(X)
  }
  
  # compute the standard deviation
  s <- apply(res, 1, stats::sd, na.rm = TRUE)
  s.median <- stats::median(s, na.rm = TRUE)
  
  if(verbose) 
  {
    print(paste("The standard deviation distribution is:"))
    print(signif(summary(s)),2)
  }
  
  if(threshold == "auto")
  {
    # Compute distribution of median(s)/x
    th <- seq(1:1000)
    x <- c()
    for(i in 1:1000)
    {
      x <- rbind(x, table(s < s.median / th[i]))
    }
    # ind <- (c(x[,2],0)-c(0,x[,2]))[1001] *-1 # first derivate
    ind <- (c(x[,2],0,0) + c(0,0,x[,2]) - 2*c(0,x[,2],0)) # second derivate
    d2 <- max(ind[-c(1,2)])
    thresh <- stats::median(s) / d2
    if(verbose) print(paste("The best automatic threshold is", thresh))
  }else if(is.numeric(threshold))
  {
    thresh <- threshold
  }else
  {
    # only the empty variables
    thresh <- 0
  }
  
  res <- res[which(s > thresh),]
  if(side!=1) # retranspose
  {
    res <- data.frame(t(res))
  }
  return(res)
}


#' cleanPopulation
#'
#' @description Looks for invalid predomics objects in a population and removes them. 
#' @param pop: is population (a list) of predomics objects
#' @param clf: the classifier object
#' @return a population of predomics objects
#' @export
cleanPopulation <- function(pop, clf)
{
  if(is.null(clf))
  {
    stop("cleanPopulation: please provide a valid clf object")
  }
  
  if(is.null(pop)) # if not null
  {
    if(clf$params$warnings) warning("cleanPopulation: empty population.")  
    return(NULL)
  }
  
  # Test which models are real and which are not and delete them otherwise
  tokeep <- unlist(lapply(pop, isModel))
  pop <- pop[tokeep]
  
  # don't continue if nothing left
  if(length(pop)==0)
  {
    return(NULL)
  }
  
  # Destroy any model if it contains
  # indices == 0
  indices.tags <- unlist(lapply(populationGet_X(element2get = "indices_", toVec = FALSE, na.rm = FALSE)(pop), length))
  if(any(names(table(indices.tags))=="0"))
  {
    pop <- pop[-which(indices.tags==0)]
  }
  
  # coefficients == 0
  coeffs.tags <- unlist(lapply(populationGet_X(element2get = "coeffs_", toVec = FALSE, na.rm = FALSE)(pop), length))
  if(any(names(table(coeffs.tags))=="0"))
  {
    pop <- pop[-which(coeffs.tags==0)]
  }
  
  # coefficients == 0
  names.tags <- unlist(lapply(populationGet_X(element2get = "names_", toVec = FALSE, na.rm = FALSE)(pop), length))
  if(any(names(table(names.tags))=="0"))
  {
    pop <- pop[-which(names.tags==0)]
  }
  
  # recompute everything to make sure we have the same numbers
  indices.tags <- unlist(lapply(populationGet_X(element2get = "indices_", toVec = FALSE, na.rm = FALSE)(pop), length))
  names.tags <- unlist(lapply(populationGet_X(element2get = "names_", toVec = FALSE, na.rm = FALSE)(pop), length))
  coeffs.tags <- unlist(lapply(populationGet_X(element2get = "coeffs_", toVec = FALSE, na.rm = FALSE)(pop), length))
  sparsity.tags <- populationGet_X(element2get = "eval.sparsity", toVec = TRUE, na.rm = FALSE)(pop)
  
  # if there is a difference in number of elements among indices, names, coeffs and sparsity
  checkEmntNb <- function(v){
    if(length(unique(v)==1))
    {
      return(TRUE)
    }else
    {
      return(FALSE)
    }
  }
  
  # make a data frame to test all
  df <- data.frame(indices.tags, 
                   names.tags, 
                   coeffs.tags, 
                   sparsity.tags)
  
  df.test <- apply(df,1,checkEmntNb)
  if(any(!df.test))
  {
    pop <- pop[-which(!df.test)]
  }
  
  return(pop)
}


################################################################
# DIGESTING RESULTS
################################################################

#' Summarize the results from an experiment object
#'
#' @description Sumarizes the results of an experiment object of the type 
#' `obj$classifier` and `obj$crossval`. This is different from the digestMC(),
#' which sumarizes a model collection obj$models
#' @import ggplot2
#' @param obj: The experiment object resulting from the learning process `fit()`
#' @param penalty: A coefficient between 0 and 1, which is applied to penalize 
#' the performance of models as a consequence of model-size. We use this to select
#' the best model of the population of models (default:NULL)
#' @param best.cv: Should we chose the best model based on information learnerd 
#' cross validation (default:TRUE). This will work if the crossvalidation data is 
#' available. If not the best model will be selected with empirical results.
#' @param best.k: If we do not wish to let the algorithm select the model size, 
#' we can fix this by setting the best.k with an integer indicating the number of 
#' variables in the model (default:NULL).
#' @param plot: Should the digested results be plotted ? (default:FALSE)
#' @param omit.na: Omit data with empty results (default:TRUE) 
#' @return an object with digested information such as the best models for each 
#' model-size, their respective scores, the best model.
#' @export
digest <- function(obj, 
                   penalty = NULL, 
                   best.cv = TRUE, 
                   best.k = NULL, 
                   plot = FALSE,
                   omit.na = TRUE)
{
  if(!isExperiment(obj))
  { 
    stop("digest: The object did not pass the sanity check for an Experiment object!") 
  }
  
  # so that it works with below
  if(is.null(penalty)) 
  {
    penalty <- 0 
  }
  
  if(length(penalty) > 1)
  {
    stop("digest: Please provide a single value for penalty!") 
  }
  
  res                       <- list()
  res$learner               <- obj$classifier$learner
  
  # Empirical
  res$best                  <- list()
  res$best$models           <- getNBestModels(obj = obj, 
                                              significance = TRUE, 
                                              by.k.sparsity = TRUE,
                                              k.penalty = penalty,
                                              n.best = 1,
                                              single.best = FALSE,
                                              single.best.cv = best.cv,
                                              single.best.k = best.k,
                                              max.min.prevalence = FALSE,
                                              X = NULL,
                                              verbose = FALSE, 
                                              evalToOrder = "fit_",
                                              return.population = TRUE, # population
                                              unique.control = TRUE
  )
  
  # sanity check
  if(is.null(res$best$models))
  {
    warning("digest: no models are found. Returning empty handed.")
    return(NULL)
  }
  res$best$model            <- getNBestModels(obj = obj,
                                              significance = TRUE,
                                              by.k.sparsity = TRUE,
                                              k.penalty = penalty,
                                              n.best = 5,
                                              single.best = TRUE, # give best
                                              single.best.cv = best.cv, # based on CV
                                              single.best.k = best.k,
                                              max.min.prevalence = FALSE,
                                              X = NULL,
                                              verbose = FALSE,
                                              evalToOrder = "fit_",
                                              return.population = FALSE # population
  )
  
  res$best$scores           <- list()
  res$best$scores$fit_      <- populationGet_X(element2get = "fit_", toVec = TRUE, na.rm = FALSE)(pop = res$best$models)
  res$best$scores$auc_      <- populationGet_X(element2get = "auc_", toVec = TRUE, na.rm = FALSE)(pop = res$best$models)
  res$best$scores$accuracy_ <- populationGet_X(element2get = "accuracy_", toVec = TRUE, na.rm = FALSE)(pop = res$best$models)
  res$best$scores$precision_<- populationGet_X(element2get = "precision_", toVec = TRUE, na.rm = FALSE)(pop = res$best$models)
  res$best$scores$recall_   <- populationGet_X(element2get = "recall_", toVec = TRUE, na.rm = FALSE)(pop = res$best$models)
  res$best$scores$f1_       <- populationGet_X(element2get = "f1_", toVec = TRUE, na.rm = FALSE)(pop = res$best$models)
  res$best$scores$cor_      <- populationGet_X(element2get = "cor_", toVec = TRUE, na.rm = FALSE)(pop = res$best$models)
  res$best$scores$ser_      <- populationGet_X(element2get = "ser_", toVec = TRUE, na.rm = FALSE)(pop = res$best$models)
  res$best$scores$rsq_      <- populationGet_X(element2get = "rsq_", toVec = TRUE, na.rm = FALSE)(pop = res$best$models)
  
  
  # Cross Validation (if activated)
  res$cv <- list()
  if(!is.null(obj$crossVal))
  {
    res$cv                  <-   obj$crossVal
    crossval                <- TRUE
  }else
  {
    warning(paste("crossval information unavailable for", obj$learner))
    crossval                <- FALSE
  }
  
  maj.class <- NA # default value for maj.class
  
  # Majoritary class for all folds
  if(!is.null(obj$classifier$lfolds))
  {
    if(obj$classifier$params$objective == "auc")
    {
      if(obj$classifier$params$evalToFit == "accuracy_")
      {
        maj.class           <- rep(NA, length(obj$classifier$lfolds))
        # if accuracy
        for(i in 1:length(maj.class))
        {
          maj.class[i]      <- max(table(obj$classifier$data$y[-obj$classifier$lfolds[[i]]])/length(obj$classifier$data$y[-obj$classifier$lfolds[[i]]])) 
        }
      }else
      {
        # if AUC
        maj.class           <- 0.5
      }
    } # if correlation it is NA
  }
  
  ##########################################################
  # Plotting the results
  ##########################################################
  if(plot)
  {
    # set the limits
    ylim <- c(0,1)
    
    # make an empty plot in case it does not work
    g.empty <- ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(ylim) + 
      theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
      ylab("") +
      xlab("Model parismony sparse") +
      ggtitle("") +
      geom_hline(aes(yintercept=1), lty=2, col="lightgray") + 
      geom_hline(aes(yintercept=0.5), lty=2, col="lightgray")
    
    # if correlation
    if(obj$classifier$params$objective=="cor")
    {
      if(crossval)
      {
        #-----------------------------------------------------
        # CORRELATION
        #-----------------------------------------------------
        # overall training accuracy learner results
        dat <- res$cv$scores$empirical.cor
        if(omit.na)
        {
          dat <- dat[rowSums(!is.na(dat))!=0,]  
        }
        df <- data.frame(parsimony = rownames(dat),dat)
        df.melt <- melt(df, id.vars = "parsimony")
        df.melt$parsimony <- as.factor(df.melt$parsimony)
        df.melt$parsimony <- factor(df.melt$parsimony,
                                    levels = levels(df.melt$parsimony)[order(as.numeric(gsub("k_","",levels(df.melt$parsimony))))]
        )
        g.cor.emp.cv <- ggplot(data = df.melt, aes(y = value, x = parsimony)) +
          geom_point(aes(color = variable), position=position_jitterdodge(dodge.width=0.9), size = 1, alpha = 0.5) +
          geom_boxplot(notch = FALSE, outlier.colour = NA, position = position_dodge(width=0.9), alpha = 0.3) +
          ylab("cor_") +
          xlab("Model parsimony") +
          ggtitle("Training performance (CV)") +
          ylim(ylim) +
          theme_bw() +
          #geom_hline(yintercept = unique(maj.class), col = "gray", linetype = "dashed") +
          theme(legend.position="bottom", legend.direction="horizontal") +
          guides(colour = "none")
        
        
        # overall testing accuracy learner results
        dat <- res$cv$scores$generalization.cor
        if(omit.na)
        {
          dat <- dat[rowSums(!is.na(dat))!=0,]  
        }
        df <- data.frame(parsimony = rownames(dat),dat)
        df.melt <- melt(df, id.vars = "parsimony")
        df.melt$parsimony <- as.factor(df.melt$parsimony)
        df.melt$parsimony <- factor(df.melt$parsimony,
                                    levels = levels(df.melt$parsimony)[order(as.numeric(gsub("k_","",levels(df.melt$parsimony))))]
        )
        g.cor.gen.cv <- ggplot(data = df.melt, aes(y = value, x = parsimony)) +
          geom_point(aes(color = variable), position=position_jitterdodge(dodge.width=0.9), size = 1, alpha = 0.5) +
          geom_boxplot(notch = FALSE, outlier.colour = NA, position = position_dodge(width=0.9), alpha = 0.3) +
          ylab("cor_") +
          xlab("Model parsimony") +
          ggtitle("Testing performance (CV)") +
          ylim(ylim) +
          theme_bw() +
          #geom_hline(yintercept = unique(maj.class), col = "gray", linetype = "dashed") +
          theme(legend.position="bottom", legend.direction="horizontal") +
          guides(colour = "none")

      }else
      {
        g.cor.emp.cv <- g.empty
        g.cor.gen.cv <- g.empty
      }
      
      # RHO Empirical
      v <- res$best$scores$cor_
      df <- data.frame(value = v, parsimony = names(v))
      df$parsimony <- factor(df$parsimony,
                             levels = levels(df$parsimony)[order(as.numeric(gsub("k_","",levels(df$parsimony))))]
      )
      g.cor.emp <- ggplot(data = df, aes(x = parsimony, y = value, group = 1)) +
        geom_line(aes(color = "gray")) +
        geom_point(size = 2, alpha = 1) +
        scale_color_manual(values = "gray") +
        ylab("cor_") +
        xlab("Model parsimony") +
        ggtitle("Emprical performance") +
        labs(subtitle = paste("L:",obj$classifier$learner,"|F:",signif(res$best$model$cor_,4),"|k:",length(res$best$model$indices_), sep="")) +
        ylim(ylim) +
        theme_bw() +
        geom_hline(yintercept = mean(maj.class), col = "gray", linetype = "dashed") +
        theme(legend.position="bottom", legend.direction="horizontal") +
        guides(colour = "none")
      
      # RSQ Empirical (R squared)
      v <- res$best$scores$rsq_
      df <- data.frame(value = v, parsimony = names(v))
      df$parsimony <- factor(df$parsimony,
                             levels = levels(df$parsimony)[order(as.numeric(gsub("k_","",levels(df$parsimony))))]
      )
      g.rsq.emp <- ggplot(data = df, aes(x = parsimony, y = value, group = 1)) +
        geom_line(aes(color = "gray")) +
        geom_point(size = 2, alpha = 1) +
        scale_color_manual(values = "gray") +
        ylab("rsq_") +
        xlab("Model parsimony") +
        ggtitle("Emprical performance") +
        labs(subtitle = paste("L:",obj$classifier$learner,"|F:",signif(res$best$model$rsq_,4),"|k:",length(res$best$model$indices_), sep="")) +
        ylim(ylim) +
        theme_bw() +
        geom_hline(yintercept = mean(maj.class), col = "gray", linetype = "dashed") +
        theme(legend.position="bottom", legend.direction="horizontal") +
        guides(colour = "none")
      
      # SER Empirical (Standar error of the regression)
      v <- res$best$scores$ser_
      df <- data.frame(value = v, parsimony = names(v))
      df$parsimony <- factor(df$parsimony,
                             levels = levels(df$parsimony)[order(as.numeric(gsub("k_","",levels(df$parsimony))))]
      )
      g.ser.emp <- ggplot(data = df, aes(x = parsimony, y = value, group = 1)) +
        geom_line(aes(color = "gray")) +
        geom_point(size = 2, alpha = 1) +
        scale_color_manual(values = "gray") +
        ylab("ser_") +
        xlab("Model parsimony") +
        ggtitle("Emprical performance") +
        labs(subtitle = paste("L:",obj$classifier$learner,"|F:",signif(res$best$model$ser_,4),"|k:",length(res$best$model$indices_), sep="")) +
        #ylim(ylim) +
        theme_bw() +
        geom_hline(yintercept = mean(maj.class), col = "gray", linetype = "dashed") +
        theme(legend.position="bottom", legend.direction="horizontal") +
        guides(colour = "none")
      
      
      grid.arrange(g.cor.emp.cv, g.cor.gen.cv, 
                   # empricial
                   g.cor.emp, g.rsq.emp, 
                   g.ser.emp,
                   ncol = 2)

    }else # if classification
    {
      if(crossval)
      {
        # make an empty plot in case it does not work
        g.empty <- ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(ylim) + 
          theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
          ylab("") +
          xlab("Model parismony sparse") +
          ggtitle("") +
          geom_hline(aes(yintercept=1), lty=2, col="lightgray") + 
          geom_hline(aes(yintercept=0.5), lty=2, col="lightgray")
        
        #-----------------------------------------------------
        # ACCURACY
        #-----------------------------------------------------
        # overall training accuracy learner results
        dat <- res$cv$scores$empirical.acc
        if(omit.na)
        {
          dat <- dat[rowSums(!is.na(dat))!=0,]  
        }
        df <- data.frame(parsimony = rownames(dat),dat)
        df.melt <- melt(df, id.vars = "parsimony")
        df.melt$parsimony <- as.factor(df.melt$parsimony)
        df.melt$parsimony <- factor(df.melt$parsimony,
                                    levels = levels(df.melt$parsimony)[order(as.numeric(gsub("k_","",levels(df.melt$parsimony))))]
        )
        g.accuracy.emp.cv <- ggplot(data = df.melt, aes(y = value, x = parsimony)) +
            geom_point(aes(color = variable), position=position_jitterdodge(dodge.width=0.9), size = 1, alpha = 0.5) +
            geom_boxplot(notch = FALSE, outlier.colour = NA, position = position_dodge(width=0.9), alpha = 0.3) +
            ylab("accuracy_") +
            xlab("Model parsimony") +
            ggtitle("Training performance (CV)") +
            ylim(ylim) +
            theme_bw() +
            geom_hline(yintercept = unique(maj.class), col = "gray", linetype = "dashed") +
            theme(legend.position="bottom", legend.direction="horizontal") +
            guides(colour = "none")

        # overall testing accuracy learner results
        dat <- res$cv$scores$generalization.acc
        if(omit.na)
        {
          dat <- dat[rowSums(!is.na(dat))!=0,]  
        }
        df <- data.frame(parsimony = rownames(dat),dat)
        df.melt <- melt(df, id.vars = "parsimony")
        df.melt$parsimony <- as.factor(df.melt$parsimony)
        df.melt$parsimony <- factor(df.melt$parsimony,
                                    levels = levels(df.melt$parsimony)[order(as.numeric(gsub("k_","",levels(df.melt$parsimony))))]
        )
        g.accuracy.gen.cv <- ggplot(data = df.melt, aes(y = value, x = parsimony)) +
            geom_point(aes(color = variable), position=position_jitterdodge(dodge.width=0.9), size = 1, alpha = 0.5) +
            geom_boxplot(notch = FALSE, outlier.colour = NA, position = position_dodge(width=0.9), alpha = 0.3) +
            ylab("accuracy_") +
            xlab("Model parsimony") +
            ggtitle("Testing performance (CV)") +
            ylim(ylim) +
            theme_bw() +
            geom_hline(yintercept = unique(maj.class), col = "gray", linetype = "dashed") +
            theme(legend.position="bottom", legend.direction="horizontal") +
            guides(colour = "none")

        #-----------------------------------------------------
        # AUC
        #-----------------------------------------------------
        # overall training accuracy learner results
        dat <- res$cv$scores$empirical.auc
        if(omit.na)
        {
          dat <- dat[rowSums(!is.na(dat))!=0,]  
        }
        df <- data.frame(parsimony = rownames(dat),dat)
        df.melt <- melt(df, id.vars = "parsimony")
        df.melt$parsimony <- as.factor(df.melt$parsimony)
        df.melt$parsimony <- factor(df.melt$parsimony,
                                    levels = levels(df.melt$parsimony)[order(as.numeric(gsub("k_","",levels(df.melt$parsimony))))]
        )
        g.auc.emp.cv <- ggplot(data = df.melt, aes(y = value, x = parsimony)) +
            geom_point(aes(color = variable), position=position_jitterdodge(dodge.width=0.9), size = 1, alpha = 0.5) +
            geom_boxplot(notch = FALSE, outlier.colour = NA, position = position_dodge(width=0.9), alpha = 0.3) +
            ylab("auc_") +
            xlab("Model parsimony") +
            ggtitle("Training performance (CV)") +
            ylim(ylim) +
            theme_bw() +
            geom_hline(yintercept = 0.5, col = "gray", linetype = "dashed") +
            theme(legend.position="bottom", legend.direction="horizontal") +
            guides(colour = "none")
        
        # overall testing accuracy learner results
        dat <- res$cv$scores$generalization.auc
        if(omit.na)
        {
          dat <- dat[rowSums(!is.na(dat))!=0,]  
        }
        df <- data.frame(parsimony = rownames(dat),dat)
        df.melt <- melt(df, id.vars = "parsimony")
        df.melt$parsimony <- as.factor(df.melt$parsimony)
        df.melt$parsimony <- factor(df.melt$parsimony,
                                    levels = levels(df.melt$parsimony)[order(as.numeric(gsub("k_","",levels(df.melt$parsimony))))]
        )
        g.auc.gen.cv <- ggplot(data = df.melt, aes(y = value, x = parsimony)) +
            geom_point(aes(color = variable), position=position_jitterdodge(dodge.width=0.9), size = 1, alpha = 0.5) +
            geom_boxplot(notch = FALSE, outlier.colour = NA, position = position_dodge(width=0.9), alpha = 0.3) +
            ylab("auc_") +
            xlab("Model parsimony") +
            ggtitle("Testing performance (CV)") +
            ylim(ylim) +
            theme_bw() +
            geom_hline(yintercept = 0.5, col = "gray", linetype = "dashed") +
            theme(legend.position="bottom", legend.direction="horizontal") +
            guides(colour = "none")
      }else
      {
        g.accuracy.emp.cv <- g.empty
        g.accuracy.gen.cv <- g.empty
        g.auc.emp.cv <- g.empty
        g.auc.gen.cv <- g.empty
      }
      
      if(all(is.na(unlist(res$best$scores))))
      {
        g.accuracy.emp <- g.empty
        g.auc.emp <- g.empty
        g.recall.emp <- g.empty
        g.precision.emp <- g.empty
      }else
      {
        # ACCURACY Empirical
        v <- res$best$scores$accuracy_
        df <- data.frame(value = v, parsimony = names(v))
        df$parsimony <- as.factor(df$parsimony)
        df$parsimony <- factor(df$parsimony,
                                    levels = levels(df$parsimony)[order(as.numeric(gsub("k_","",levels(df$parsimony))))]
        )
        g.accuracy.emp <- ggplot(data = df, aes(x = parsimony, y = value, group = 1)) +
          geom_line(aes(color = "gray")) +
          geom_point(size = 2, alpha = 1) +
          scale_color_manual(values = "gray") +
          ylab("accuracy_") +
          xlab("Model parsimony") +
          ggtitle("Emprical performance") +
          labs(subtitle = paste("L:",obj$classifier$learner,"|F:",signif(res$best$model$accuracy_,4),"|k:",length(res$best$model$indices_), sep="")) +
          ylim(ylim) +
          theme_bw() +
          geom_hline(yintercept = mean(maj.class), col = "gray", linetype = "dashed") +
          theme(legend.position="bottom", legend.direction="horizontal") +
          guides(colour = "none")
          
        
        # AUC Empirical
        v <- res$best$scores$auc_
        df <- data.frame(value = v, parsimony = names(v))
        df$parsimony <- as.factor(df$parsimony)
        df$parsimony <- factor(df$parsimony,
                               levels = levels(df$parsimony)[order(as.numeric(gsub("k_","",levels(df$parsimony))))]
        )
        g.auc.emp <- ggplot(data = df, aes(x = parsimony, y = value, group = 1)) +
          geom_line(aes(color = "gray")) +
          geom_point(size = 2, alpha = 1) +
          scale_color_manual(values = "gray") +
          ylab("auc_") +
          xlab("Model parsimony") +
          ggtitle("Emprical performance") +
          labs(subtitle = paste("L:",obj$classifier$learner,"|F:",signif(res$best$model$auc_,4),"|k:",length(res$best$model$indices_), sep="")) +
          ylim(ylim) +
          theme_bw() +
          geom_hline(yintercept = mean(maj.class), col = "gray", linetype = "dashed") +
          theme(legend.position="bottom", legend.direction="horizontal") +
          guides(colour = "none")
        
        # RECALL Empirical
        v <- res$best$scores$recall_
        df <- data.frame(value = v, parsimony = names(v))
        df$parsimony <- as.factor(df$parsimony)
        df$parsimony <- factor(df$parsimony,
                               levels = levels(df$parsimony)[order(as.numeric(gsub("k_","",levels(df$parsimony))))]
        )
        g.recall.emp <- ggplot(data = df, aes(x = parsimony, y = value, group = 1)) +
          geom_line(aes(color = "gray")) +
          geom_point(size = 2, alpha = 1) +
          scale_color_manual(values = "gray") +
          ylab("recall_") +
          xlab("Model parsimony") +
          ggtitle("Emprical performance") +
          labs(subtitle = paste("L:",obj$classifier$learner,"|F:",signif(res$best$model$recall_,4),"|k:",length(res$best$model$indices_), sep="")) +
          ylim(ylim) +
          theme_bw() +
          geom_hline(yintercept = 0.5, col = "gray", linetype = "dashed") +
          theme(legend.position="bottom", legend.direction="horizontal") +
          guides(colour = "none")
        
        # PRECISION Empirical
        v <- res$best$scores$precision_
        df <- data.frame(value = v, parsimony = names(v))
        df$parsimony <- as.factor(df$parsimony)
        df$parsimony <- factor(df$parsimony,
                               levels = levels(df$parsimony)[order(as.numeric(gsub("k_","",levels(df$parsimony))))]
        )
        g.precision.emp <- ggplot(data = df, aes(x = parsimony, y = value, group = 1)) +
          geom_line(aes(color = "gray")) +
          geom_point(size = 2, alpha = 1) +
          scale_color_manual(values = "gray") +
          ylab("precision_") +
          xlab("Model parsimony") +
          ggtitle("Emprical performance") +
          labs(subtitle = paste("L:",obj$classifier$learner,"|F:",signif(res$best$model$precision_,4),"|k:",length(res$best$model$indices_), sep="")) +
          ylim(ylim) +
          theme_bw() +
          geom_hline(yintercept = 0.5, col = "gray", linetype = "dashed") +
          theme(legend.position="bottom", legend.direction="horizontal") +
          guides(colour = "none")
        
      } # end existing scores
      
      grid.arrange(g.accuracy.emp.cv, g.accuracy.gen.cv, 
                   g.auc.emp.cv, g.auc.gen.cv, 
                   # empricial
                   g.accuracy.emp, g.auc.emp, 
                   g.recall.emp, g.precision.emp, 
                   ncol = 2)
      
    } # end classification
  }
  
  # return digested results
  return(res)
}


#' Summarize the results from a given modelCollection object
#'
#' @title digestModelCollection
#' @description Sumarizes the results of a model collection object of the type clf_res$models. This is different from the digest() which sumarizes an experiment clf_res$classifier$models
#' @param obj: a modelCollection object
#' @param X: the dataset (default = NULL)
#' @param clf: the classifier object
#' @param k.penalty: the penalty to apply for sparsity (default:0).
#' @param mmprev: activate the max.min.prevalence selector (default:FALSE)
#' @return an object with sumarized results such as the best models for each k_sparse, their respective scores, the best model, etc
#' @export
digestModelCollection <- function(obj, X = NULL, clf, k.penalty = 0, mmprev = FALSE)
{
  if(!isModelCollection(obj))
  {
    # STOP for debug
    stop("digestModelCollection: The object is not a modelCollection")
  }
  
  res <- list()
  res$learner <- clf$learner
  
  # extract the best models
  # TODO include getTheBestModel & getBestIndividual in getBestModels
  if(mmprev & !is.null(X))
  {
    res$best_models           <- getNBestModels(obj = obj, 
                                                significance = TRUE, 
                                                by.k.sparsity = TRUE,
                                                k.penalty = k.penalty,
                                                n.best = 1,
                                                single.best = FALSE,
                                                single.best.cv = FALSE,
                                                single.best.k = NULL,
                                                max.min.prevalence = TRUE, # max min prevalence
                                                X = NULL,
                                                verbose = FALSE, 
                                                evalToOrder = "fit_",
                                                return.population = TRUE, # population
                                                unique.control = TRUE
    )
    res$best_model            <- getNBestModels(obj = obj,
                                                significance = TRUE,
                                                by.k.sparsity = TRUE,
                                                k.penalty = k.penalty,
                                                n.best = 1,
                                                single.best = TRUE, # !!! give best
                                                single.best.cv = FALSE, # based on CV
                                                single.best.k = NULL,
                                                max.min.prevalence = TRUE, # max min prevalence
                                                X = NULL,
                                                verbose = FALSE,
                                                evalToOrder = "fit_",
                                                return.population = FALSE # this will be a model
    )
  }else
  {
    res$best_models           <- getNBestModels(obj = obj, 
                                                significance = TRUE, 
                                                by.k.sparsity = TRUE,
                                                k.penalty = k.penalty,
                                                n.best = 1,
                                                single.best = FALSE,
                                                single.best.cv = FALSE,
                                                single.best.k = NULL,
                                                max.min.prevalence = FALSE,
                                                X = NULL,
                                                verbose = FALSE, 
                                                evalToOrder = "fit_",
                                                return.population = TRUE, # population
                                                unique.control = TRUE
    )
    res$best_model            <- getNBestModels(obj = obj,
                                                significance = TRUE,
                                                by.k.sparsity = TRUE,
                                                k.penalty = k.penalty,
                                                n.best = 1,
                                                single.best = TRUE, # !!! give best
                                                single.best.cv = FALSE, # based on CV
                                                single.best.k = NULL,
                                                max.min.prevalence = FALSE,
                                                X = NULL,
                                                verbose = FALSE,
                                                evalToOrder = "fit_",
                                                return.population = FALSE # this will be a model
    )
  }
  res$best_models.fit       <- populationGet_X(element2get = "fit_", toVec = TRUE, na.rm = FALSE)(pop = res$best_models)
  res$best_models.accuracy  <- populationGet_X(element2get = "accuracy_", toVec = TRUE, na.rm = FALSE)(pop = res$best_models)
  res$best_models.precision <- populationGet_X(element2get = "precision_", toVec = TRUE, na.rm = FALSE)(pop = res$best_models)
  res$best_models.recall    <- populationGet_X(element2get = "recall_", toVec = TRUE, na.rm = FALSE)(pop = res$best_models)
  res$best_models.f1_       <- populationGet_X(element2get = "f1_", toVec = TRUE, na.rm = FALSE)(pop = res$best_models)
  res$best_models.cor_      <- populationGet_X(element2get = "cor_", toVec = TRUE, na.rm = FALSE)(pop = res$best_models)
  
  return(res)
}


#' Merge a list of empirical scores form digest results
#' 
#' @title mergeMeltScoreEmpirical
#' @description mergeMeltScoreEmpirical returns a data frames that contain the performance of each digest in the list with their sparsity.
#' @param list.results.digest: a list of digest objects one for each learner used. For example, list(res.terda.digest, res.terga.digest, res.terbeam.digest) 
#' @param k_catalogue: the k_catalogue that will serve to build the result matrix
#' @param score: which score is to be used for value (default: fit_)
#' @return a data.frame
#' @export
mergeMeltScoreEmpirical <- function(list.results.digest, k_catalogue, score = "fit_")
{
  if(is.null(list.results.digest))
  {
    stop("mergeMeltScoreEmpirical: list results is empty!")
  }
  
  if(is.null(list.results.digest[[1]]$best$scores[[score]]))
  {
    warning(paste("mergeMeltScoreEmpirical: results from this score are not found!",score))
    return(NULL)
    
  }else
  {
    tmp.value <- c(); tmp.method <- c()
    for(i in 1:length(list.results.digest))
    {
      values <- list.results.digest[[i]]$best$scores[[score]][k_catalogue]; 
      method <- rep(names(list.results.digest)[i], length(k_catalogue))
      if(!is.null(values))
      {
        names(values) <- k_catalogue  
        tmp.value <- c(tmp.value, values)
        tmp.method <- c(tmp.method, method)
      }
    }
    res <- data.frame(tmp.method, names(tmp.value), tmp.value)
    colnames(res) <- c("method","variable","value"); rm(tmp.method, tmp.value, values)
    res$variable  <- factor(res$variable, levels = k_catalogue)
    res$method    <- factor(res$method, levels = names(list.results.digest))
  }
  
  return(res)
}


#' Merge a list of cross validation scores form digest results
#' 
#' @import reshape2
#' @title mergeMeltScoreCV
#' @description mergeMeltScoreCV returns a list of data frames that contain the performance of each digest in the list with their sparsity.
#' @param list.results.digest: a list of digest objects one for each learner used. For example, list(res.terda.digest, res.terga.digest, res.terbeam.digest) 
#' @param k_catalogue: the k_catalogue that will serve to build the result matrix
#' @param generalization: get results from CV generalization (if TRUE) or empirical otherwise (default: TRUE)
#' @param score: the name of the score that needs to be used for the whole dataset visualization.
#' @return a list of two data.frames
#' @export
mergeMeltScoreCV <- function(list.results.digest, k_catalogue, generalization = TRUE, score = "auc_")
{
  
  fuseCvResults <- function(list.results.digest, inds.folds, ind, mask)
  {
    tmp <- mask # empty matrix
    indr <- intersect(rownames(mask), rownames(list.results.digest[[1]]$cv$scores[[ind]]))
    indc <- paste("fold", length(inds.folds[[1]]), sep="_")
    tmp[indr, ] <- list.results.digest[[1]]$cv$scores[[ind]][indr, inds.folds[[1]]] # initialize
    cv.res <- tmp
    
    if (length(list.results.digest)>=2)
    {
      for(i in 2:length(list.results.digest))
      { 
        tmp <- mask # empty matrix
        indr <- intersect(rownames(mask), rownames(list.results.digest[[i]]$cv$scores[[ind]]))
        indc <- paste("fold", length(inds.folds[[1]]), sep="_")
        tmp[indr, ] <- list.results.digest[[i]]$cv$scores[[ind]][indr, inds.folds[[i]]] # initialize
        cv.res <- data.frame(cv.res, tmp) 
      }
    }
    
    cv.res <- t(cv.res)
    return(cv.res)
  }
  
  e.nb      <- length(list.results.digest) # number of classifiers
  e.name    <- names(list.results.digest)
  
  # get the number of k-folds for each experiment, which can be differnt from one to another.
  nbfolds <- rep(NA, e.nb)
  names(nbfolds) <- e.name
  for(i in 1:e.nb)
  {
    if(!is.null(ncol(list.results.digest[[i]]$cv$scores$empirical.auc)))
    {
      nbfolds[i] <- ncol(list.results.digest[[i]]$cv$scores$empirical.auc)  
    }else
    {
      nbfolds[i] <- NA
    }
  }
  
  # get the same indices for the folds (minimal number)
  inds.folds <- NULL
  if(min(nbfolds, na.rm = TRUE) != max(nbfolds, na.rm = TRUE))
  {
    inds.folds <- list()
    for(i in 1:e.nb)
    {
      if(!is.na(nbfolds[i]))
      {
        set.seed(1234)
        inds.folds[[i]] <- sample(x = seq_len(nbfolds[i]), size = min(nbfolds), replace = FALSE)  
      }else
      {
        inds.folds[[i]] <- NULL
      }
    }
  }else
  { # else just get the increasing index
    inds.folds <- list()
    for(i in 1:e.nb)
    {
      if(!is.na(nbfolds[i]))
      {
        inds.folds[[i]] <- seq_len(nbfolds[i])
      }else
      {
        inds.folds[[i]] <- NULL
      }
    }
  }
  
  if(!is.null(list.results.digest[[1]]$best$scores[[score]]))
  {
    # Whole datasest empirical
    tmp.value <- c()
    tmp.method <- c()
    all.folds <- c()
    
    for(i in 1:e.nb) # for all the classifiers
    {
      values <- list.results.digest[[i]]$best$scores[[score]][k_catalogue]; 
      method <- rep(e.name[i], length(k_catalogue))
      if(!is.null(values))
      {
        names(values) <- k_catalogue  
        tmp.value <- c(tmp.value, values)
        tmp.method <- c(tmp.method, method)
        
        # focus on the inds.folds index
        all.folds <- c(all.folds, colnames(list.results.digest[[i]]$cv$scores$empirical.auc[,inds.folds[[i]]]))
      }
    }
    # put results
    res <- data.frame(tmp.method, names(tmp.value), tmp.value)
    colnames(res) <- c("method","variable","value"); rm(tmp.method, tmp.value, values)
    res$variable  <- factor(res$variable, levels = k_catalogue)
    res$method    <- factor(res$method, levels = e.name)
  }
  
  mask            <- data.frame(matrix(nrow = length(k_catalogue), ncol = min(nbfolds, na.rm = TRUE)))
  rownames(mask)  <- k_catalogue
  colnames(mask)  <- paste("fold",seq_len(ncol(mask)),sep="_")
  
  # CROSS VALIDATION 
  if(score=="auc_")
  {
    # auc
    if(generalization) # generalization results
    {
      cv.res <- fuseCvResults(list.results.digest, inds.folds, ind="generalization.auc", mask)
    }else
    {
      cv.res <- fuseCvResults(list.results.digest, inds.folds, ind="empirical.auc", mask)
    }
  }else if(score=="accuracy_")
  {
    # accuracy
    if(generalization) # generalization results
    {
      cv.res <- fuseCvResults(list.results.digest, inds.folds, ind="generalization.acc", mask)
    }else
    {
      cv.res <- fuseCvResults(list.results.digest, inds.folds, ind="empirical.acc", mask)
    }
  }else if(score=="recall_")
  {
    # recall
    if(generalization) # generalization results
    {
      cv.res <- fuseCvResults(list.results.digest, inds.folds, ind="generalization.rec", mask)
    }else
    {
      cv.res <- fuseCvResults(list.results.digest, inds.folds, ind="empirical.rec", mask)
    }
  }else if(score=="precision_")
  {
    # precision
    if(generalization) # generalization results
    {
      cv.res <- fuseCvResults(list.results.digest, inds.folds, ind="generalization.prc", mask)
    }else
    {
      cv.res <- fuseCvResults(list.results.digest, inds.folds, ind="empirical.prc", mask)
    }
  }else if(score=="f1_")
  {
    # f1-score
    if(generalization) # generalization results
    {
      cv.res <- fuseCvResults(list.results.digest, inds.folds, ind="generalization.f1s", mask)
    }else
    {
      cv.res <- fuseCvResults(list.results.digest, inds.folds, ind="empirical.f1s", mask)
    }
  }else if(score=="cor_")
  {
    # regression
    if(generalization) # generalization results
    {
      cv.res <- fuseCvResults(list.results.digest, inds.folds, ind=6, mask)
    }else
    {
      cv.res <- fuseCvResults(list.results.digest, inds.folds, ind=5, mask)
    }
  }else
  {
    warning(paste("mergeMeltScoreCV: The score:", score, "is not a valid one."))
  }
  
  if(all(is.na(cv.res)))
  {
    warning(paste("mergeMeltScoreCV: The score:", score, "does not contain any data."))
    return(NULL)
  }
  # add a column method
  method <- c()
  for(i in 1:e.nb) # for all the classifiers
  { 
    method <- c(method, rep(e.name[i], min(nbfolds, na.rm = TRUE))) 
  }
  
  # melt things together
  if(score=="auc_" | score=="accuracy_" | score=="recall_" | score=="precision_" | score=="f1_" | score=="cor_")
  {
    # auc
    cv.res                      <- data.frame(cv.res, method)
    cv.res$method               <- factor(cv.res$method, levels = e.name) # reorder the factor
    cv.res.long                 <- melt(cv.res, id = c('method')) #from package reshape2
    cv.res.long$method          <- factor(cv.res.long$method, levels = e.name) # reorder the factor
    # compute for each method and variable the standard error
    cv.res.long.summary         <- summarySE(data = cv.res.long, measurevar="value", groupvars=c("method", "variable"), na.rm = TRUE)
    cv.res.long.summary$method  <- factor(cv.res.long.summary$method, levels = e.name) # reorder the factor
  }else
  {
    warning(paste("mergeMeltScoreCV: The score:",score,"is not a valid one."))
  }
  
  res       <- list()
  res$cv    <- cv.res.long
  res$cv.se <- cv.res.long.summary
  
  return(res)
}



#' Merge a list of cross validation scores form digest results
#' 
#' @import reshape2
#' @title mergeMeltBestScoreCV
#' @description mergeMeltBestScoreCV returns a list of data frames that contain the best performance of the different learners without any focus on sparsity.
#' @param list.results.digest: a list of digest objects one for each learner used. For example, list(res.terda.digest, res.terga.digest, res.terbeam.digest) 
#' @param k_catalogue: the k_catalogue that will serve to build the result matrix (default:NULL)
#' @param score: the name of the score that needs to be used for the whole dataset visualization.
#' @param min.kfold.nb: wether we should restrict all experiments in the smallest number of k-folds of a comparative analyses (default = FALSE)
#' @return a list of two data.frames
#' @export
mergeMeltBestScoreCV <- function(list.results.digest, k_catalogue = NULL, score="auc_", penalty = 0, min.kfold.nb = FALSE)
{
  
  # if the best.k is set
  if(!is.null(penalty))
  {
    if(length(penalty) != length(list.results.digest) & length(penalty) != 1)
    {
      stop("mergeMeltBestScoreCV: Please provide a vector of penalty of the same length as there are experiments. 
           When the penalty is to be disabled for a given experiment, please provide 0 for that one!")
    }
    if(any(!is.numeric(penalty)))
    {
      stop("mergeMeltBestScoreCV: Please provide a list penalty that are numeric!")
    }
  }else
  {
    penalty <- rep(0, length(list.results.digest)) 
    names(penalty) <- names(list.results.digest)
  }
  
  if(score=="auc_" | score=="accuracy_" | score=="recall_" | score=="precision_" | score=="f1_" | score=="cor_")
  {
    e.nb     <- length(list.results.digest) # number of classifiers
    e.name <- names(list.results.digest) # names of the classifiers
    
    # get the number of k-folds for each experiment, which can be differnt from one to another.
    nbfolds <- rep(NA, e.nb)
    names(nbfolds) <- e.name
    for(i in 1:e.nb)
    {
      if(!is.null(ncol(list.results.digest[[i]]$cv$scores$empirical.auc)))
      {
        nbfolds[i] <- ncol(list.results.digest[[i]]$cv$scores$empirical.auc)  
      }else
      {
        nbfolds[i] <- NA
      }
    }
    
    inds.folds <- NULL
    if(min(nbfolds, na.rm = TRUE) != max(nbfolds, na.rm = TRUE) & min.kfold.nb)
    {
      inds.folds <- list()
      for(i in 1:e.nb)
      {
        if(!is.na(nbfolds[i]))
        {
          set.seed(1234)
          inds.folds[[i]] <- sample(x = seq_len(nbfolds[i]), size = min(nbfolds), replace = FALSE)  
        }else
        {
          inds.folds[[i]] <- NULL
        }
      }
    }else
    { # else just get the increasing index
      inds.folds <- list()
      for(i in 1:e.nb)
      {
        if(!is.na(nbfolds[i]))
        {
          inds.folds[[i]] <- seq_len(nbfolds[i])
        }else
        {
          inds.folds[[i]] <- NULL
        }
      }
    }
    
    # Sanity check
    if(length(penalty) == 1)
    {
      penalty <- rep(penalty, e.nb)      
    }
    
    # Abreviations for the results
    key <- data.frame(real = c("auc_","accuracy_","recall_","precision_","f1_","cor_"),
                      abbrv = c("auc","acc","rec","prc","f1s","cor"))
    rownames(key) <- key$real
    
    # build a list with results on each classifier.
    res.split <- list() 
    if(is.null(k_catalogue))
    {
      k_catalogue <- rownames(list.results.digest[[1]]$cv$k$acc)
    }
    
    for(i in 1:e.nb) # for each classifier
    {
      emp.name <- paste("empirical", as.character(key[score,]$abbrv), sep=".")
      gen.name <- paste("generalization", as.character(key[score,]$abbrv), sep=".")
      # for each classifier
      emp.data <- list.results.digest[[i]]$cv$scores[[emp.name]][k_catalogue, inds.folds[[i]]]
      gen.data <- list.results.digest[[i]]$cv$scores[[gen.name]][k_catalogue, inds.folds[[i]]]
      # plot for debug
      # par(mfrow=c(2,1)); image(as.matrix(t(emp.data))); image(as.matrix(t(gen.data)))
      
      if(!is.null(emp.data) & !is.null(gen.data))
      {
        # compute the penalty data
        spar      <- as.numeric(gsub("k_","",k_catalogue))
        mean_score <- rowMeans(emp.data, na.rm = TRUE)
        mean_score_penalty  <- mean_score - penalty[i] * spar
        
        if(isModel(obj = list.results.digest[[i]]$best$model))
        {
          k_sparse  <- paste0("k_",list.results.digest[[i]]$best$model$eval.sparsity)  
          #plot(mean_score, mean_score_penalty, ylim=c(0,1),xlim=c(0.5,1))
          ind <- match(k_sparse, k_catalogue)
          best_empirical <- rep(NA,ncol(emp.data))
          best_generalization <- rep(NA,ncol(gen.data))
          best_empirical <- as.numeric(emp.data[ind,])
          best_generalization <- as.numeric(gen.data[ind,])
        }else
        {
          warning(paste("mergeMeltBestScoreCV: No best model was identified for this learner based on the constrains decided."))
          # plot(mean_score, mean_score_penalty, ylim=c(0,1),xlim=c(0.5,1))
          # ind <- match(k_sparse, k_catalogue)
          best_empirical <- rep(NA,ncol(emp.data))
          best_generalization <- rep(NA,ncol(gen.data))
          # best_empirical <- as.numeric(emp.data[ind,])
          # best_generalization <- as.numeric(gen.data[ind,])
        }
        
        # for GLMNET models
        if(all(isModelSotaGLMNET(list.results.digest[[i]]$best$model)))
        {
          ind.fold <- apply(!is.na(emp.data),2, function(x) {return(max(spar[x]))})
          ind.fold <- match(ind.fold, spar)
          
          for(k in 1:ncol(emp.data))
          {
            best_empirical[k] <- as.numeric(emp.data[ind.fold[k],k])
            best_generalization[k] <- as.numeric(gen.data[ind.fold[k],k])
          }
          # This is the number of features selected in the best model whole dataset. 
          # It will be different form the generalization due to the fact that experiments are selected in a subset of sparsity
          k_sparse <- paste0("k_",list.results.digest[[i]]$best$model$eval.sparsity)
        }
        
        # for TERDA models
        if(all(isModelTerda(list.results.digest[[i]]$best$model)))
        {
          # ind.fold <- apply(!is.na(emp.data),2, function(x) {return(max(spar[x]))})
          # ind.fold <- match(ind.fold, spar)
          # 
          # for(k in 1:ncol(emp.data))
          # {
          #   best_empirical[k] <- as.numeric(emp.data[ind.fold[k],k])
          #   best_generalization[k] <- as.numeric(gen.data[ind.fold[k],k])
          # }
          
          # This is the number of features selected in the best model whole dataset. 
          # It will be different form the generalization due to the fact that experiments are selected in a subset of sparsity
          # k_sparse <- paste0("k_",list.results.digest[[i]]$best$model$eval.sparsity)
          
          # find the closest in sparsity
          df <- as.data.frame(cbind(sparproxy = abs(spar - list.results.digest[[i]]$best$model$eval.sparsity),
                                    nbnona = rowSums(!is.na(emp.data))))
          k_sparse <- k_catalogue[which.min(df$sparproxy)]
          #k_sparse <- k_catalogue[which.min(abs(spar - list.results.digest[[i]]$best$model$eval.sparsity))]
          ind <- match(k_sparse, k_catalogue)
        }
        res.split[[i]] <- data.frame(best_empirical, best_generalization, k_sparse)
      }
      else
      {
        res.split[[i]] <- data.frame(best_empirical = NA, best_generalization = NA)[-1,]
      }
    }
    names(res.split) <- e.name
    
    if(all(unlist(lapply(res.split, nrow))==0))
    {
      warning(paste("mergeMeltBestScoreCV: The score:",score,"does not contain any data."))
      return(NULL)
    }
    
    # melt the list together
    #melt(res.split, level="method")
    cv.res.melt         <- melt(res.split, level="method", measure.vars = c("best_empirical","best_generalization"))
    cv.res.melt$Lmethod <- factor(cv.res.melt$Lmethod, levels = e.name) # reorder the factor
    # empirical
    cv.res.emp          <- cv.res.melt[cv.res.melt$variable=="best_empirical",]
    cv.res.emp.summary  <- summarySE(data=cv.res.emp, measurevar="value", groupvars=c("Lmethod"), na.rm = TRUE)
    # reorder the levels
    ind.order <- match(names(res.split), cv.res.emp.summary$Lmethod)
    cv.res.emp.summary <- cv.res.emp.summary[ind.order,]
    # generalization
    cv.res.gen          <- cv.res.melt[cv.res.melt$variable=="best_generalization",]
    cv.res.gen.summary  <- summarySE(data=cv.res.gen, measurevar="value", groupvars=c("Lmethod"), na.rm = TRUE)
    # reorder the levels
    cv.res.gen.summary <- cv.res.gen.summary[ind.order,]
    
    res                 <- list()
    res$split           <- res.split
    res$melt            <- cv.res.melt
    res$summary.emp     <- cv.res.emp.summary
    res$summary.gen     <- cv.res.gen.summary
    
    return(res)
  }
  else
  {
    warning(paste("mergeMeltBestScoreCV: The score:", score, "is not a valid one."))
    return(NULL)
  }
}




#' Merge a list of Scores form a digest results
#' 
#' @title mergeResults
#' @description mergeResults returns a list of data frames that contain the performance of each digest in the list with their sparsity.
#' @param list.results: a list of Experiment objects one for each learner used. For example, list(res.terda, res.terga, res.terbeam) 
#' @param sparsity: Sometimes a given method will have results with somehow different sparsity. This param will allow to set the catalogue of sparsity
#' @param best.k: a vector defining wether a given k should be used to set the best model selection (default:NULL).
#' @param colors: a vector defining the colors to be used in the graphics. If not specified they will be set by default. (default:NULL).
#' @param pch: a vector defining the shape of the points to be used in the graphics. If not specified they will be set by default. (default:NULL).
#' @return list of data.frames and lists
#' @export
mergeResults <- function(list.results, sparsity = NULL, penalty = 0.001, best.k = NULL, colors = NULL, pch = NULL)
{
  # sanity checks
  if(!all(unlist(lapply(list.results, isExperiment))))
  {
    stop("mergeResults: Please provide a list of experiments")
  }
  
  # if no names set them
  if(is.null(names(list.results)))
  {
    names(list.results) <- paste("Experiment",1:length(list.results))
  }
  
  e.nb <- length(list.results)
  e.name <- names(list.results)
  
  # if the best.k is set
  if(!is.null(best.k))
  {
    if(length(best.k) != length(list.results))
    {
      stop("mergeResults: Please provide a vector best.k of the same length as there are experiments. 
           When the best.k is to be disabled for a given experiment, please provide NA for that one!")
    }
    # if(any(!is.numeric(best.k)))
    # {
    #   stop("mergeResults: Please provide a list best.k that are numeric!")
    # }
  }else
  {
    best.k <- rep(NA, e.nb)
    names(best.k) <- e.name
  }
  
  # if the best.k is set
  if(!is.null(penalty))
  {
    if(length(penalty) != e.nb)
    {
      stop("mergeResults: Please provide a vector of penalty of the same length as there are experiments. 
           When the penalty is to be disabled for a given experiment, please provide 0 for that one!")
    }
    if(any(!is.numeric(penalty)))
    {
      stop("mergeResults: Please provide a list penalty that are numeric!")
    }
  }else
  {
    penalty <- rep(0, e.nb) 
    names(penalty) <- e.name
  }
  
  # if the best.k is set
  if(!is.null(colors))
  {
    if(length(colors) != e.nb)
    {
      stop("mergeResults: Please provide a vector of colors of the same length as there are experiments.")
    }
  }else
  {
    # if not specified
    colors <- c(brewer.pal(9,"Set1"),brewer.pal(9, "Set3"), brewer.pal(9, "Spectral"),brewer.pal(8, "Dark2"))[seq_len(e.nb)]
    names(colors) <- e.name
  }
  
  # if the best.k is set
  if(!is.null(pch))
  {
    if(length(pch) != e.nb)
    {
      stop("mergeResults: Please provide a vector of pch of the same length as there are experiments.")
    }
  }else
  {
    # if not specified
    pch <- rep(19, e.nb)
    names(pch) <- e.name
  }
  
  # digest all the elements for the empirical learners
  list.results.digest <- list()
  sparsity.all <- c()
  for(i in 1:length(list.results))
  {
    list.results.digest[[i]]  <- digest(obj = list.results[[i]], penalty = penalty[i], best.k = best.k[i], plot = FALSE)
    sparsity.all              <- c(sparsity.all, list.results[[i]]$classifier$params$sparsity)
  }
  # to close set the NULL as elements if they are at the end of the list, add and omit one element
  list.results.digest[[i+1]]  <- NA; list.results.digest[[i+1]] <- NULL
  names(list.results.digest)  <- e.name # give the names
  
  # omit the empty ones (not valid) and restrict results to the valid ones;
  ind.tokeep                  <- !unlist(lapply(list.results.digest, is.null))
  list.results                <- list.results[ind.tokeep]
  list.results.digest         <- list.results.digest[ind.tokeep]
  penalty                     <- penalty[ind.tokeep]
  best.k                      <- best.k[ind.tokeep]
  colors                      <- colors[ind.tokeep]
  pch                         <- pch[ind.tokeep]
  e.nb                        <- sum(ind.tokeep)
  e.name                      <- e.name[ind.tokeep]
  
  # get the sparsities found in the data
  sparsity.all                <- sort(unique(sparsity.all))
  sparsity.full               <- seq_len(max(sparsity.all, na.rm = TRUE))
  
  # fix k_catalogue
  if(!is.null(sparsity))
  {
    k <- sparsity
  }else
  {
    k <- sparsity.full
  }
  k_catalogue <- paste("k", k, sep="_") # set a catalogue of k_sparse individuals
  
  
  # CREATE AN OBJECT CONTAINING SCORE RESULTS
  list.scores <- list()
  list.scores$empirical <- list()
  
  # EMPIRICAL long version of the best_models scores by method and variable (k_sparse)
  list.scores$empirical$fit_        <- mergeMeltScoreEmpirical(list.results.digest, k_catalogue, score = "fit_")
  list.scores$empirical$auc_        <- mergeMeltScoreEmpirical(list.results.digest, k_catalogue, score = "auc_")
  list.scores$empirical$accuracy_   <- mergeMeltScoreEmpirical(list.results.digest, k_catalogue, score = "accuracy_")
  list.scores$empirical$recall_     <- mergeMeltScoreEmpirical(list.results.digest, k_catalogue, score = "recall_")
  list.scores$empirical$precision_  <- mergeMeltScoreEmpirical(list.results.digest, k_catalogue, score = "precision_")
  list.scores$empirical$f1_         <- mergeMeltScoreEmpirical(list.results.digest, k_catalogue, score = "f1_")
  list.scores$empirical$cor_        <- mergeMeltScoreEmpirical(list.results.digest, k_catalogue, score = "cor_")
  
  
  # CROSS VALIDATION results per k_sparsity
  #auc
  list.scores$cv.empirical.auc      <- mergeMeltScoreCV(list.results.digest, k_catalogue, generalization = FALSE, score = "auc_")
  list.scores$cv.generalization.auc <- mergeMeltScoreCV(list.results.digest, k_catalogue, generalization = TRUE, score = "auc_")
  # accuracy
  list.scores$cv.empirical.acc      <- mergeMeltScoreCV(list.results.digest, k_catalogue, generalization = FALSE, score = "accuracy_")
  list.scores$cv.generalization.acc <- mergeMeltScoreCV(list.results.digest, k_catalogue, generalization = TRUE, score = "accuracy_")
  # recall
  list.scores$cv.empirical.rec      <- mergeMeltScoreCV(list.results.digest, k_catalogue, generalization = FALSE, score = "recall_")
  list.scores$cv.generalization.rec <- mergeMeltScoreCV(list.results.digest, k_catalogue, generalization = TRUE, score = "recall_")
  # precision
  list.scores$cv.empirical.prc      <- mergeMeltScoreCV(list.results.digest, k_catalogue, generalization = FALSE, score = "precision_")
  list.scores$cv.generalization.prc <- mergeMeltScoreCV(list.results.digest, k_catalogue, generalization = TRUE, score = "precision_")
  # f1-score
  list.scores$cv.empirical.f1s      <- mergeMeltScoreCV(list.results.digest, k_catalogue, generalization = FALSE, score = "f1_")
  list.scores$cv.generalization.f1s <- mergeMeltScoreCV(list.results.digest, k_catalogue, generalization = TRUE, score = "f1_")
  # correlation
  list.scores$cv.empirical.cor      <- mergeMeltScoreCV(list.results.digest, k_catalogue, generalization = FALSE, score = "cor_")
  list.scores$cv.generalization.cor <- mergeMeltScoreCV(list.results.digest, k_catalogue, generalization = TRUE, score = "cor_")
  
  
  # BEST MODEL scores regardless of k_sparse
  list.scores$best.auc              <- mergeMeltBestScoreCV(list.results.digest, k_catalogue, score="auc_", penalty = penalty)
  list.scores$best.acc              <- mergeMeltBestScoreCV(list.results.digest, k_catalogue, score="accuracy_", penalty = penalty)
  list.scores$best.rec              <- mergeMeltBestScoreCV(list.results.digest, k_catalogue, score="recall_", penalty = penalty)
  list.scores$best.prc              <- mergeMeltBestScoreCV(list.results.digest, k_catalogue, score="precision_", penalty = penalty)
  list.scores$best.f1s              <- mergeMeltBestScoreCV(list.results.digest, k_catalogue, score="f1_", penalty = penalty)
  list.scores$best.cor              <- mergeMeltBestScoreCV(list.results.digest, k_catalogue, score="cor_", penalty = penalty)
  
  
  # Make a final object to return
  list.scores$list.results.digest   <- list.results.digest
  list.scores$k_catalogue           <- k_catalogue
  
  # add the majoritary distribution. We compute this for all experiments in the case that we compare different experiments
  maj.class <- rep(NA,length(list.results))
  for(i in 1:length(list.results))
  {
    if(!is.null(list.results[[i]]$classifier$data))
    {
      maj.table <- table(list.results[[i]]$classifier$data$y)  
    }else
    {
      maj.table <- table(list.results[[i]]$classifier$y)
    }
    
    maj.class[i] <- max(maj.table/sum(maj.table))
  }
  
  list.scores$maj.class <- maj.class
  list.scores$penalty   <- penalty
  list.scores$best.k    <- best.k
  list.scores$colors    <- colors
  list.scores$pch       <- pch
  
  return(list.scores)
}



#' Merge a list of cross validation scores form digest results
#' 
#' @import reshape2
#' @import ggplot2
#' @title mergeMeltImportanceCV
#' @description mergeMeltImportanceCV returns a list of data frames that contain the feature importance of the different learners without any focus on sparsity.
#' @param list.results.digest: a list of digest objects one for each learner used. For example, list(res.terda.digest, res.terga.digest, res.terbeam.digest) 
#' @param filter.cv.prev: filter variable for each learner based on the appearence prevalence in the cross validation.
#' @param min.kfold.nb: wether we should restrict all experiments in the smallest number of k-folds of a comparative analyses (default = FALSE)
#' @param type: the type of importance "mda (mean decreased accuracy)" or "pda (prevalence decreased accuracy)" (default = mda) 
#' @param learner.grep.pattern: select a subset of learners using a grep pattern (default:"*")
#' @param nb.top.features: the number of top features to focus on the plot
#' @param feature.selection: the names of the features to be selected (default:NULL)
#' @param fixed.order: if the order of features in the plot should follow the feature selection one (default = FALSE)
#' @param scaled.importance: the scaled importance is the importance multipied by the prevalence in the folds. If (default = TRUE) this will be used, the mean mda 
#' will be scaled by the prevalence of the feature in the folds and ordered subsequently 
#' @param make.plot: make a plot for all the learners
#' @param main: should add the title to the graph for correct alignment (default:FALSE)
#' @param cv.prevalence: wether or not to plot the distribution of the prevalence of the feature in the top-models for each k-fold in the graph (default:FALSE)
#' @return a list of several data.frames and a ggplot object
#' @export
mergeMeltImportanceCV <- function(list.results, filter.cv.prev = 0.5, min.kfold.nb = FALSE,
                                  type = "mda",
                                  learner.grep.pattern = "*", 
                                  nb.top.features = 25, feature.selection = NULL, fixed.order = FALSE, 
                                  scaled.importance = TRUE, 
                                  make.plot = TRUE, main = FALSE, 
                                  cv.prevalence = TRUE)
{
  # sanity checks
  valid.efip <- unlist(lapply(list.results, function(x){!is.null(x$crossVal$fip)}))
  if(!any(valid.efip))
  {
    warning(paste("mergeMeltImportanceCV: Some learners does not have feature importance data!",paste(names(valid.efip)[valid.efip], collapse = ", ")))
  }
  
  if(length(which(!valid.efip)) == length(list.results))
  {
    warning("mergeMeltImportanceCV: stopping since thre were no data on feature importance for none of the learners")
    return(NULL)
  }
  
  if(filter.cv.prev < 0 | filter.cv.prev > 1)
  {
    stop("mergeMeltImportanceCV: please provide a valid filter prevalence value percentage [0, 1]")
  }
  
  # restrict to the valid efip learners
  list.results      <- list.results[valid.efip]
  
  e.nb              <- length(list.results) # number of classifiers
  e.name            <- names(list.results) # names of the classifiers
  
  # get the number of k-folds for each experiment, which can be differnt from one to another.
  nbfolds           <- rep(NA, e.nb)
  names(nbfolds)    <- e.name
  for(i in 1:e.nb)
  {
    nbfolds[i]      <- ncol(list.results[[i]]$crossVal$scores$empirical.auc)
  }
  
  inds.folds        <- NULL
  if(min(nbfolds) != max(nbfolds) & min.kfold.nb)
  {
    inds.folds <- list()
    for(i in 1:e.nb)
    {
      set.seed(1234)
      inds.folds[[i]]           <- sample(x = seq_len(nbfolds[i]), size = min(nbfolds), replace = FALSE)
    }
  }else
  { # else just get the increasing index
    inds.folds      <- list()
    for(i in 1:e.nb)
    {
      inds.folds[[i]]           <- seq_len(nbfolds[i])
    }
  }
  
  # build a list with results on each classifier.
  res.split         <- list()
  res.split.fprev   <- list()
  for(i in 1:e.nb) # for each classifier
  {
    if(!is.null(feature.selection))
    {
      ind <- match(feature.selection, rownames(list.results[[i]]$crossVal$fip$mda))
      if(any(is.na(ind)))
      {
        print("mergeMeltImportanceCV: the following features  are not found in the importance data... skiping")
        print(feature.selection[is.na(ind)])
        if(sum(is.na(ind)) == length(feature.selection)) 
        {
          warning(paste("mergeMeltBestScoreCV: no feature found, returning NULL"))
          return(NULL)
        }
        # features <- rownames(list.results[[i]]$crossVal$fip$mda[features, inds.folds[[i]]])
        features <- feature.selection[!is.na(ind)]
      }else
      {
        features <- feature.selection
      }
    }else
    {
      features <- rownames(list.results[[i]]$crossVal$fip$mda)
    }
    
    # for each classifier
    if(type == "mda")
    {
      imp.data                    <- list.results[[i]]$crossVal$fip$mda[features, inds.folds[[i]]]  
    }
    
    if(type == "pda")
    {
      imp.data                    <- list.results[[i]]$crossVal$fip$pda[features, inds.folds[[i]]]  
    }
    
    if(type != "pda" & type != "mda")
    {
      stop("mergeMeltImportanceCV: please provide a valid type: either pda or mda")
    }
    
    # filter by feature prevalence in the crossval final populations
    imp.data.prev               <- rowSums(!is.na(imp.data)) / ncol(imp.data)
    imp.data.prev               <- imp.data.prev[order(imp.data.prev, decreasing = TRUE)]
    # filter to 50 % of the crossval
    if(any(imp.data.prev/ncol(imp.data) > filter.cv.prev))
    {
      imp.data.prev             <- imp.data.prev[imp.data.prev / ncol(imp.data) >= filter.cv.prev]  
    }else{
      print("mergeMeltImportanceCV: filter.cv.prev too strong, canceling the filtering step...")
    }
    
    imp.data                    <- imp.data[names(imp.data.prev),]
    colnames(imp.data)          <- paste("fold_", seq_len(ncol(imp.data)), sep="")
    
    # scale the mda by mutiplying each value of by the prevalence in the folds
    imp.data.sc <- imp.data * imp.data.prev
    
    if(scaled.importance)
    {
      # use the scaled importance
      imp.data <- imp.data.sc
    }
    
    imp.data.melt <- melt(data.frame(feature = rownames(imp.data), imp.data), measure.vars = colnames(imp.data))
    colnames(imp.data.melt)     <- c("feature","folds","value")
    res.split[[i]]              <- data.frame(imp.data.melt)
    res.split.fprev[[i]]        <- imp.data.prev
  }
  names(res.split)              <- e.name
  names(res.split.fprev)        <- e.name
  
  
  if(all(unlist(lapply(res.split, nrow)) == 0))
  {
    warning(paste("mergeMeltBestScoreCV: all the results are empty, returning NULL"))
    return(NULL)
  }
  
  if(any(unlist(lapply(res.split, nrow)) == 0))
  {
    warning(paste("mergeMeltBestScoreCV: some results are empty, omitting them"))
    res.split                   <- res.split[-which(unlist(lapply(res.split, nrow)) == 0)]
  }
  
  # melt the list together
  if(length(grep(learner.grep.pattern, names(res.split))) == 0)
  {
    warning(paste("mergeMeltBestScoreCV: pattern",learner.grep.pattern,"did not find any results"))
    return(NULL)
  }
  
  # PREVALENCE data. 
  # Similarly melt the data on prevalence on the folds, which will be added to the same plot
  # select a subset of learners based on pattern and melt them together
  if(length(list.results) == 1)
  {
    fprev.melt                    <- melt(res.split.fprev[grep(learner.grep.pattern, names(res.split.fprev))], level="feature")
    fprev.melt                    <- data.frame(rownames(fprev.melt), fprev.melt)
    # rename melted data frame
    colnames(fprev.melt)          <- c("feature","value","method")
    fprev.melt$method             <- factor(fprev.melt$method, levels = names(list.results))
    # Transform in percentage
    fprev.melt$value              <- fprev.melt$value * 100
    
  }else
  {
    res.split.fprev               <- lapply(res.split.fprev, as.data.frame)
    for(i in 1:length(res.split.fprev))
    {
      res.split.fprev[[i]]        <- data.frame(feature = rownames(res.split.fprev[[i]]), res.split.fprev[[i]])
    }
    fprev.melt                    <- melt(res.split.fprev[grep(learner.grep.pattern, names(res.split.fprev))], level="feature")
    fprev.melt                    <- fprev.melt[,-2]
    # rename melted data frame
    colnames(fprev.melt)          <- c("feature","value","method")
    fprev.melt$method             <- factor(fprev.melt$method, levels = names(list.results))
    # Transform in percentage
    fprev.melt$value              <- fprev.melt$value * 100
  }
  
  
  # IMPORTANCE data
  # select a subset of learners based on pattern and melt them together
  cv.res.melt                   <- melt(res.split[grep(learner.grep.pattern, names(res.split))], level="feature")
  # rename melted data frame
  colnames(cv.res.melt)         <- c("feature","folds","variable","value","method")
  # omit a column that is not needed
  cv.res.melt                   <- cv.res.melt[,-3]
  # reset the factor order
  cv.res.melt$method            <- factor(cv.res.melt$method, levels = names(list.results))
  # Transform in percentage
  cv.res.melt$value             <- cv.res.melt$value * 100
  # summarize by feature and method to prepare for the graphic
  cv.res.summary                <- summarySE(data=cv.res.melt, measurevar="value", groupvars=c("feature","method"), na.rm = TRUE)
  # reorder features by value
  cv.res.summary.ord            <- summarySE(data=cv.res.melt, measurevar="value", groupvars=c("feature"), na.rm = TRUE)
  
  # Set the order of the features in the plot
  if(fixed.order)
  {
    # get the order by feature.selection
    lev.ord                       <- match(features, cv.res.summary.ord$feature)
    names(lev.ord)                <- features
  }else
  {
    # get the order by mean importance
    lev.ord                       <- order(cv.res.summary.ord$value, decreasing = TRUE)
    names(lev.ord)                <- as.character(cv.res.summary.ord$feature)  
  }
  
  # get the names of the top features
  if(is.null(nb.top.features))
  {
    nb.top.features <- NA
  }
  
  # if the feature mask is specified than we can not select a subset based on number
  if(!is.null(feature.selection) & !is.null(nb.top.features))
  {
    nb.top.features <- NA
  }
  
  features.tokeep               <- as.character(cv.res.summary.ord$feature[lev.ord])[1:min(nb.top.features, length(lev.ord), na.rm = TRUE)]
  # restrict data to these top features
  cv.res.summary                <- cv.res.summary[as.character(cv.res.summary$feature) %in% features.tokeep, ]
  # reorder the factor by the top
  cv.res.summary$feature        <- factor(cv.res.summary$feature, levels = rev(features.tokeep))
  
  cv.res.summary$sign           <- as.character(getSign(X = list.results[[1]]$classifier$data$X[as.character(cv.res.summary$feature),], 
                                                        y = list.results[[1]]$classifier$data$y))
  
  # same for prevalence
  # restrict data to these top features
  fprev.melt                    <- fprev.melt[as.character(fprev.melt$feature) %in% features.tokeep, ]
  # reorder the factor by the top
  fprev.melt$feature            <- factor(fprev.melt$feature, levels = rev(features.tokeep))
  
  res                           <- list()
  res$split                     <- res.split
  res$melt                      <- cv.res.melt
  res$summary                   <- cv.res.summary
  res$fprev                     <- fprev.melt
  
  if(make.plot)
  {
    pd <- position_dodge(0.3) # move them .05 to the left and right
    
    # fix color space
    col <- c("deepskyblue1","firebrick1")
    col.n <- c("-1","1")
    tab.v <- table(cv.res.summary$sign)
    if(length(tab.v) < 2) col <- col[col.n %in% names(tab.v)]
    
    if(cv.prevalence)
    {
      g <- ggplot(cv.res.summary, aes(x=feature, y=value)) + 
        # the prevalence data on the bottom
        geom_bar(data = fprev.melt, stat="identity", position="identity", aes(fill = "1")) +
        geom_hline(yintercept = min(0, cv.res.summary$value, na.rm = TRUE), col="gray") +
        ylab("Feature importance & prevalence (CV)") +
        xlab("") +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        theme_bw() +
        coord_flip() +
        scale_fill_manual("Dataset", values = c("gray90","gray90")) +
        scale_color_manual("Dataset", values = col) +
        # the importance data 
        geom_errorbar(aes(ymin = value - se, ymax = value + se, color = sign), width=.1, position=pd) +
        #geom_line(aes(group = feature, color = sign), position=pd) +
        geom_point(position = pd, size=2, shape=19, aes(color = sign)) + # 21 is filled circle
        guides(colour = "none", fill = "none")
    }else
    {
      g <- ggplot(cv.res.summary, aes(x=feature, y=value, color = sign)) + 
        # the prevalence data on the bottom
        #geom_bar(data = fprev.melt, stat="identity", position="identity", aes(fill = "1", color = "1")) +
        geom_hline(yintercept = min(0, cv.res.summary$value, na.rm = TRUE), col="gray") +
        # the importance data 
        geom_errorbar(aes(ymin = value - se, ymax = value + se), width=.1, position=pd) +
        #geom_line(aes(group = feature, color = "1"), position=pd) +
        geom_point(position = pd, size=2, shape=19) + # 21 is filled circle
        ylab("Feature importance") +
        xlab("") +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        theme_bw() +
        coord_flip() +
        scale_color_manual("Dataset", values = col) +
        guides(colour = "none", fill = "none")
    }
    
    
    # if only one learner
    if(length(list.results) == 1)
    {
      # add or not the title
      if(!main) 
      {
        if(type == "mda")
        {
          g <- g + ggtitle("mda|cv")
        }
        if(type == "pda")
        {
          g <- g + ggtitle("pda|cv")
        }
      }
    }else
    {
      g <- g + facet_grid(. ~ method, scales = "free")
    }
    
    # return the graph
    res$g                         <- g
  }  
  return(res)
}


################################################################
# FILE MANIPULATION FUNCTIONS
################################################################

#' Save a population to a file
#' 
#' @description You can use this function to save a population to a file on you're disk (it will be in your working directory)
#' @param pop: The population to be saved
#' @param fileName: The name of the file were you want to save the population
#' @import yaml
#' @export
savePopulation <- function(pop, fileName, compress = TRUE)
{
  if(length(pop) > 0)
  {
    if(length(pop[[1]][[1]]) > 1) # If the population is structured by it's sparsity
    {
      # we transform it into a non-structured list of models
      pop <- modelCollectionToPopulation(pop)
    }
    
    saveResults(pop, paste(fileName, "yml", sep = "."), compress = compress)
  }
}

#' Load a population from a file
#' 
#' @description This function is used to load a population from a file on your disk (it must be in your working directory)
#' @param fileName: The name of the file were the population is stored
#' @return a population object
#' @import yaml
#' @export
loadPopulation <- function(fileName)
{
  return(yaml.load(readChar(file(fileName), file.info(fileName)$size)))
}

#' Save the results of the fit function
#' @export
saveResults <- function(fitResults, fileName, compress = TRUE)
{
  if(compress)
  {
    write(as.yaml(fitResults), gzfile(paste(fileName, "gz", sep = ".")))
  } else
  {
    write(as.yaml(fitResults), file(fileName))
  }
}

#' Load the results of a fit
#' @export
loadResults <- function(fileName)
{
  return(yaml.load(readChar(gzfile(fileName), file.info(fileName)$size)))
}


################################################################
# EXPERIMENT HANDLERS
################################################################

#' Function used to create the counter for building clf$experiment$id
make.counter <- function()
{
  count <- 0
  f <- function()
  {
    count <<- count + 1
    return(count)
  }
  return(f)
}

#' The counter for the experiment id (used in the clf builders)
counter <- make.counter()


################################################################
# Comparative results (Blaise + Edi)
################################################################

# Computing the ranks of the models that are equivalent to the best
compare.perf <- function(result, paired = FALSE, verbose = TRUE)
{
  res <- apply(result, c(1, 2), mean, na.rm=T, na.last=FALSE, verbose = FALSE)
  res.rank <- res; res.rank[rep(TRUE,length(res))] <- NA
  res.best <- res; res.best[rep(TRUE,length(res))] <- 0
  res.pval <- res.best
  
  # for each dataset
  for(i in 1:ncol(res))
  { 
    # which one is the best learner for this dataset
    best.learner <- which.max(res[,i])
    if(verbose) print(paste("The best learner is",names(best.learner),"; index:", best.learner))
    
    # copy the orger
    # ord <- 1 + length(res[,i]) - rank(res[,i], ties.method="min")
    # res.rank[,i] <- ord
    
    res.rank[,i] <- rank(-res[,i])
    
    # for each algorithm
    for(j in 1:nrow(res))
    {
      
      if(j == best.learner) # When finding the best rank put best to 1
      { 
        res.best[j, i] <- 1 
        res.pval[j, i] <- 0 
      }else # test it
      {
        pval <- NA # init pval to NA
        try(pval <- t.test(x = result[j, i, ], 
                           y = result[best.learner, i, ], 
                           paired = paired)$p.value, silent = TRUE)
        
        # add pval to the pval graph
        res.pval[j, i] <- pval
        if(!is.na(pval))
        { 
          if(pval >= 0.05) # if not significant
          {
            res.best[j, i] <- 1
          } else # if significant
          {
            res.best[j, i] <- 0
          }
        }
        #if(verbose) print(paste("algorithm",j, ", dataset", i, ", max rank =",which.max(res[,i]), ", pval =",pval))
      }
    } # end algorithm
    
  } # end dataset
  list(rank=res.rank,
       best=res.best,
       pval=res.pval)
}


# Compute the ranks of the models that are equivalent to the best
compare.perf2dim <- function(result, paired = TRUE, verbose = TRUE)
{
  res <- apply(result, c(1), mean, na.rm=T, na.last=FALSE, verbose = FALSE)
  res.rank <- res; res.rank[rep(TRUE,length(res))] <- NA
  res.best <- res; res.best[rep(TRUE,length(res))] <- 0
  res.pval <- res.best
  
  
  # which one is the best learner for this dataset
  best.learner <- which.max(res)
  if(verbose) print(paste("The best learner is",names(best.learner),"; index:", best.learner))
  
  # copy the orger
  # ord <- 1 + length(res[,i]) - rank(res[,i], ties.method="min")
  # res.rank[,i] <- ord
  
  res.rank <- rank(-res)
  
  # for each algorithm
  for(j in 1:length(res))
  {
    
    if(j == best.learner) # When finding the best rank put best to 1
    { 
      res.best[j] <- 1 
      res.pval[j] <- 0 
    }else # test it
    {
      pval <- NA # init pval to NA
      try(pval <- t.test(x = result[j, ], 
                         y = result[best.learner, ], 
                         paired = paired)$p.value, silent = TRUE)
      
      # add pval to the pval graph
      res.pval[j] <- pval
      if(!is.na(pval))
      { 
        if(pval >= 0.05) # if not significant
        {
          res.best[j] <- 1
        } else # if significant
        {
          res.best[j] <- 0
        }
      }
      #if(verbose) print(paste("algorithm",j, ", dataset", i, ", max rank =",which.max(res[,i]), ", pval =",pval))
    }
  } # end algorithm
  
  
  list(rank=res.rank,
       best=res.best,
       pval=res.pval)
}



# Returns the smallest k with an equivalent accuracy to the best model
find.best.k <- function(res)
{
  best.k <- which.max( apply(res, 1, mean, na.rm=TRUE) )
  for (i in 1:nrow(res)) 
  {
    #print(i)
    if(i==best.k) {return (i)}
    pval <- NaN
    try(pval <- t.test(as.numeric(res[i,]),as.numeric(res[best.k,]),paired=T)$p.value, silent = TRUE)
    
    if(is.nan(pval)) 
    {
      print(paste("warning : t-test unpaired in i=",i))
      try(pval <- t.test(as.numeric(res[i,]),as.numeric(res[best.k,]),paired=F)$p.value, silent = TRUE)
      
      if(is.nan(pval)) 
      {
        print(paste("warning : pas de t-test en i=",i))
        pval <- 1
      }
    }
    
    if(pval<0.01) 
    {
      return (i)
    }
  }
}



# # r = the error rate ; n = the number of examples on which the error rate is estimated; p the significant threshold
# conf.inter <- function (r, n, p=0.05)
# {
#   t = -qnorm(p/2)
#   s <- sqrt(r * (1 - r)/n)
#   ci <- 0.5/n + t * s
#   return(rbind(r - ci, r, r + ci))
# }

confInterBinomial <- function(accuracy, n, p = 0.05)
{
  # compute the statistic
  t = -qnorm(p/2)
  # compute the variance
  s <- sqrt(accuracy * (1 - accuracy)/n)
  # compute the confidence interval
  ci <- 0.5/n + t * s
  
  res <- c(accuracy - ci, accuracy, accuracy + ci)
  names(res) <- c("inf","accuracy","sup")
  
  return(res)
}


#' Select the top significant best part of the population
#' 
#' @description This function allows to select the best part of a population that is significantly not different from the best model
#' @param pop: a list of model objects
#' @param score: the attribute of the model to be used for the evaluation
#' @param p: the p-value threshold
#' @param k_penalty: the penalty to apply to the score based on the k_sparsity (default:0)
#' @param k_max: select the best population below a given threshold. If (default:0) no selection is performed.
#' @return a sub part of the population
#' @export
selectBestPopulation <- function(pop, score = "fit_", p = 0.05, k_penalty = 0, k_max = 0)
{
  if(!isPopulation(obj = pop))
  {
    stop("selectBestPopulation: please provide a valid population.")
  }
  
  pop       <- sortPopulation(pop, evalToOrder = score)
  eval      <- populationGet_X(score, toVec = TRUE, na.rm = FALSE)(pop)
  spar      <- populationGet_X("eval.sparsity", toVec = TRUE, na.rm = FALSE)(pop)
  
  # apply the k_penalty
  eval <- eval - (k_penalty * spar)
  #plot(eval, eval - (k_penalty * spar))
  
  #eval.ci   <- as.numeric(prop.test(eval[1], n = length(pop), conf.level=0.95, correct = FALSE)$conf.int)
  eval.ci   <- confInterBinomial(accuracy = eval[1], n = length(pop), p = p)
  
  # OVERRIDE  
  # for correlation use the standard error of the mean and reverse
  # if(pop[[1]]$objective == "cor")
  # {
  #   eval.th <- eval.ci["sup"] # the highest value
  #   if(k_max == 0)
  #   {
  #     ind     <- eval < eval.th
  #   }else
  #   {
  #     ind     <- (eval < eval.th) & (spar <= k_max)
  #   }
  # }else
  # {
  #   eval.th <- eval.ci["inf"] # the lowest value
  #   if(k_max == 0)
  #   {
  #     ind     <- eval > eval.th
  #   }else
  #   {
  #     ind     <- (eval > eval.th) & (spar <= k_max)
  #   }
  # }
  
  
  eval.th <- eval.ci["inf"] # the lowest value
  if(k_max == 0)
  {
    ind     <- eval > eval.th
  }else
  {
    ind     <- (eval > eval.th) & (spar <= k_max)
  }
  
  
  # eval.max  <- prop.test(eval[1], n = length(pop), conf.level=0.95, correct = FALSE)
  # eval.min  <- confInterBinomial(accuracy = eval[1], n = length(pop), p = 0.05)["inf"]
  
  ind[is.na(ind)] <- FALSE
  
  if(!isPopulation(pop[ind]))
  {
    warning("selectBestPopulation: no models were found after selection")
    return(NULL)
  }else
  {
    return(pop[ind])  
  }
}

