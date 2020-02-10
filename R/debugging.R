
#' @export
debugging_terda_YC <- function() {
  
  # le code ci-dessous est repris de test_cirrhose_k_species.R, de predomics_testing
  #setwd('/Users/Yann/Dropbox/Eclipse/predomics/predomics_testing/analyses/2.db_cirrhose_k_species_stage1')
  #library(momr)
  # set.seed(1)
  # 
  # load("../../data/pasolli_2016/cirrhosis_stage_1_known_species.rda")
  # features <- rownames(X)
  # features.simple <- paste("SP_",1:length(features),sep = "")
  # rownames(X) <- features.simple
  # 
  # X <- filterNoSignal(X, side = 1, threshold = "auto")
  load("testing/test_cirrhose_species_stage1.rda")
  
  # Global options
  sparsity = 1:30
  nCores = 1
  seed = 1234
  useCustomLanguage = FALSE
  experiment.save = "nothing"
  print_ind_method = "short"
  #evalToFit = "fit_"
  
  # specific terda
  terda.nbIter = length(sparsity)
  terda.nbRR = 1000
  
  clf.terda   <- terda(method = "glmnetRR", 
                       sparsity = sparsity, 
                       verbose = TRUE, 
                       nIterations = terda.nbIter, 
                       nCores = nCores, 
                       nRR = terda.nbRR, 
                       seed = seed,
                       #useCustomLanguage = useCustomLanguage,
                       experiment.id = "terda",
                       print_ind_method = print_ind_method,
                       experiment.save = experiment.save,
                       language='bininter',
                       alpha = 1)
  
  lfolds <- create.folds(y, k = 5, list = TRUE, returnTrain = FALSE, seed = seed)
  
  res.terda.cv   <- fit(X, y, clf.terda, cross.validate = TRUE, lfolds = lfolds)
  return(res.terda.cv)
}


#' @export
test_commputeIntercept <- function(){
  load("../data/testing_data/test_intercept.rda")
  # loading mod, X, y
  
  myAssert(condition = isModel(obj = mod), message = "isModel")
  
  myAssert(isclose(computeIntercept(mod$score_, y, verbose=TRUE, sign="auto")$intercept, 0.0006854003), 
           message="auto sign")
  
  myAssert(isclose(computeIntercept(mod$score_, y, verbose=TRUE, sign="<")$intercept, 0.3720596), 
           message="sign <")
  
  myAssert(isclose(computeIntercept(mod$score_, y, verbose=TRUE, sign=">")$intercept, 0.0006854003), 
           message="sign >")
  
  myAssert(isclose(computeIntercept(mod$score_[1], y[1], verbose=TRUE, sign="auto")$intercept, 2.830001e-05), 
           message="one value")
  
  myAssert(isclose(computeIntercept(mod$score_, abs(y), verbose=TRUE, sign="auto")$intercept, 0.3720596), 
           message="only positive")
  
  myAssert(isclose(computeIntercept(mod$score_, -abs(y), verbose=TRUE, sign="auto")$intercept, -1.47e-05), 
           message="only negative")
  
  print("computeIntercept: OK")
}




test_evaluateModel <- function()
{
  load("../data/testing_data/test_population.rda")
  # loading mod, X, y
  
  myAssert(condition = isModel(obj = mod), message = "isModel")
  myAssert(evaluateModel(mod, X, y, clf)$accuracy_,  0.7237569)
  myAssert(evaluateModel(mod, X, y, clf, force.re.evaluation = TRUE)$accuracy_,  0.7237569)
  myAssert(evaluateModel(mod, X, y, clf, force.re.evaluation = TRUE)$auc_,  0.7520285)
  myAssert(evaluateModel(mod, X[,-c(1:10)], y[-c(1:10)], clf, force.re.evaluation = TRUE, mode = "train")$accuracy_,  0.7237569)
  myAssert(evaluateModel(mod, X[,-c(1:10)], y[-c(1:10)], clf, force.re.evaluation = TRUE, mode = "train")$accuracy_,  0.7237569)
  
  # > evaluateModel(mod, X, y, clf)$accuracy_
  # [1] 0.6187845
  # > evaluateModel(mod, X, y, clf, force.re.evaluation = TRUE)$accuracy_
  # [1] 0.6187845
  # > evaluateModel(mod, X, y, clf, force.re.evaluation = TRUE, mode = "train")$accuracy_
  # [1] 0.6187845
  # > evaluateModel(mod, X, y, clf, force.re.evaluation = TRUE, mode = "test")$accuracy_
  # [1] 0.6187845
  
}

#debugging_terda_YC()



# #### This is a test function, maybe it should go somewhere else ?
# test_population_of_vectors <- function() {
#   #source('mainTerDa.R')
#   clf <- terda()
#   
#   l = list()
#   l[[1]] = c(0.0000000001,1.0,0.0,-5.9)
#   l[[2]] = c(3,2.0,0.0,0.0)
#   l[[3]] = c(0.0,1.5,1.0,-1.0)
#   
#   load("../data/DATAMETA1.rda")
#   y <- DATAMETA1$class
#   X <- t(DATAMETA1[,-1])
#   
#   return(populationOfVectorsToModel(clf,X,y,l))
# }



