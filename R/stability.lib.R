
#' analyse stability of models from digest
#'
#' @description This function analyses prevalence of features of best model of different sparsity in crossval (here still k-folds)
#' @param X: dataset to classify
#' @param y: variable to predict
#' @param clf: an object containing the different parameters of the classifier
#' @param digested.result: the digest result from digest
#' @param method: wether to compute the stability of the best compared to the best in the folds (exact), or the top best (fuzzy)
#' @return an object with first a list of feature presence tables for each k_sparsity and a list of feature presence frequency
#' @export
bestModelStability <- function(X, y, clf, digested.result, method = "exact")
{

  # TODO sanity checks
  if(!(method == "exact" | method == "fuzzy"))
  {
    stop("bestModelStability: the method is unknown, should be one of c(exact/fuzzy)")
  }

  res <- list()
  # for all the k_sparse best models
  for(j in 1:length(digested.result$best$models))
  {
    if (isModel(digested.result$best$models[[j]])) {
      best.model <- digested.result$best$models[[j]]

    } else{
      best.model <- digested.result$best$models[[j]][[1]] #get the first best
    }

    acc <- best.model$accuracy_
    ci <- conf.inter(acc,dim(X)[2])
    if(clf$params$verbose) printModel(mod = best.model, method = clf$params$print_ind_method, score = "fit_")

    # convert to dense format
    best.model.dense <- modelToDenseVec(natts = nrow(X), mod = best.model)

    if(method == "exact")
    {
      k_sparse <- paste("k", best.model$eval.sparsity, sep = "_")
      # for all n-folds
      selected <- rep(0,length(digested.result$cv$nfold))
      for (i in 1:length(digested.result$cv$nfold))
      {
        mod.list <- digested.result$cv$nfold[[i]]$results[[k_sparse]]
        n <- length(mod.list)
        mod.list <- lapply(mod.list,function(x){
          if(x$fit_ > conf.inter(best.model$fit_,n))
            return(x)
        })
        mod.list <- mod.list[!unlist(lapply(mod.list, is.null))]
        for(mod in mod.list){
          #mod <- digested.result$cv$nfold[[i]]$results[[k_sparse]][[1]]
          if(!is.null(mod)){
            mod <- modelToDenseVec(natts = nrow(X), mod = mod)
            if(length(which(abs(best.model.dense)-abs(mod)==0))==nrow(X))
              selected[i] <- 1
          }
        }

      }
    }else # if fuzzy
    {
      k_sparse <- paste("k", best.model$eval.sparsity, sep = "_")
      # for all n-folds
      selected <- c()
      p<-1
      for (i in 1:length(digested.result$cv$nfold))
      {
        # get all of them
        for (s in 1:length(digested.result$cv$nfold[[i]]$results[[k_sparse]]))
        {
          mod <- digested.result$cv$nfold[[i]]$results[[k_sparse]][[s]]
          if(!is.null(mod)){
            mod <- modelToDenseVec(natts = nrow(X), mod = mod)
            if (length(which(abs(best.model.dense) - abs(mod) == 0)) == nrow(X)) {
              selected[p] <- 1
            } else{
              selected[p] <- 0
            }}else{
              selected[p] <- 0
            }
          p <- p + 1
        }
      }
    }
    res[[j]] <- selected
  } # end k_sparsity

  s <- unlist(lapply(res,function(x)sum(x)/length(x)))

  return(list(frequency=s, presence=res))
}



#' analyse stability of models from digest
#'
#' @description This function analyses prevalence of features of best model of different sparsity in crossval (here still k-folds) 
#' @param X: dataset to classify
#' @param y: variable to predict
#' @param clf: an object containing the different parameters of the classifier
#' @param digested.result: the digest result from digest
#' @param method: wether to compute the stability of the best compared to the best in the folds (exact), or the top best (fuzzy)
#' @return an object with first a list of feature presence tables for each k_sparsity and a list of feature presence frequency
#' @export
bestModelFeatureStability <- function(X, y, clf, digested.result, method = "fuzzy") 
{
  
  # TODO sanity checks
  if(!(method == "exact" | method == "fuzzy"))
  {
    stop("bestModelFeatureStability: the method is unknown, should be one of c(exact/fuzzy)")
  }
  
  res <- list()
  # for all the k_sparse best models
  for(j in 1:length(digested.result$best$models)) 
  {
    best.model <- digested.result$best$models[[j]] #get the first best
    if(clf$params$verbose) printModel(mod = best.model, method = clf$params$print_ind_method, score = "fit_")
    
    # convert to dense format
    best.model.dense <- modelToDenseVec(natts = nrow(X), mod = best.model)
    
    if(method == "exact")
    {
      k_sparse <- paste("k", best.model$eval.sparsity, sep = "_")
      # for all n-folds
      selected <- list()
      for (i in 1:length(digested.result$cv$nfold)) 
      {
        selected[[(length(selected)+1)]] <- digested.result$cv$nfold[[i]]$results[[k_sparse]][[1]]
      }
    }else # if fuzzy
    {
      k_sparse <- paste("k", best.model$eval.sparsity, sep = "_")
      # for all n-folds
      selected <- list()
      for (i in 1:length(digested.result$cv$nfold)) 
      {
        # get all of them
        for (s in 1:length(digested.result$cv$nfold[[i]]$results[[k_sparse]])) 
        {
          selected[[(length(selected)+1)]] <- digested.result$cv$nfold[[i]]$results[[k_sparse]][[s]]
          #print(s)
        }
      }
    }
    
    # delete null if any (potentially in some folds)
    best.in.nfolds <- selected[!unlist(lapply(selected, is.null))]
    #plotPopulation(best.in.nfolds, X, y)
    # transform to dense version
    best.in.folds.dense <- listOfModelsToListOfDenseVec(clf, X, y, list.models = best.in.nfolds)
    # make a matrix
    best.in.folds.dense <- t(do.call(rbind.data.frame, best.in.folds.dense))
    rownames(best.in.folds.dense) <- rownames(X)
    
    # get the best model features in the best k-folds.
    best.in.folds.dense[best.model.dense != 0,]
    #best.in.folds.dense <- which(best.in.folds.dense!=0)
    res[[j]] <- best.in.folds.dense[best.model.dense != 0,]
    
    if(k_sparse == "k_1")
    {
      res[[j]] <- t(as.matrix(res[[j]], nrow = 1, byrow = TRUE))
      rownames(res[[j]]) <- best.model$names_
    }
  } # end k_sparsity
  
  # # compute frequency
  # res[[1]] <- t(as.matrix(res[[1]], nrow = 1, byrow = TRUE))
  # rownames(res[[1]]) <- digested.result$best$models[[1]]$names_
  s <- list()
  for(i in 1:length(res))
  {
    s[[i]] <- sort(apply(res[[i]], 1 , function(x) sum(x!=0)/length(x)), decreasing = TRUE)
  }
    
  return(list(frequency=s, presence=res))
}


#' Compute the cross-validation of leave one out for test stability
#'
#' @description Compute the cross-validation emprirical and generalization scores.
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the response vector
#' @param clf: clf
#' @param lfolds: leave one out folds for cross-validation
#' @param return.all: return all results from the crossvalidation for feature stability testing
#' @return a list containing generalisation scores for each fold as well as a matrix with the mean values.
#' @export
LPO_best_models <- function(X, y, clf, p=1, lfolds=NULL, return.all=FALSE,nk=20)
{
  
  # create the folds if they are not given
  if(is.null(lfolds))
  {
    lfolds = create.folds(y, k = -p, list = TRUE, returnTrain = FALSE) 
    nfolds <- nk
  } else {
    nfolds <- length(lfolds)
  }
  
  # create the empty object structure that will be returned
  res.crossval                                <- list()
  res.crossval$nfold                          <- list()
  #res.crossval$k                              <- list()
  res.crossval$scores                         <- list()
  # res.crossval$scores$empirical.auc           <- as.data.frame(matrix(nrow=length(clf$params$sparsity), ncol=nfolds))
  # rownames(res.crossval$scores$empirical.auc) <- c(paste("k",clf$params$sparsity,sep="_"))
  # colnames(res.crossval$scores$empirical.auc) <- paste("fold",1:nfolds,sep="_")
  # # add others using the same model
  # res.crossval$scores$generalization.auc      <- res.crossval$scores$empirical.auc
  res.crossval$scores$empirical.acc           <- as.data.frame(matrix(nrow=length(clf$params$sparsity), ncol=nfolds))
  rownames(res.crossval$scores$empirical.acc) <- c(paste("k",clf$params$sparsity,sep="_"))
  colnames(res.crossval$scores$empirical.acc) <- paste("fold",1:nfolds,sep="_")
  res.crossval$scores$generalization.acc      <- res.crossval$scores$empirical.acc # accuracy
  
  res.crossval_2                               <- res.crossval
  # for each fold compute the best models
  registerDoSNOW(cl <- makeCluster(clf$params$nCores, type = "SOCK",outfile=''))
  #cl <- makeCluster(clf$params$nCores, type = "FORK",outfile='LOG.TXT')
  #registerDoParallel(cl)
  clf$params$cluster <- cl
  if(clf$params$parallelize.folds)
  {
    print("Starting cross validation in parallel")
    # execute each crossVal in //
    
    
    res.all <- foreach (i = 1:nfolds) %dorng%
    {
      
      #printClassifier(obj = clf)
      # if (clf$params$verbose) {cat("===> k-fold\t",i,"\n")}
      
      # prepare the datasets
      # training dataset
      x_train = X[,-lfolds[[i]]]
      y_train = y[-lfolds[[i]]]
      # testing dataset
      x_test = X[,lfolds[[i]]]
      y_test = y[lfolds[[i]]]
      
      runClassifier(X = x_train, y =  y_train, clf = clf)
      
    } # end of folds loop (foreach)
    
  }else # no parallel
  {
    # execute each crossval in serial
    res.all <- list()
    for (i in 1:nfolds)
    {
      #print("Starting crossvalidation not in parallel")
      if (clf$params$verbose) {cat("===> k-fold\t",i,"\n")}
      
      # prepare the datasets
      # training dataset
      x_train = X[,-lfolds[[i]]]
      y_train = y[-lfolds[[i]]]
      # testing dataset
      x_test = X[,lfolds[[i]]]
      y_test = y[lfolds[[i]]]
      
      res.all[[i]]                      <- runClassifier(X = x_train, y =  y_train, clf = clf)
      
    } # end of folds loop (for)
  } # end else parallelize.folds
  
  
  # Dispatch the results in the custom output structure
  for(i in 1:length(res.all))
  {
    
    # prepare the datasets
    # training dataset
    x_train = X[,-lfolds[[i]]]
    y_train = y[-lfolds[[i]]]
    # testing dataset
    x_test = X[,lfolds[[i]]]
    y_test = y[lfolds[[i]]]
    
    res_train <- res.all[[i]]
    
    # res_train                         <- runClassifier(X = x_train, y =  y_train, clf = clf)
    res_train.digest                    <- digestModelCollection(obj = res_train, X = x_train, clf, mmprev=FALSE)
    res_train.digest_2                  <- digestModelCollection(obj = res_train, X = x_train, clf, mmprev=TRUE)
    if(!is.null(res_train.digest))
    {
      # for all the best models of each k-sparse (create empty matrix) for auc
      # res.crossval$k$auc                <- as.data.frame(matrix(nrow=max(clf$params$sparsity), ncol=2))
      # rownames(res.crossval$k$auc)      <- c(paste("k",c(1:max(clf$params$sparsity)), sep="_")) 
      # colnames(res.crossval$k$auc)      <- c("empirical","generalization")
      # add another table for accuracy
      #res.crossval$k$acc                <- as.data.frame(matrix(nrow=max(clf$params$sparsity), ncol=2))
      # add another table for correlation
      # res.crossval$k$cor                <- res.crossval$k$auc
      
      # for all k-sparse BEST models
      for(k in 1:length(res_train.digest$best_models))
      {
        k_sparse.name       <- names(res_train.digest$best_models)[k]
        mod                 <- res_train.digest$best_models[[k_sparse.name]]
        mod.train           <- evaluateModel(mod=mod, X = x_train, y = y_train, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE, mode='train')
        mod.test            <- list() #evaluateModel(mod=mod, X=x_test, y=y_test, clf=clf, eval.all = TRUE, force.re.evaluation = TRUE, mode='test')
        scorelist           <- getModelScore(mod = mod, X = as.matrix(x_test), clf = clf, force.re.evaluation = TRUE)
        mod$score_          <- scorelist$score_
        mod$pos_score_      <- scorelist$pos_score_
        mod$neg_score_      <- scorelist$neg_score_
        
        mod.test$accuracy_  <- evaluateAccuracy(mod = mod, X = x_test, y = y_test, clf = clf)$accuracy_
        
        # Empirical fitting score
        #res.crossval$scores$empirical.auc[k_sparse.name,i]            <- mod.train$auc_
        res.crossval$scores$empirical.acc[k_sparse.name,i]            <- mod.train$accuracy_
        #res.crossval$scores$empirical.cor[k_sparse.name,i]            <- mod.train$cor_
        
        # Generalization fitting score
        #res.crossval$scores$generalization.auc[k_sparse.name,i]       <- mod.test$auc_
        res.crossval$scores$generalization.acc[k_sparse.name,i]       <- mod.test$accuracy_
        #res.crossval$scores$generalization.cor[k_sparse.name,i]       <- mod.test$cor_
        
        # store by k
        # AUC
        #res.crossval$k$auc[k_sparse.name,"empirical"]                 <- mod.train$auc_
        #res.crossval$k$auc[k_sparse.name,"generalization"]            <- mod.test$auc_
        # Accuracy
        #res.crossval$k$acc[k_sparse.name,"empirical"]                 <- mod.train$accuracy_
        #res.crossval$k$acc[k_sparse.name,"generalization"]            <- mod.test$accuracy_
        # Regression
        #res.crossval$k$cor[k_sparse.name,"empirical"]                 <- mod.train$cor_
        #res.crossval$k$cor[k_sparse.name,"generalization"]            <- mod.test$cor_
      } # end of k_sparse loop
      
      if(return.all)
      {
        res.crossval$nfold[[i]]             <- list(results = res_train$models, resultsDigest = res_train.digest)
      }
      
      # if saving move one level up
      if(!(clf$params$popSaveFile=="NULL"))
      {
        setwd("..")
      }
    }
    if(!is.null(res_train.digest_2))
    {
      # for all the best models of each k-sparse (create empty matrix) for auc
      # res.crossval$k$auc                <- as.data.frame(matrix(nrow=max(clf$params$sparsity), ncol=2))
      # rownames(res.crossval$k$auc)      <- c(paste("k",c(1:max(clf$params$sparsity)), sep="_")) 
      # colnames(res.crossval$k$auc)      <- c("empirical","generalization")
      # add another table for accuracy
      #res.crossval$k$acc                <- as.data.frame(matrix(nrow=max(clf$params$sparsity), ncol=2))
      # add another table for correlation
      # res.crossval$k$cor                <- res.crossval$k$auc
      
      # for all k-sparse BEST models
      for(k in 1:length(res_train.digest_2$best_models))
      {
        k_sparse.name       <- names(res_train.digest_2$best_models)[k]
        mod                 <- res_train.digest_2$best_models[[k_sparse.name]]
        mod.train           <- evaluateModel(mod = mod, X = x_train, y = y_train, clf = clf, eval.all = TRUE, force.re.evaluation = TRUE, mode='train')
        mod.test            <- list() #evaluateModel(mod=mod, X=x_test, y=y_test, clf=clf, eval.all = TRUE, force.re.evaluation = TRUE, mode='test')
        scorelist           <- getModelScore(mod = mod, X = as.matrix(x_test), clf = clf, force.re.evaluation = TRUE)
        mod$score_          <- scorelist$score_
        mod$pos_score_      <- scorelist$pos_score_
        mod$neg_score_      <- scorelist$neg_score_
        
        mod.test$accuracy_  <- evaluateAccuracy(mod = mod, X = x_test, y = y_test, clf = clf)$accuracy_
        
        # Empirical fitting score
        #res.crossval$scores$empirical.auc[k_sparse.name,i]            <- mod.train$auc_
        res.crossval_2$scores$empirical.acc[k_sparse.name,i]            <- mod.train$accuracy_
        #res.crossval$scores$empirical.cor[k_sparse.name,i]            <- mod.train$cor_
        
        # Generalization fitting score
        #res.crossval$scores$generalization.auc[k_sparse.name,i]       <- mod.test$auc_
        res.crossval_2$scores$generalization.acc[k_sparse.name,i]       <- mod.test$accuracy_
        #res.crossval$scores$generalization.cor[k_sparse.name,i]       <- mod.test$cor_
        
        # store by k
        # AUC
        #res.crossval$k$auc[k_sparse.name,"empirical"]                 <- mod.train$auc_
        #res.crossval$k$auc[k_sparse.name,"generalization"]            <- mod.test$auc_
        # Accuracy
        #res.crossval$k$acc[k_sparse.name,"empirical"]                 <- mod.train$accuracy_
        #res.crossval$k$acc[k_sparse.name,"generalization"]            <- mod.test$accuracy_
        # Regression
        #res.crossval$k$cor[k_sparse.name,"empirical"]                 <- mod.train$cor_
        #res.crossval$k$cor[k_sparse.name,"generalization"]            <- mod.test$cor_
      } # end of k_sparse loop
      
      if(return.all)
      {
        res.crossval_2$nfold[[i]]             <- list(results = res_train$models, resultsDigest = res_train.digest_2)
      }
      
      # if saving move one level up
      if(!(clf$params$popSaveFile=="NULL"))
      {
        setwd("..")
      }
    }
  }
  
  # should we return all
  if(return.all)
  {
    names(res.crossval$nfold)             <- paste("n", 1:nfolds, sep = "_") # they should have the same length
    names(res.crossval_2$nfold)             <- paste("n", 1:nfolds, sep = "_") # they should have the same length
  }
  
  reorderByRownamesNumeric <- function(mat)
  {
    ind <- as.numeric(gsub("k_","",rownames(mat)))
    mat <- mat[order(ind),]
  }  
  # reorder results
  res.crossval$scores$empirical.acc       <- reorderByRownamesNumeric(res.crossval$scores$empirical.acc)
  res.crossval$scores$generalization.acc  <- reorderByRownamesNumeric(res.crossval$scores$generalization.acc)

  # accuracy
  res.crossval$scores$mean.acc            <- data.frame(cbind(rowMeans(res.crossval$scores$empirical.acc, na.rm = TRUE),
                                                              rowMeans(res.crossval$scores$generalization.acc, na.rm = TRUE)))
  colnames(res.crossval$scores$mean.acc)  <- c("empirical","generalization")
  
  res.crossval_2$scores$empirical.acc       <- reorderByRownamesNumeric(res.crossval_2$scores$empirical.acc)
  res.crossval_2$scores$generalization.acc  <- reorderByRownamesNumeric(res.crossval_2$scores$generalization.acc)
  res.crossval_2$scores$mean.acc            <- data.frame(cbind(rowMeans(res.crossval_2$scores$empirical.acc, na.rm = TRUE),
                                                              rowMeans(res.crossval_2$scores$generalization.acc, na.rm = TRUE)))
  colnames(res.crossval_2$scores$mean.acc)  <- c("empirical","generalization")
  ret <- c()
  ret$res.crossval <- res.crossval
  ret$res.crossval_2 <- res.crossval_2

  stopCluster(cl)
  
  return(ret)
}


#' analyse stability of models from digest
#'
#' @description This function analyses prevalence of features of best model of different sparsity in crossval (here still k-folds) 
#' @param X: dataset to classify
#' @param y: variable to predict
#' @param clf: an object containing the different parameters of the classifier
#' @param tmp: the digested result object from digest
#' @return a list of each sparsity the frequency of each feature of empirical best model in k-folds cross validation
#' @export
AnalyseStableModels_LOO <- function(X, y, clf, tmp, loo)
{
  #loo <- LOO_best_models(X,y,clf,return.all = TRUE)
  res <- list()
  for(j in 1:length(tmp$best$models)) {
    best.model <- tmp$best$models[[j]]
    printModel(mod = best.model, method = clf$params$print_ind_method, score = "fit_")
    # |-40|+295|-304|-373|F=0.8306|K=4
    best.model.dense <-
      modelToDenseVec(natts = nrow(X), mod = best.model)
    
    best.in.nfolds  <- list()
    n <- 1
    k_sparse <- paste("k", best.model$eval.sparsity, sep = "_")
    
    for (i in 1:length(loo$nfold))
    {
      selected <-
        ceiling(length(loo$nfold[[i]]$results[[k_sparse]]) * 0.1)
      for (s in 1:selected) {
        best.in.nfolds[[n]] <- loo$nfold[[i]]$results[[k_sparse]][[s]]
        
        n <- n + 1
      }
    }
    
    best.in.nfolds <-
      best.in.nfolds[!unlist(lapply(best.in.nfolds, is.null))]
    #plotPopulation(best.in.nfolds, X, y)
    best.in.folds.dense <- listOfModelsToListOfDenseVec(clf, X, y, list.models = best.in.nfolds)
    
    best.in.folds.dense <-
      t(do.call(rbind.data.frame, best.in.folds.dense))
    rownames(best.in.folds.dense) <- rownames(X)
    
    # get the best model features in the best k-folds.
    best.in.folds.dense[best.model.dense != 0,]
    #best.in.folds.dense <- which(best.in.folds.dense!=0)
    res[[j]] <- best.in.folds.dense[best.model.dense != 0,]
    
  }
  res[[1]] <- t(as.matrix(res[[1]], nrow = 1, byrow = TRUE))
  rownames(res[[1]]) <- tmp$best$models[[1]]$names_
  s<- list()
  for(i in 1:length(res))
    s[[i]]<-sort(apply(res[[i]],1,function(x)length(which(x!=0))/length(x)),decreasing = TRUE)
  return(list(freq=s,origin=res)) 
  
}

#' getGraph 
#'
#' @description This function gets a graph of the result of analyseStableModels
#' @param X: dataset to classify
#' @param mat: AnalyseStableModels()$origin
#' @param threshold: Used to select a number of edgets. By default this is zero.
#' @return a graph
#' @export
getGraph <- function(mat, X, threshold=0){
  links <- array(0, dim=c(nrow(X),nrow(X)+1))
  rownames(links) <- rownames(X)
  colnames(links) <- c(rownames(X),'nodes')
  for(i in 1:length(mat))
  {
    stas <- apply(mat[[i]],1,function(x)length(which(x!=0))/length(x))
    for(j in 1:length(stas))
      links[names(stas[j]),'nodes'] <- links[names(stas[j]),'nodes'] + as.numeric(stas[j])
    if(i>1){
      comb <- combn(names(stas), 2)
      for(j in 1:dim(comb)[2]){
        links[comb[1,j],comb[2,j]] <- links[comb[1,j],comb[2,j]]+1
        links[comb[2,j],comb[1,j]] <- links[comb[2,j],comb[1,j]]+1
      }
    }
  }
  g <- graph.adjacency(links[,1:nrow(X)], mode="undirected", weighted=TRUE)
  V(g)$size   <- links[,'nodes']
  E(g)$width  <- E(g)$weight
  
  iso         <- V(g)[degree(g)==0]
  g2 <- delete.vertices(g, iso)
  
  iso_2       <- V(g2)[V(g2)$size==0] 
  g2 <- delete.vertices(g2, iso_2)
  
  co <- layout_nicely(g2)
  plot(g2, layout=co, vertex.label.cex=0.5)
  return(g2)
}

#' compare stability of different modeles (intra k)
#'
#' @description This function compares stability of different modeles (intra k)
#' @param X: dataset to classify
#' @param tmp: the digested result from digest
#' @return a num
#' @export
sim_intra <- function(tmp, X)
{
  n <- dim(X)[1]
  s <- list()
  #t <- list()
  for(i in 1:length(tmp$best$models)){
    s[[i]] <- tmp$best$models[[i]]$indices_
    #t[[i]] <- tmp$cv$nfold$n_1$resultsDigest$best_models[[i]]$indices_
  }
  
  c <- length(s)
  sim <- 0
  #sim2 <- 0
  for(i in 1:(c-1))
    for(j in (i+1):c){
      sim <- sim +  (length(intersect(s[[i]],s[[j]]))-length(s[[i]])*length(s[[j]])/n)/(min(length(s[[i]]),length(s[[j]]))-max(0,length(s[[i]])+length(s[[j]])-n))
      #sim2 <- sim2 +  (length(intersect(t[[i]],t[[j]]))-length(t[[i]])*length(t[[j]])/n)/(min(length(t[[i]]),length(t[[j]]))-max(0,length(t[[i]])+length(t[[j]])-n))
    }
  
  return(sim*2/(c*(c-1)))#list(sim*2/(c*(c-1)),sim2*2/(c*(c-1))))
}

conf.inter <- function (r, n)
{
  s <- sqrt(r * (1 - r)/n)
  ci <- 0.5/n + 1.96 * s
  return(r - ci)#1.96*s)
}
# r : error rate of the best classifier
# n : number of test examples
# return the upper bound of the error rate 



#' compare stability of different modeles (inter k)
#'
#' @description This function compares stability of different modeles (inter k)
#' @param X: dataset to classify
#' @param tmp: the digested result from digest
#' @return a num
#' @export
sim_inter <- function(tmp, X)
{
  n <- dim(X)[1]
  kunch <- c()
  kunch_oppo <- c()
  kalousis <- c()
  dunne <- c()
  cw <- c()
  cw_rel <- c()
  auc_ <- c()
  acc_ <- c()
  language <- list()
  language$ratio <- list()
  language$terinter <- list()
  language$bininter <- list()
  lang_cv<-list()
  best_spar <- names(tmp$best$models)
  
  # compute the model in the conf_inter
  accuracy <- tmp$best$scores$accuracy_
  Accuracy <- max(accuracy)
  ci <- conf.inter(Accuracy, dim(X)[2])
  
  for(k in best_spar){
    
    t <- list()
    coef <- list()
    auc <- c()
    acc <- c()
    lang <- c()
    for (i in 1:length(tmp$cv$nfold)) {
      mod <- tmp$cv$nfold[[i]]$resultsDigest$best_models[[k]]
      if(length(mod)==0 || mod$accuracy_ < ci ){
        t[[i]] <- NA
        coef[[i]] <- NA
        auc[i] <- NA
        acc[i] <- NA
        lang[i] <- NA    
      }else{
        t[[i]] <- tmp$cv$nfold[[i]]$resultsDigest$best_models[[k]]$indices_
        if(is.null(tmp$cv$nfold[[i]]$resultsDigest$best_models[[k]]$coeffs_)){
          coef[[i]] <- 0
        }else{coef[[i]] <- tmp$cv$nfold[[i]]$resultsDigest$best_models[[k]]$coeffs_}
        
        auc[i] <- tmp$cv$nfold[[i]]$resultsDigest$best_models[[k]]$auc_
        acc[i] <- tmp$cv$nfold[[i]]$resultsDigest$best_models[[k]]$accuracy_
        lang[i] <- tmp$cv$nfold[[i]]$resultsDigest$best_models[[k]]$language
      }
    }
    
    
    # remove NA
    t <- t[!is.na(t)]
    coef <- coef[!is.na(coef)]
    auc <- auc[!is.na(auc)]
    acc <- acc[!is.na(acc)]
    
    c <- length(t)
    if(c==0 ||  c==1){
      kunch <- c(kunch,0)
      kunch_oppo <- c(kunch_oppo ,0)
      kalousis <- c( kalousis,0)
      dunne  <- c(dunne ,0)
      cw <- c(cw,0)
      cw_rel <- c(cw_rel,0)
      auc_ <- c(auc_,0)
      acc_ <- c(acc_,0)
    }else{
      sim_1 <- 0
      sim_2 <- 0
      sim_3 <- 0
      sim_4 <- 0
      for(i in 1:(c-1))
        for(j in (i+1):c){
          a <- coef[[i]]
          b <- coef[[j]]
          oppo <- sum(unname(a[intersect(names(a),names(b))]!=b[intersect(names(a),names(b))]))
          
          inter <- intersect(t[[i]],t[[j]])
          sim_1 <- sim_1 +  (length(inter)-length(t[[i]])*length(t[[j]])/n)/(min(length(t[[i]]),length(t[[j]]))-length(t[[i]])*length(t[[j]])/n)#max(0,length(t[[i]])+length(t[[j]])-n))
          sim_2 <- sim_2 +  length(inter)/length(union(t[[i]],t[[j]]))
          sim_3 <- sim_3 + (length(setdiff(t[[i]],t[[j]]))+length(setdiff(t[[j]],t[[i]])))/n
          sim_4 <- sim_4 +  (length(inter) - oppo/2 -length(t[[i]])*length(t[[j]])/n)/(min(length(t[[i]]),length(t[[j]]))-length(t[[i]])*length(t[[j]])/n)
          
        }
      kunch <- c(kunch,sim_1*2/(c*(c-1)))
      kunch_oppo <- c(kunch_oppo,sim_4*2/(c*(c-1)))
      kalousis <- c(kalousis,sim_2*2/(c*(c-1)))
      dunne <- c(dunne,sim_3*2/(c*(c-1)))
      tab_t <- table(unlist(t))
      
      cw <- c(cw,sum(apply(tab_t,1,function(x) x*(x-1)/((c-1)*sum(tab_t)))))
      D <- sum(tab_t) %%  n
      H <- sum(tab_t) %% c
      cw_rel <- c(cw_rel,(n*(sum(tab_t)-D+sum(apply(tab_t,1,function(x) x*(x-1))))-sum(tab_t)*sum(tab_t)+D*D)/(n*(H*H+c*(sum(tab_t)-H)-D)-sum(tab_t)*sum(tab_t)+D*D))
      
      auc_ <- c(auc_,mean(auc))
      acc_ <- c(acc_,mean(acc))
    }
    if(length(t[which(lang=='ratio')])!=0){
      language$ratio[(length(language$ratio)+1):(length(language$ratio)+length(t[which(lang=='ratio')]))] <- t[which(lang=='ratio')]
    }
    if(length(t[which(lang=='terinter')])!=0){
      language$terinter[(length(language$terinter)+1):(length(language$terinter)+length(t[which(lang=='terinter')]))] <- t[which(lang=='terinter')]
    }
    if(length(t[which(lang=='bininter')])!=0){
      language$bininter[(length(language$bininter)+1):(length(language$bininter)+length(t[which(lang=='bininter')]))] <- t[which(lang=='bininter')]
    }
    lang_cv[[k]] <- table(lang)
    
  }
  rel <- function(t,n){
    tab_t <- table(unlist(t))
    c <- length(t)
    D <- sum(tab_t) %%  n
    H <- sum(tab_t) %% c
    cw_rel <- (n*(sum(tab_t)-D+sum(apply(tab_t,1,function(x) x*(x-1))))-sum(tab_t)*sum(tab_t)+D*D)/(n*(H*H+c*(sum(tab_t)-H)-D)-sum(tab_t)*sum(tab_t)+D*D)
    return(cw_rel)
  }
  
  #names(s) <-  names(tmp$best$models)
  res <- c()
  res$sparsity <- unlist(lapply(best_spar, function(x) as.numeric(unlist(strsplit(x, "k_"))[2])))
  res$kunch <- kunch
  res$kunch_oppo <- kunch_oppo
  res$kalousis <- kalousis
  res$dunne <- dunne
  res$cw <- cw
  res$cw_rel <- cw_rel
  res$lang_cv <- lang_cv
  res$cw_rel_languages <- unlist(lapply(language, function(x) rel(x,n)))
  res$language <-unname( unlist(lapply(tmp$best$models,function(x)x$language)))
  res$auc_ <- auc_
  res$acc_ <- acc_
  return(res)
}

pop_better <- function(mod.col,eval='fit_',k_penalty=0.01){
  if(isPop(mod.col)){
    pop <- mod.col
    mod.col <- listOfModels2ModelCollection(mod.col)
  }else{
    pop <- modelCollectionToPopulation(mod.col)
  }
  bests   <- getBestIndividual(mod.col,evalToOrder = eval)
  ev      <- populationGet_X(element2get = eval, toVec = TRUE, na.rm = TRUE)(bests)
  names   <- names(ev)
  ev_2    <- shift(ev,1)
  ev_2[1] <- 0.5
  names(ev_2) <- names
  new_pop <- lapply(pop,function(x){
    
    if(as.numeric(x[eval]) > unname(ev_2[paste('k',x$eval.sparsity,sep = '_')] + k_penalty)){
      return(x)
    }
  })
  new_pop      <- new_pop[unlist(lapply(new_pop,is.null))==FALSE]
  mod.res      <- listOfModels2ModelCollection(new_pop)
  return(new_pop)
}



#ggplot()+geom_point(aes(x=res$auc_,y=res$,color=res$language))
                    
f<-function(tmp,X){
  l <- lapply(tmp$best$models,function(x)x$indices_)
  all <- sort(unique(unlist(unname(l))))
  mat <- matrix(0,nrow=length(all),ncol=length(l))
  rownames(mat) <- all
  for(i in 1:length(l)){
    mat[as.character(unlist(unname(l[i]))) ,i] <- 1
  }
  #image(t(mat))
  heatmap(mat)
}                    

# load("../../2.db_cirrhose_k_species_stage1/6.metal_languages_1K_models/results.metal_spar_v2m_new1_to_30.rda")
# load("../../2.db_cirrhose_k_species_stage1/3.terga2_languages_3K_models/results.terga2.terinter_spar_1_to_30.rda")
# load("../../../data/pasolli_2016/cirrhosis_stage_1_known_species.rda")
# 
# 
# load("../../2.db_cirrhose_k_species_stage2/6.metal_languages_1K_models/results.metal_spar_v2m_new1_to_30.rda")
# load("../../../data/pasolli_2016/cirrhosis_stage_2_known_species.rda")
# 
# load("../../2.db_colorectal_k_species/6.metal_languages_1K_models/results.metal_spar_v2m_new1_to_30.rda")
# load("../../../data/pasolli_2016/colorectal_known_species.rda")
# 
# load("../../2.db_ibd_k_species/6.metal_languages_1K_models/results.metal_spar_v2m_new1_to_30.rda")
# load("../../../data/pasolli_2016/ibd_known_species.rda")
# 
# load("../../2.db_obesity_k_species/6.metal_languages_1K_models/results.metal_spar_v2m_new1_to_30.rda")
# load("../../../data/pasolli_2016/obesity_known_species.rda")
# 
# load("../../2.db_t2d_k_species/6.metal_languages_1K_models/results.metal_spar_v2m_new1_to_30.rda")
# load("../../../data/pasolli_2016/t2d_known_species.rda")
# 
# load("../../2.db_wt2d_k_species/6.metal_languages_1K_models/results.metal_spar_v2m_new1_to_30.rda")
# load("../../../data/pasolli_2016/wt2d_known_species.rda")
# X <- X_filt
# tmp <- digest(res.metal.cv.v2m_new)
# res <- sim_inter(tmp,X)
# f(tmp,X)
# #ggplot for res$lang_cv
# 
# ggplot(melt(res$lang_cv)) +aes(lang, value,fill=factor(lang),color=factor(lang))+
#   geom_bar(stat='identity') +
#   facet_wrap(~L1)
# 
# #plot for result of sim_inter
# ggplot()+geom_line(aes(x=c(1:length(res$kunch)),y=res$kunch,color=factor("kunch")))+
#   geom_line(aes(x=c(1:length(res$kunch)),y=res$kunch_oppo,color=factor('kunch_oppo')))+
#   geom_line(aes(x=c(1:length(res$kunch)),y=res$kalousis,color=factor('kalousis')))+
#   geom_line(aes(x=c(1:length(res$kunch)),y=res$dunne,color=factor('dunne')))+
#   geom_line(aes(x=c(1:length(res$kunch)),y=res$cw,color=factor('cw')))+
#   geom_line(aes(x=c(1:length(res$kunch)),y=res$acc_,color=factor('acc')))+
#   geom_line(aes(x=c(1:length(res$kunch)),y=res$auc_,color=factor('auc')))+
#   xlab('sparsity')+ylab('stability')+labs(colour = "measures")
# 
# ggplot()+geom_line(aes(x=res$auc_,y=res$kunch,color=factor("kunch")))+
#   geom_line(aes(x=res$auc_,y=res$kalousis,color=factor('kalousis')))+
#   geom_line(aes(x=res$auc_,y=res$dunne,color=factor('dunne')))+
#   geom_line(aes(x=res$auc_,y=res$cw,color=factor('cw')))+
#   geom_line(aes(x=res$auc_,y=res$acc_,color=factor('acc')))+
#   xlab('auc')+ylab('stability')+labs(colour = "measures")+
#   geom_text(aes(label = names(res$auc_)))
# 
# 
# auc<-lapply(res.metal.cv.v2m_new$classifier$models,function(x) lapply(x,function(y)y$auc_))
# 
# ind <- unname(lapply(res.metal.cv.v2m_new$classifier$models,function(x) lapply(x,function(y)y$indice)))


getGraph_features <- function(ind,X){
  library('binhf')
  dig<-digest(res.metal.cv.v2m_new)
  auc<-dig$best$scores$auc_
  name<-names(auc)
  #unname(unlist(lapply(dig$best$models,function(x)x$auc_)))
  auc<-unname(shift(auc,1))
  names(auc)<-name
  auc[1] <- 0
  ind <- list()
  score <- list()
  n   <- 1
  for(i in 1:length(res.metal.cv.v2m_new$classifier$models)){
    t<- unique(res.metal.cv.v2m_new$classifier$models[[i]])
    for(j in 1:length(t)){
      if(t[[j]]$auc_ >= auc[i]){
        ind[[n]] <- t[[j]]$coeffs_
        score[n] <- t[[j]]$accuracy_
        n<-n+1
      }
    }
  }
  
  feat <- array(0, dim=c(length(unique(names(unlist(ind)))),30))
  rownames(feat)<-unique(names(unlist(ind)))
  for(i in 1:length(ind)){
    for(j in 1:length(ind[[i]])){
      feat[names(ind[[i]])[j],length(ind[[i]])] <-feat[names(ind[[i]])[j],length(ind[[i]])] +1 
    }
  }
  heatmap.2(feat,Rowv=NA,Colv = NA)
  
  links <- array(0, dim=c(nrow(X),nrow(X)+2))
  rownames(links) <- rownames(X)
  colnames(links) <- c(rownames(X),'nodes','scores')
  for(i in 1:length(ind)){
    if(length(ind[[i]])==1){
      links[names(ind[[i]]),'nodes'] <- links[names(ind[[i]]),'nodes']+1
    }else{
    pairs <- combn(names(ind[[i]]),2)
    for(j in 1:length(pairs[1,])){
      links[pairs[1,j],pairs[2,j]]<-links[pairs[1,j],pairs[2,j]]+1
      links[pairs[2,j],pairs[1,j]]<-links[pairs[2,j],pairs[1,j]]+1
       }
    for(j in 1:length(ind[[i]])){
      links[names(ind[[i]])[j],'nodes'] <- links[names(ind[[i]])[j],'nodes']+1
      links[names(ind[[i]])[j],'scores'] <- links[names(ind[[i]])[j],'scores']+as.numeric(score[[i]])
       }
     }
  }
  links[,'scores']<-links[,'scores']/links[,'nodes']
  links[is.na(links)] <- 0
  g <- graph.adjacency(links[,1:nrow(X)], mode="undirected", weighted=TRUE)
  V(g)$size   <- links[,'scores']*5
  E(g)$width  <- E(g)$weight
  
  iso         <- V(g)[degree(g)==0]
  g2 <- delete.vertices(g, iso)
  
  iso_2       <- V(g2)[V(g2)$size==0] 
  g2 <- delete.vertices(g2, iso_2)
  # 
  iso_3       <- E(g2)[E(g2)$width< max(E(g2)$weight)-5] 
  g2 <- delete.edges(g2, iso_3)
  
  iso_4         <- V(g2)[degree(g2)==0]
  g3 <- delete.vertices(g2, iso_4)
  
  co <- layout_nicely(g2)
  co <- layout.fruchterman.reingold(g3)
  plot(g3, layout=co, vertex.label.cex=0.8)
}
# m <- melt(ind)
# names(m) <- c('features','num','k')
# lapply(1:30, function(x)m[m$k>x,])
# 
# #inds <- lapply(ind,function(x)lapply(x,function(y)lapply(y,unlist)) )
# #inds <- lapply(ind,function(x)unlist(lapply(x,function(y) unlist(y)))) 
# inds <- list()
# n <- 1
# for(i in 1:length(ind)){
#   #inds[[n]] <- 
#   for(j in 1:length(ind[[i]])){
#     
#     inds[[n]]<-ind[[i]][[j]]
#     n<-n+1
#   }
# }
# 
# d<-digest_mmprev(clf.metal.v2m_new,X,res.metal.cv.v2m_new)
