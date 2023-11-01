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
# @author: Lucas Robin
# @author: DAO Van Sang
# @author: Yann Chevaleyre
# @date: August 2016                                                    
################################################################

######################################## CORE FUNCTIONS OF TERDA ###############

#' Solve with GLMNET and create models
#' @description Create Models by applying randomized roundings on the a solution given by GLMNET
#' @importFrom glmnet glmnet
glmnetRR <- function(clf, X, y)
{
  p <- clf$params
  check.X_y_w(X,y) #sanity check
  
  family <- "binomial"
  
  switch( tolower(p$language) ,
          "logreg" = {
            intercept    = TRUE
            lower.limits = -Inf
            upper.limits = Inf
            lambdas = NULL
          },
          "bininter" = {
            lower.limits = 0
            upper.limits = 1
            intercept    = TRUE

            lambdas = exp(-c(0:p$nblambdas)/(28*sqrt(p$nblambdas/150)))
          },
          "bin" = {
            lower.limits = 0
            upper.limits = 1
            intercept    = FALSE
            lambdas = exp(-c(0:p$nblambdas)/(28*sqrt(p$nblambdas/150)))
          },
          "terinter" = {
            lower.limits = -1
            upper.limits = 1
            intercept    = TRUE
            lambdas = exp(-c(0:p$nblambdas)/(28*sqrt(p$nblambdas/150)))
          },
          "ter" = {
            lower.limits = -1
            upper.limits = 1
            intercept    = FALSE
            lambdas = exp(-c(0:p$nblambdas)/(28*sqrt(p$nblambdas/150)))
          },
          "regression" = {
            lower.limits = -1
            upper.limits = 1
            intercept    = FALSE
            family <- "gaussian"
            lambdas = exp(-c(0:p$nblambdas)/(28*sqrt(p$nblambdas/150)))
          },
          {
            stop('Unknown language! Please verify the input parameters.')
          }
  )
  
  if(clf$params$verbose) print("... ... language set and lamdas computed")

  ## Modification par Yann: - je mets standardize=TRUE. C'etait ce qui empechait d'avoir plus d'une feature par modele
  ##                        - j'enleve  nlambda = 1000, car (peut etre un bug de glmnet), pour nlambda=100 (valeur par defaut) j'obtiens plein de coefs,
  ##                          et pour nlambda = 1000 je n'obtiens qu'un coef non nul.
  ## Nouvelle modification par Yann:
  ## Finalement, je donne une liste de lambdas : exp(-c(0:100)/28)/5
  ## cette liste correspond a peu pres aux lambdas qu'il trouve lui-meme quand il marche bien

  # nblambdas = 400
  # lambdas = exp(-c(0:nblambdas)/(28*sqrt(nblambdas/150)))
  if(clf$params$verbose) print("... ... running glmnet")
  system.time(glmmod <- glmnet(x = t(X), 
                               y = y, 
                               alpha = p$alpha, 
                               family = family, 
                               lower.limits = lower.limits, 
                               upper.limits = upper.limits, 
                               intercept = intercept, 
                               lambda = lambdas, 
                               standardize = TRUE)
              )
  if(clf$params$verbose) print("... ... glmnet object is created")
  
  if(clf$params$plot)
  {
    if(clf$params$verbose) print("... ... creating plots")
    par(mfrow=c(2,1))
    #plot(glmmod)
    plot(glmmod, xvar="lambda")
    plot(glmmod, xvar="dev")
  }
  
  #plot(glmmod); plot(glmmod, xvar="dev"); plot(glmnet, xvar="lambda")
  pop <- list()
  lobj <- list()
  # for(i in 1:length(pop))
  # {
  #   lobj[[i]] <- glmmod
  # }

  # initialize results
  glmmod.sparsity <- rep(0, length(glmmod$lambda))
  glmmod.sumcoefs <- rep(0, length(glmmod$lambda))
  glmmod.stddev   <- rep(0, length(glmmod$lambda))
  
  # for all lamdas get the coeffs
  for(j in 1:length(glmmod$lambda))
  {
    coef_j             <- coef(glmmod, s = glmmod$lambda[j])  # vecteur de coefs du jieme modele
    coef_j             <- coef_j[-1]          # changement YANN - enleve l'intercept
    glmmod.sparsity[j] <- sum(coef_j!=0)      # nombre de coefs non nuls pour chaque lambda possible
    
    # changement YANN - glmmod.sumcoefs[j] <- sum(abs(coef_j[2:nrow(coef_j)]))    # somme des val.absolue des coefs pour chaque lambda possible
    glmmod.sumcoefs[j] <- sum(abs(coef_j))    # somme des val.absolue des coefs pour chaque lambda possible
    glmmod.stddev[j]   <- 0
    if (p$nRR >= 1) 
    {
      glmmod.stddev[j]   <- sqrt(sum(abs(coef_j)*(1-abs(coef_j))))
    }
  } # end coeffs
  
  if(clf$params$verbose) print("... ... preparing coefficient calculation")

  #print(c('Sparsity of models found by GLMNet (before rounding): ',intersect(glmmod.sparsity, p$sparsity)))
  if(clf$params$debug)
  {
    cat(paste('Sparsity of models found by GLMNet (before rounding):\n',paste(unique(glmmod.sparsity), collapse = ","), "\n",sep=""))
  }

  p.sparsity.min  = min(p$sparsity)
  p.sparsity.max  = max(p$sparsity)

  lamdas.zone <- rep(FALSE,length(glmmod.sumcoefs))
  
  # for the different lambdas
  for(j in seq_along(glmmod.sumcoefs)) 
  {
    sumcoefs      = glmmod.sumcoefs[j]
    sdev          = glmmod.stddev[j]
    lambda        = glmmod$lambda[j]

    # On ne traite que les modeles dont la sparsite apres rounding a une chance d'etre dans le range de p$sparsity.
    # Note: la sparsite APRES rounding sera dans l'intervale [sumcoefs-3*sdev,sumcoefs+3*sdev] avec p>99%
    # car P[moy-3.stddev < x < moy+3.stddev ] > 99%
    # Donc on traite si [sumcoefs-3*sdev,sumcoefs+3*sdev] intersecte avec [p.sparsity.min, p.sparsity.max]
    SparsityInRange = (sumcoefs-3*sdev < p.sparsity.max) && (sumcoefs+3*sdev > p.sparsity.min) && (sumcoefs > 0)

    # Edi added this to solve a bug when SparsityInRange is NA
    #if(is.na(SparsityInRange)) SparsityInRange <- FALSE
      
    if(SparsityInRange || p$nRR==0) # if in the sparsity range or no roundings
    {
      if(clf$params$verbose & !clf$params$debug)
      {
        cat(paste(".")) # progress
      }
      
      if(clf$params$debug) # give more information when debugging
      {
        cat(paste("glmnetRR...\t","j:",j,"\tlambda:",signif(lambda,3),"\n"))
      }

      glmmod.coefs <- as.matrix(coef(glmmod, s = lambda))
      glmmod.coefs <- glmmod.coefs[2:(nrow(glmmod.coefs)),]
      
      if(!all(glmmod.coefs==0)) # if the model is not empty
      {
        # Randomized Rounding
        if(p$nRR >= 1) 
        {
          set.seed(clf$params$seed) # Different "seed" give us different result. # want different results
          res           <- multipleRR(clf = clf, X = X, y = y, w = glmmod.coefs, n = p$nRR) # use the no parallel version due to more global paralleliszaition
        }else 
        {
          res           <- list(glmmod.coefs)
        }
        
        if(clf$params$language == "logreg")
        {
          mod           <- denseVecToModel(X = X, y = y, v = res[[1]], clf = clf, obj=glmmod)
          if(!is.null(mod))
          {
            mod$lambda    <- lambda  
          }
          pop           <- c(pop, list(mod)) # add it to the pop
        }else{
          mod           <- denseVecToModel(X = X, y = y, v = res[[1]], clf = clf)
          pop           <- c(pop, list(mod)) # add it to the pop
        }
        lamdas.zone[j] <- TRUE # to plot the zone when it is zommed
      }
    }else
    {
      if(clf$params$debug) 
      {
        print(paste(j, "This model has no coeffs. lambda =",signif(lambda,4)))
        # Exclude the fact when the model has no coefficients.
      }
    }
  } # end loop for all coeffs
  
  # clean the population (NULL) models
  pop <- cleanPopulation(pop = pop, clf = clf)
  # delete the duplicated ones
  pop <- pop[!duplicated(populationGet_X("indices_", toVec = FALSE)(pop))]
  
  if(clf$params$verbose) cat("\n")
  if(clf$params$verbose) print("... ... coefficients are computed")
  
  if(clf$params$plot)
  {
    par(mfrow=c(1,1))
    plot(lambdas, col="white", pch="*", main="The lambdas yielding valid coefficients")
    abline(h=lambdas[lamdas.zone], col="gray")
    points(lambdas, col="red", pch="*")
  }
  
  if(length(pop)==0)
  {
    if(clf$params$verbose) print("... ... no model is found returning empty handed")
    return(NULL)
  }
  
  if(clf$params$verbose) print("... ... models are found and created")
  
  pop2add <- list()
  # IF TERDA mode
  if(clf$params$language != "logreg")
  {
    # FILL-up coefficients for a given sparsity
    sparsity.models   <- populationGet_X("eval.sparsity")(pop)
    #sparsity.models   <- unlist(lapply(pop, function(x){sum(x!=0)}))
    sparsity.ok       <- as.numeric(names(table(sparsity.models)))
    sparsity.missing  <- clf$params$sparsity[!clf$params$sparsity %in% sparsity.ok]
    # Create another population to be added to the GLMNET one. This will be seeded from the pop. 
    # The number of individuals to add will be that of the mean prevalence
    nb.individuals    <- round(mean(table(sparsity.models)))
    
    for(i in seq_along(sparsity.missing))
    {
      # get the current sparsity
      current.sparsity <- sparsity.missing[i]
      # compute the difference, and focus on the positive (greater than current sparsity) to draw smaller ones.
      ds <- sparsity.ok-current.sparsity
      # best matching sparsity
      bms <- sparsity.ok[min(which(ds>0))]
      
      # no higher sparsities found
      if(is.na(bms) | bms < 1)
      {
        # the feature space to draw from will be all the models
        features.todraw <- unique(populationGet_X("indices_", toVec = TRUE, na.rm = TRUE)(pop))
        #features.todraw <- unique(unlist(lapply(pop, function(x){which(x!=0)})))
      }else
      {
        # get available models from whom to draw features
        models.todraw <- pop[which(sparsity.models==bms)]
        # get available features from this sparsity
        features.todraw <- unique(populationGet_X("indices_", toVec = TRUE, na.rm = TRUE)(models.todraw))
        #features.todraw <- unique(unlist(lapply(models.todraw, function(x){which(x!=0)})))
      }
      
      # create new individuals
      for(j in 1:nb.individuals)
      {
        set.seed(clf$params$current_seed) + i + j
        pop2add <- c(pop2add, 
                     list(sample(features.todraw, size = current.sparsity, replace = TRUE)))
      }
    } # end of adding individuals
    
    if(clf$params$verbose) print("... ... created new models with missing sparsities")
  }
  
  # transform the population of sparse vectors in list of models
  if(length(pop2add) > 0)
  {
    pop2add.mod <- listOfSparseVecToListOfModels(X = X, y = y, clf = clf, v = pop2add)
    pop.final <- c(pop, pop2add.mod)
  }else
  {
    pop.final <- pop
  }
  
  # clean the population (NULL) models
  pop.final <- cleanPopulation(pop = pop.final, clf = clf)
  # delete the duplicated ones
  pop.final <- pop.final[!duplicated(populationGet_X("indices_", toVec = FALSE)(pop.final))]
  pop.final <- evaluatePopulation(X = X, y = y, clf = clf, pop = pop.final, eval.all = TRUE, force.re.evaluation = TRUE)
  
  # convert the final population to a model collection object
  model.collection <- listOfModels2ModelCollection(pop = pop.final)
  
  if(clf$params$verbose) print("... ... converting to modelCollection")
  
  return(model.collection)
}




#BUG ICI:
#> source('mainTerDa.R')
#> clf <- terda(nIterations=10,method="myRRSkbest",sparsity=c(3,15))
#> res.terda <- terda_fit(X,y,clf)
#[1] "sparsity ....> "
#[1] 15
#[1] "sparsity ....> "
#[1] 3
#Show Traceback
#
#Rerun with Debug
#Erreur dans if ((!remove.zero.vec) | sum(abs(wrr)) > 0) { : 
#                                                            valeur manquante là où TRUE / FALSE est requis 

myRRSkbest <- function(clf, X, y) 
{
  p <- clf$params
  check.X_y_w(X,y) #sanity check
  tX = as.data.frame(t(X))
  
  #dim_x2 = dim(tX)[2] # nombre des attributes
  natts = dim(tX)[2] # nombre des attributes
  fidx <- NULL # les index des attibuts qu'on va recalculer ou nul si on recalcule
  fval <- NULL # les coefficients des attibuts avec un poid non nuls
  best_dense_vecs <- list()
  # for all the iterations
  for (i in 1:p$nIterations) {
    # fix sparsity
    if (length(p$sparsity) == 1) {  # update Edi
      sparsity <- p$sparsity
      if(clf$params$verbose) print(paste( "sparsity fixed ....> " , sparsity ))
    }else { 
      sparsity <- sample(p$sparsity, 1) 
      if(clf$params$verbose) print(paste( "sparsity random ....> " , sparsity ))
    }
    
    # TODO fix k_sparse 
    
    # run Solver
    #sparsity = p$sparsity #commented by Edi
    # build linear program
    l <- buildlp(tX, y, k.sparse = sparsity, vartype = "real", gamma = p$gamma, lb = p$lb, ub = p$ub, fcoefidx = fidx, fcoefval = fval)
    wb <- runsolver(l, tX)
    fidx <- c(100:nrow(X))
    fval <- rep(0,length(fidx))
    
    fidx <- NULL
    fval <- NULL
    sparsity=7
    wb <- buildAndSolveReducedLp(x = tX, y = y, k.sparse = sparsity, vartype = "real", gamma = p$gamma, mytimeout = 60, lb = p$lb, ub = p$ub, fcoefidx = fidx, fcoefval = fval)
    wb <- wb[-length(wb)] # enlever l'offset
    #if (length(wb) > 0) { # parfois le solveur ne marche pas :()
    #       if(clf$params$verbose){
    #         print(paste("solver did not find a model in iteration =",i))
    #       }
    res <- multipleRR(clf = clf, X = X, y = y, w = wb, n = p$nRR) # appelle nRR fois le rounding
    # TODO (shouldn't we get the bestmodels for all the modelCollection instead ?) 
    # When we explore a range of k_sparse the best of the collection will tend to be a high number. 
    # By selecting the best per k_sparse we get a larger distribution (more models)
    #UNIQUE
    #       mod <- getTheBestModel(res)
    #       #getBestModels()
    #       ter_wb <- modelToDenseVec(natts, mod)
    #       best_dense_vecs[[length(best_dense_vecs)+1]] <- ter_wb
    
    # ALL the best per each k_sparse /// Update Edi
    #best.models <- getBestModels(res) # all the best models for each k_sparse
    best.models <- getNBestModels(obj = res,
                                  significance = TRUE, 
                                  by.k.sparsity = TRUE,
                                  k.penalty = 0,
                                  n.best = 1,
                                  single.best = FALSE,
                                  single.best.cv = FALSE,
                                  single.best.k = NULL,
                                  max.min.prevalence = TRUE,
                                  X = NULL,
                                  verbose = FALSE, 
                                  evalToOrder = "unpenalized_fit_",
                                  return.population = TRUE # MC 
    )
    ter_wbs <- listOfModelsToListOfDenseVec(clf = clf, X = X, y = y, list.models = best.models)
    
    the.best.model <- getTheBestModel(res) # the best model
    ter_wb  <- modelToDenseVec(natts = natts, mod = the.best.model)
    # save best models of a given collection
    best_dense_vecs <- c(best_dense_vecs, ter_wbs)
    #    }
    # selectionner des attributs au hasard, dont on fixe les valeurs a cette iteration
    fidx = sample(c(1:natts), replace=FALSE, size=round(natts - p$kBest))
    fval = ter_wb[fidx]
  }
  
  # keep only unique models
  best_dense_vecs <- unique(best_dense_vecs)
  
  #if(clf$params$verbose) print(best_dense_vecs)
  return(models = listOfDenseVecToModelCollection(clf = clf, X = X, y = y, v = best_dense_vecs)) 
  #return(  list(clf=clf,models=listOfDenseVecToModelCollection(clf = clf,X = X,y = y,v = best_dense_vecs )  )  )
}


#' multipleRR
#'
#' @description computes multiple randomized rounding for a given vector of wi
#' @param clf: the classifier parameter object
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the response vector
#' @param w: a vector of wi coefficients
#' @param n: number of round roundings to compute
#' @param remove.zero.vec: whether to remove the zero vectors
#' @return an population of dense vectors
multipleRR <- function(clf, X, y, w, n, remove.zero.vec = TRUE){
  
  if (length(w)==0) 
  {
    stop("multipleRR: rounding vectors of lenght 0 does not make sense")
  }
  
  # make list of randomized rounded coef vectors
  listW <- list()

  # MODIF YANN: Rounder uniquement les coefs qui ne sont pas deja a zero ou un
  #idxn1 = which(w==-1)
  idx0 = which(w==0)
  idx1 = which(w==1)
  w_reduced <- w[-c(idx0,idx1)]
  #w_reduced <- w[-c(idxn1,idx1)]
  # END MODIF YANN
    
  for(i in 1:n)
  {
    # MODIF YANN
    # AVANT, IL Y AVAIT :wrr <- terDA.rr(w)
    if(length(w_reduced) > 0)
    {
      wrr_reduced <- terDA.rr(w_reduced)
    }
    wrr <- w
    if(length(w_reduced) > 0)
    {
      wrr[-c(idx0, idx1)] <- wrr_reduced
    }
    # FIN MODIF YANN    
    
    if(clf$params$debug)
    {
      cat(paste("\tmultipleRR...\t","rounding:",i,"\n"))
    }
    
    if ((!remove.zero.vec) | sum(abs(wrr)) > 0) 
    {
      listW[[ length(listW)+1 ]] <- wrr
    }
  } # end for loop
  
  if(length(listW)==0)
  {
    warning("No model found with this configuration") #updated Edi
    return(NULL)
  }else
  {
    #models <- listOfDenseVecToModelCollection(clf = clf, X = X, y = y, v = listW)
    models <- listW # we don't transform this in a model collection yet since no needed here
    return(models) #updated Edi
  }
}


#' multipleRR_par
#'
#' @description computes in parallel multiple randomized rounding for a given vector of wi
#' @param clf: the classifier parameter object
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the response vector
#' @param w: a vector of wi coefficients
#' @param n: number of round roundings to compute
#' @param remove.zero.vec: whether to remove the zero vectors
#' @return an population of dense vectors
multipleRR_par <- function(clf, X, y, w, n, remove.zero.vec = TRUE)
{
  sparsite = sum(w != 0)
  #stopifnot(sparsite>0)
  myAssert(sparsite>0, message = "multipleRR_par: null sparsity!",stop = TRUE)
  # si w est de sparsite 1, par exemple w=(0,0,0.2,0,0,...,0) alors le seul vecteur a explorer est (0,0,1,0,0,...,0).
  # Donc inutile de faire des tas de RR
  if(sparsite == 1) 
  {
    wi = which(w!=0)
    if (w[wi] > 0) # if positive coefficient
    {
      w[wi] = 1.0
    } else # if negative coefficient
    {
      w[wi] = -1.0
    }
    models <-  multipleRR(clf, X, y, w, 1)
  }
  
  # if par is False, then call the non-parallel version of RR
  if(!clf$params$parallel) # if not in parallel
  {
    models <-  multipleRR(clf, X, y, w, n) # run multipleRR not in //
  } else 
  {
    if (length(w)==0) 
    {
      stop("multipleRR_par: trying randomized rounding (function multipleRR) with vectors of lenght 0")
    }
    
    # make list of randomized rounded coef vectors
    tmpW <- foreach(i = 1:n, .combine = cbind, .export = "terDA.rr") %dorng% 
    {
      terDA.rr(w) # Prend du temps !
    }
    
    listW <- list()
    
    for(i in 1:n)
    {
      if(n==1)
      {
        if((!remove.zero.vec) | sum(abs(tmpW)) > 0) 
          listW[[ length(listW)+1 ]] <- tmpW
      }else
      {
        if((!remove.zero.vec) | sum(abs(tmpW[,i])) > 0) 
          listW[[ length(listW)+1 ]] <- tmpW[,i]
      }
      
    }
    
    if(length(listW)==0) 
    {
      if(clf$params$verbose)
      {
        warning("No model found with this configuration") #updated Edi
      }
      models = NULL
      
    }else
    {
      #models = listOfDenseVecToModelCollection(clf = clf, X = X, y = y, v = listW) #updated Edi
      models = listW #updated Edi
    }
  }
  
  return(models)
}


# Solve over the weights
#
# This method is to test PredOmics without Randomized Rounding method.
# @param x matrix or dataframe of predictors, of dimension n*p; each row is an observation vector.
# @param y response variable (1 or -1)
# @param nfolds k-fold cross-validation. Default value is 1.
# @param gamma is the hinge loss parameter.. Defines the margin
# @param lb is the lower bound of coefficients
# @param ub is the upper bound of coefficients
# @param k.sparse is the sparsity ( non-negative real value). Default value is \code{k.sparse = NULL} - no constraint.
# @param vartype is the type of coefficients : \code{"integer", "binary", "real"}. Default \code{vartype = "integer"}
# @param type "AUC" or "error"
# @param plotROC TRUE or FALSE 
# @return Return a list of objects.
norr <- function(x, y, nfolds = 1, k.sparse = NULL, gamma = 1.0, vartype = "integer", lb = -1.0, ub = 1.0, type = "AUC"){
  check.input(x = x, y = y, nfolds = nfolds, k.sparse = k.sparse, gamma = gamma,  vartype = vartype)
  set.seed(42202210) # Different "seed" give us different result.
  temps <- proc.time()
  this.call = match.call()
  print("Running ... ")
  L = list()
  # In this cas: 1-fold (nfolds == 1)
  if (nfolds == 1) {
    l <- buildlp(x, y, k.sparse = k.sparse, vartype = vartype, gamma = gamma, lb = lb, ub = ub)
    wb <- runsolver(l, x)
    w_no_zero = length(which(wb[-length(wb)] != 0, arr.ind = TRUE)) # DF 
    
    # confusion matrix
    # mconf =  matrix.conf(x = x, y = y, wb =  wb, vIdxCM = c(length = length(y)))
    
    if (type == "error") {
      # Calcul l'erreur
      best_score = round(error.rate(x, y, wb),4)
      best_model = colnames(x)[which(wb[-length(wb)] != 0, arr.ind = TRUE)]
      result = list(Call = this.call, Df = w_no_zero, best_model = best_model, coeff_best_model = wb[-length(wb)], best_score = best_score, Time = round(((proc.time() - temps)[3][[1]]),4))
    } else if (type == "AUC") {
      # calculer l'AUC 
      if ( max(abs(wb[-length(wb)])) != 0 ) {
        best_score = terDA.AUC(wb, x, y)
        best_model = colnames(x)[which(wb[-length(wb)] != 0, arr.ind = TRUE)]
        # LR = terDA.LR(x, wb, vartype)
        result = list(Call = this.call, Df = w_no_zero, best_model = best_model, coeff_best_model = wb[-length(wb)], best_score = best_score, Time = round(((proc.time() - temps)[3][[1]]),4))
        print("DONE.")
      } else {
        print("Can't solve this case. All of coefficients = 0 !")
        result = list(Call = this.call, Df = NA, best_model = best_model, coeff_best_model = wb[-length(wb)], best_score = best_score, Time = round(((proc.time() - temps)[3][[1]]),4))
      }
    }
  } else if (nfolds > 1) {
    print("Updating ... ")
    # moved to file terDA_unused.R in section UPDATING
  } else {print("Parameters are wrong!")}
  return(result)
}
# END
# 
# 


# Round and Solve over all weights : RRallW()
# @description This method is to test PredOmics with Randomized Rounding.
# @param x matrix or dataframe of predictors, of dimension n*p; each row is an observation vector.
# @param y response variable (1 or -1)
# @param nfolds k-folds cross-validation (nfolds = 1 , 2, 3, etc.). Default value is 1.
# @param k.sparse is the sparsity (non-negative real value). Default value is \code{k.sparse = NULL} - no constraint.
# @param gamma is the hinge loss parameter. Defines the margin
# @param lb is the lower bound of coefficients.
# @param ub is the upper bound of coefficients.
# @param nRR Number of iterations using Randomized Rounding method. Default value is 20.
# @param vartype is the type of coefficients : \code{"integer", "binary", "real"}. Default \code{vartype = "real"}
# @return A list of objects.
RRallW <- function(x, y, nfolds = 1, k.sparse = NULL, nRR = 20, gamma = 1.0, vartype = "real", type = "AUC", lb = -1.0, ub = 1.0) {
  check.input(x = x, y = y, nfolds = nfolds, gamma = gamma, nRR = nRR, k.sparse = k.sparse, vartype = vartype)
  this.call = match.call()
  set.seed(42202210) # Different "seed" give us different result.
  print("Running ...")
  temps <- proc.time()
  dim_x2 = dim(x)[2] # nombre des attributes
  if (nfolds == 1) {
    best_score = c() ; AUC.tmp = c()
    wb.tmp = list()
    
    # Solver lp
    l <- buildlp(x, y, k.sparse = k.sparse, vartype = vartype, gamma = gamma, lb = lb, ub = ub)
    wb <- runsolver(l, x)
    
    # Index et valeurs correspondant de tous les coefs soient deja entiers
    fidx = which(is.wholenumber(wb[-length(wb)]))
    fval = round(wb[-length(wb)][fidx],0)
    
    # index et valeurs coresspondant de tous les coefs restants ( pas encore entiers)
    if (length(fidx) > 0) {
      fidx_res = c(1:dim_x2)[-fidx]
    } else {fidx_res = c(1:dim_x2)}

    if (length(fidx_res) != 0) { # It means that : NOT all of coefficients already are integers.
      fval_res = wb[fidx_res]
      for (m in 1:nRR) {
        fval_res_rr = terDA.rr(fval_res) # Applying method RR on all of fval_res
        wb[c(fidx,fidx_res)] = c(fval, fval_res_rr) # Deja tester beaucoup de fois. C'est toujours vrai.
        wb.tmp[[m]] = wb
      
        if (type == "error") {
          best_score[m] = error.rate(x, y, wb) 
        } else if (type == "AUC") {
          if ( max(abs(wb[-length(wb)])) != 0 ) {
            AUC.tmp[m] = terDA.AUC(wb, x, y)
          } else { 
            print("Suggest: Increase nRR.")
            AUC.tmp[m] = NA
          }
        }
      }
      
      # type : AUC or error 
      if (type == "error") {
        best_score.extremum = min(best_score, na.rm = TRUE)
        idx.best_score = which(best_score == best_score.extremum,  arr.ind = TRUE)[1] 
      } else if (type == "AUC") {
        best_score.extremum = max(AUC.tmp, na.rm = TRUE)
        idx.best_score = which(AUC.tmp == best_score.extremum,  arr.ind = TRUE)[1] 
      }
      wb.res =  wb.tmp[[idx.best_score]]
      coeff_best_model = wb.res[-length(wb.res)]
      w_no_zero = length(which(coeff_best_model != 0, arr.ind = TRUE)) 
      best_model = colnames(x)[which(coeff_best_model != 0, arr.ind = TRUE)]
      result = list(Call = this.call, Df = w_no_zero, best_model = best_model, coeff_best_model = coeff_best_model, best_score = best_score.extremum, Time = round(((proc.time() - temps)[3][[1]]),4))
    } else {
      if ( type == "error") {
        best_score = error.rate(x, y, wb) 
      } else if (type == "AUC") {
        best_score = terDA.AUC(wb, x, y) 
      }
      w_no_zero =  length(which(wb[-length(wb)] != 0, arr.ind = TRUE)) # DF
      coeff_best_model = wb[-length(wb)]
      best_model = colnames(x)[which(coeff_best_model != 0, arr.ind = TRUE)]
      result = list(Call = this.call, Df = w_no_zero, best_model = best_model, coeff_best_model = wb[-length(wb)], best_score = best_score, Time = round(((proc.time() - temps)[3][[1]]),4))
    }
  } else if (nfolds > 1) {
    print("Updating ... ")
  }
  return(result)
}

#   # moved to file terDA_unused.R in section ELSE




############################## RRkBest() ###########################################################################
##### This function is the best ! 
##### But it's hard to implement correctly, because we need recall the solver many time, but there is k.sparse !!! 
##### We can't sure that : sum|weights| <= k.sparse when we use Randomized Rounding. Ignore these cases.
##### Ca perde beaucoup beaucoup de temps !!!!! 
####################################################################################################################
# Repeat Round and Solve over Best k weights
#
# @description This method is to test PredOmics with forced Randomized Rounding adding the best elements selection.
# @param x matrix or dataframe of predictors, of dimension n*p; each row is an observation vector.
# @param y response variable (1 or -1)
# @param nRR Number of iterations using Randomized Rounding method. Default value is 20.
# @param kBest Number of best candidat for Randomized Rounding by step. Default value is 5.
# @param nfolds k-fold cross-validation. Default value is 1.
# @param lb is the lower bound of coefficients
# @param ub is the upper bound of coefficients
# @param k.sparse is the sparsity (non-negative real value). Default value is \code{k.sparse = NULL} - no constraint.
# @param gamma is the hinge loss parameter.. Defines the margin.
# @param vartype is the type of coefficients : \code{"integer", "binary", "real"}. Default \code{vartype = "real"}
# @return A list of objects.
RRkBest <- function(x, y, nfolds = 1, k.sparse = NULL, gamma = 1.0, nRR = 20, kBest = 5, vartype = "real",  lb = -1.0, ub = 1.0, type = "AUC") {
  check.input(x = x, y = y, nfolds = nfolds, gamma = gamma, nRR = nRR, k.sparse = k.sparse, kBest = kBest, vartype = vartype)
  this.call = match.call()
  #set.seed(42202210)
  temps <- proc.time()
  dim_x2 = dim(x)[2] # Nombre des colonnes
  print("Running ...")
  
  ######################### FIRST CAS: nfolds  = 1 #########################################
  if (nfolds == 1) {
    l <- buildlp(x, y, k.sparse = k.sparse, vartype = vartype, gamma = gamma, lb = lb, ub = ub)
    wb <- runsolver(l, x)
    
    wb_init = wb # wb initiation
    # index et valeurs correspondants sont deja INTEGER
    fidx_init = which(is.wholenumber(wb[-length(wb)]))
    fval_init = round(wb[-length(wb)][fidx_init],0)
    
    # index et valeurs correspondants  ne sont pas encore INTERGER.
    if (length(fidx_init) > 0) {
      fidx_res_init = c(1:dim_x2)[-fidx_init]
    } else {fidx_res_init = c(1:dim_x2)}
    fval_res_init = wb[fidx_res_init]
    
    if (length(fidx_res_init) > 0) {
      # longueur des coefs restants ( ne pas encore INTEGER)
      len_fidx_res_init = length(fidx_res_init)
      count = 1 ; v_err = c(); listRR = list() ; AUC.tmp = c()
      
      # Demarrer des iterarions
      for (ii in 1:nRR) {
        # Initiation
        fidx = fidx_init
        fval = fval_init
        fidx_res = fidx_res_init
        fval_res = fval_res_init
        len_fidx_res = len_fidx_res_init
        
        for (m in 1:len_fidx_res) {
          # Recalcul fidx,fval,fidx_res,fval_res avec nouveau 'wb'
          fidx = which(is.wholenumber(wb[-length(wb)]))
          fval = round(wb[-length(wb)][fidx],0)
          fidx_res = c(1:dim_x2)[-fidx]
          
          if (length(fidx_res) == 0 ) {
            break
          } else {
            fval_res = wb[fidx_res]
            # Selection un meilleur candidat: 'idxRR'
            fidx_min = terDA.mindist(fval_res, u = kBest)
            idxRR = fidx_res[fidx_min]
            # Randrounding sur les meilleurs candidats
            valRR_current = terDA.rr(fval_res[fidx_min])
            
            # Nouveau 'fidx' et nouveau 'fval'
            fidx = c(fidx,idxRR)
            fval = c(fval,valRR_current)
            
            # Executation buildlp() avec nouveau 'w'
            l <- buildlp(x,y, vartype = vartype, k.sparse = k.sparse, gamma = gamma, fcoefidx = fidx, fcoefval = fval, lb = lb , ub = ub)
            wb <- runsolver(l, x)
          }
        }
        wb[-length(wb)] = round(wb[-length(wb)],10)
        listRR[[count]] = wb
        
        if (length(wb) > 1) { 
          if (type == "error") {
            v_err[count] = error.rate(x,y,wb)
          } else if (type == "AUC") {
            if ( max(abs(wb[-length(wb)])) != 0 ) {
              AUC.tmp[count] = terDA.AUC(wb, x, y)
            } else { 
              print("Suggest: Increase nRR.")
              AUC.tmp[count] = NA
            }
          }
        } else {
          if (type == "error") {
            v_err[count] = NA 
          } else {AUC.tmp[count] = NA}
        }
        count = count  + 1
        wb = wb_init
      } # finis un nRR. Il y a encore (nRR - 1) iterations
      
      # type : AUC or error 
      if (type == "error") {
        best_score.extremum = min(v_err, na.rm = TRUE)
        idx.best_score = which(v_err == best_score.extremum,  arr.ind = TRUE)[1] 
      } else if (type == "AUC") {
        best_score.extremum = max(AUC.tmp, na.rm = TRUE)
        idx.best_score = which(AUC.tmp == best_score.extremum,  arr.ind = TRUE)[1] 
      }
      wb.res =  listRR[[idx.best_score]]
      coeff_best_model = wb.res[-length(wb.res)]
      w_no_zero = length(which(coeff_best_model != 0, arr.ind = TRUE)) 
      best_model = colnames(x)[which(coeff_best_model != 0, arr.ind = TRUE)]
      result = list(Call = this.call, Df = w_no_zero, best_model = best_model, coeff_best_model = coeff_best_model, best_score = best_score.extremum, Time = round(((proc.time() - temps)[3][[1]]),4))
      
      #       # Minimal error
      #       which_min_err = which.min(v_err)
      #       wb = listRR[[which_min_err]]
      #       emp_error = min(v_err)
      #       
      #       # AUC
      #       auc = terDA.AUC(wb, x , y)
      #       mconf = matrix.conf(x = x, y = y, wb =  wb, vIdxCM = c(length = length(y)))
    } else {
      if ( type == "error") {
        best_score = error.rate(x, y, wb) 
      } else if (type == "AUC") {
        best_score = terDA.AUC(wb, x, y) 
      }
      w_no_zero =  length(which(wb[-length(wb)] != 0, arr.ind = TRUE)) # DF
      coeff_best_model = wb[-length(wb)]
      best_model = colnames(x)[which(coeff_best_model != 0, arr.ind = TRUE)]
      result = list(Call = this.call, Df = w_no_zero, best_model = best_model, coeff_best_model = wb[-length(wb)], best_score = best_score, Time = round(((proc.time() - temps)[3][[1]]),4))
      
      
      #       emp_error = error.rate(x, y, wb)
      #       auc = terDA.AUC(wb , x , y)
      #       print(objf)
      #       mconf = matrix.conf(x = x, y = y, wb =  wb, vIdxCM = c(length = length(y)))
      #       # result = list(Call = this.call, Error = round(100*emp_error,2), Coef = wb, Time = round(((proc.time() - temps)[3][[1]]),4))
    }
    # LR = terDA.LR(x ,wb)
    # w_no_zero = which(wb[-length(wb)] != 0, arr.ind = TRUE)
    # result = list(Call = this.call, matrix.conf = mconf, fobj = objf, Df = length(w_no_zero), Coef = wb, AUC = auc, Error = round(100*emp_error,2), LR = LR, Time = round(((proc.time() - temps)[3][[1]]),4))
    # result = list(Call = this.call, matrix.conf = mconf, Coef = wb[-length(wb)], LR = LR, AUC = auc, Error = round(100*emp_error,2), Time = round(((proc.time() - temps)[3][[1]]),4))
    print("DONE.")
  }  else if (nfolds > 1) {
    print("Updating ... ")
    # moved to file terDA_unused_R in RRkbest_else section 
  } else {print("Parameters are wrong!")}
  return(result)
}
############# End of the function terDA.kbest  ################





################### RRkRand ################################
# Repeat Round and Solve over k random weights
#
# @description This method is to test PreOmics with forced Randomized Rounding.
# @param x matrix or dataframe of predictors, of dimension n*p; each row is an observation vector.
# @param y response variable (1 or -1)
# @param nRR Number of iterations using Randomized Rounding method. Default value is 20.
# @param kStep Step of Rounding. Default value is 5.
# @param k.sparse is the sparsity (non-negative real value). Default value is \code{k.sparse = NULL} - no constraint.
# @param gamma is the hinge loss parameter.. Defines the margin
# @param lb is the lower bound of coefficients
# @param ub is the upper bound of coefficients
# @param nfolds k-fold cross-validation. Default value is 1. (nfolds = 10 correspondent the case 10-fold cross-validations)
# @param vartype is the type of coefficients : \code{"integer", "binary", "real"}. Default \code{vartype = "real"}
# @return A list of objects.
RRkRand <- function(x ,y , nfolds = 1, k.sparse = NULL, gamma = 1.0, nRR = 20, kStep = 15, vartype = "real", lb = -1.0, ub = 1.0, type = "AUC"){
  check.input(x = x, y = y, nfolds = nfolds, gamma = gamma, nRR = nRR, k.sparse = k.sparse, kStep = kStep, vartype = vartype)
  this.call = match.call()
  #set.seed(42202210)
  temps <- proc.time()
  dim_x2 = dim(x)[2] # nombre de colonnes
  print("Running ... ")
  
  if (nfolds == 1) {
    # Lancer buildlp()  pour avoir 'wb' initiation
    l <- buildlp(x, y, vartype = vartype, k.sparse = k.sparse, gamma = gamma, lb = lb, ub = ub)
    wb <- runsolver(l, x)
    wb_init = wb
    
    # Index (fidx) et valeurs (fval) correspondant de tous les coefs soient deja entiers (already interger)
    fidx = which(is.wholenumber(wb[-length(wb)]))
    fval = round(wb[-length(wb)][fidx],0)
    
    # fidx_init et fval_init sont de revenir des iterations quand nRR change apres.
    # (C-a-d que Sauvgarder fidx_init et fval_init quand nRR change)
    fidx_init = fidx
    fval_init = fval
    
    # index et valeurs coresspondant de tous les coefs restants ( pas encore entiers)
    if (length(fidx) > 0) {
      fidx_res = c(1:dim_x2)[-fidx]
    } else {fidx_res = c(1:dim_x2)}
    
    if (length(fidx_res) == 0) { # it means that ALL of coefficients already are INTEGER.
      
      if ( type == "error") {
        best_score = error.rate(x, y, wb_init) 
      } else if (type == "AUC") {
        best_score = terDA.AUC(wb_init, x, y) 
      }
      w_no_zero =  length(which(wb_init[-length(wb_init)] != 0, arr.ind = TRUE)) # DF
      coeff_best_model = wb_init[-length(wb_init)]
      best_model = colnames(x)[which(coeff_best_model != 0, arr.ind = TRUE)]
      result = list(Call = this.call, Df = w_no_zero, best_model = best_model, coeff_best_model = wb[-length(wb)], best_score = best_score, Time = round(((proc.time() - temps)[3][[1]]),4))
      #       
      #       
      #       emp_ERROR = error.rate(x, y, wb_init)
      #       # AUC
      #       AUC = terDA.AUC(wb_init, x, y)
      #       # matrix confusion
      #       mconf = matrix.conf(x = x, y = y, wb =  wb_init, vIdxCM = c(length = length(y)))
      #       w_no_zero = which(wb_init[-length(wb_init)] != 0, arr.ind = TRUE)
      #       # result = list(Call = this.call, fobj = fobj.init, matrix.conf = mconf, AUC = AUC, Error = 100*emp_ERROR, Df = length(w_no_zero), Coef = wb_init, Time = round(((proc.time() - temps)[3][[1]]),4))
      #       result = list(Call = this.call, matrix.conf = mconf, AUC = AUC, Error = 100*emp_ERROR, Coef = wb_init[-length(wb_init)], Time = round(((proc.time() - temps)[3][[1]]),4))
      #       
    } else {# Starting ...
      
      fval_res = wb[fidx_res]
      len_fidx_res = length(fidx_res) # length
      
      #Start forced RR
      count = 1 ;        v_err = c() ; AUC.tmp = c();    listRR = list() # listRR contenant les listes des coefs
      
      ##############################################################################
      # THERE ARE 'nRR' ITERATIONS
      for (ii in 1:nRR) {
        # Prendre une permutation of 'len_fidx_res' composants
        sampleidx <- sample(fidx_res, len_fidx_res, replace = F)
        for (m in 1:len_fidx_res)
        {
          if (m %% kStep == 0) {
            rr_tmp_divised = terDA.rr(wb[sampleidx[m]])
            fidx = c(fidx, sampleidx[m])
            fval = c(fval, rr_tmp_divised)
            l <- buildlp(x, y, k.sparse = k.sparse, gamma = gamma, vartype = vartype, fcoefidx = fidx, fcoefval = fval, lb = lb, ub = ub)
            wb <- runsolver(l, x)
          } else {
            # Apply Randomized rounding
            rr_tmp = terDA.rr(wb[sampleidx[m]])
            fidx = c(fidx, sampleidx[m])
            fval = c(fval, rr_tmp)
            wb[fidx] = fval
          }
        }
        
        wb[-length(wb)] = round(wb[-length(wb)],10)
        listRR[[count]] = wb
        if (length(wb) > 1) { 
          if (type == "error") {
            v_err[count] = error.rate(x,y,wb)
          } else if (type == "AUC") {
            if ( max(abs(wb[-length(wb)])) != 0 ) {
              AUC.tmp[count] = terDA.AUC(wb, x, y)
            } else { 
              print("Suggest: Increase nRR.")
              AUC.tmp[count] = NA
            }
          }
        } else {
          if (type == "error") {
            v_err[count] = NA 
          } else {AUC.tmp[count] = NA}
        }
        #
        count = count + 1
        fidx = fidx_init
        fval = fval_init
        wb = wb_init
      } # finid un nRR
      
      # type : AUC or error 
      if (type == "error") {
        best_score.extremum = min(v_err, na.rm = TRUE)
        idx.best_score = which(v_err == best_score.extremum,  arr.ind = TRUE)[1] 
      } else if (type == "AUC") {
        best_score.extremum = max(AUC.tmp, na.rm = TRUE)
        idx.best_score = which(AUC.tmp == best_score.extremum,  arr.ind = TRUE)[1] 
      }
      wb.res =  listRR[[idx.best_score]]
      coeff_best_model = wb.res[-length(wb.res)]
      w_no_zero = length(which(coeff_best_model != 0, arr.ind = TRUE)) 
      best_model = colnames(x)[which(coeff_best_model != 0, arr.ind = TRUE)]
      result = list(Call = this.call, Df = w_no_zero, best_model = best_model, coeff_best_model = coeff_best_model, best_score = best_score.extremum, Time = round(((proc.time() - temps)[3][[1]]),4))
      
      #       # Minimal error
      #       which_min_err = which.min(v_err)
      #       wb = listRR[[which_min_err]]
      #       LR = terDA.LR(x, wb)
      #       mconf = matrix.conf(x = x, y = y, wb =  wb, vIdxCM = c(length = length(y)))
      #       w_no_zero = which(wb[-length(wb)] != 0, arr.ind = TRUE)
      #       auc = terDA.AUC(wb, x , y)
      #       # result = list(Call = this.call, fobj = v.fobj[which_min_err], Df = w_no_zero, matrix.conf = mconf, Error = round(min(v_err)*100,2), AUC = auc, LR = LR, Coef = listRR[[which_min_err]], Time = round(((proc.time() - temps)[3][[1]]),4))
      #       # DERNIERE CHANCE
      #       result = list(Call = this.call, matrix.conf = mconf, Error = round(min(v_err)*100,2), AUC = auc, Coef = listRR[[which_min_err]][-length(listRR[[which_min_err]])], LR = LR, Time = round(((proc.time() - temps)[3][[1]]),4))
      #       
    }
    print("DONE.")
  } else if (nfolds > 1) {
    # How to do in this case:
    # Compute min(error) between 100 errors for found a best model (It's a vector of coefficients)
    # From this vector of coefficients, we can compute the error with database test.
    
    count = 1 ;    v_err = c() ;    L = list() ;    v_err_fold = c();     v.AUC = c() ; v.fobj = c(); v.fobj.f = c()
    lfolds = create.folds(y, k = nfolds, list = TRUE, returnTrain = FALSE) # Index of 10-fold
    
    for (jj in 1:nfolds) # Normalement: nfolds == 10
    {
      # preparation le dataset
      x_train = x[-lfolds[[jj]],] ;      row.names(x_train) = NULL
      y_train = y[-lfolds[[jj]]]  ;      row.names(y_train) = NULL
      x_test = x[lfolds[[jj]],]   ;      row.names(x_test) = NULL
      y_test = y[lfolds[[jj]]]    ;      row.names(y_test) = NULL
    
      # Executer le modele
      l <- buildlp(x_train, y_train, k.sparse = k.sparse, gamma = gamma, vartype = vartype, lb = lb , ub = ub)
      wb <- runsolver(l, x_train)
      
      # Index et valeurs correspondant de tous les coefs soient deja entiers
      fidx = which(is.wholenumber(wb[-length(wb)]))
      fval = round(wb[-length(wb)][fidx], 0)
      
      # fidx_init et fval_init sont de revenir des itetations quand nRR change
      fidx_init = fidx
      fval_init = fval
      
      # index et valeurs coresspondant de tous les coefs restants ( pas encore entiers)
      if (length(fidx) > 0) {
        fidx_res = c(1:dim_x2)[-fidx]
      } else {fidx_res = c(1:dim_x2)}
      
      if (length(fidx_res) == 0) { # Tat ca cac he so da la so nguyen roi.
        matrix.conf(x = x_test, y = y_test, wb =  wb, vIdxCM = c(length = length(y_test)))
        v_err_fold[jj] = error.rate(x_test, y_test, wb)
        v.AUC[jj] = terDA.AUC(wb, x_test, y_test)
        print(paste("k = ",jj, "Emp error =", round(v_err_fold[jj]*100, 4),"%",",AUC =", v.AUC[jj]))
      } else {
        fval_res = wb[fidx_res]
        len_fidx_res = length(fidx_res) # length
        
        # nRR iterations
        for (ii in 1:nRR) {
          # Prendre une permutation of 'len_fidx_res' composant
          sampleidx <- sample(fidx_res, len_fidx_res, replace = FALSE)
          for (m in 1:len_fidx_res)
          {
            ##### Executation buildlp() avec nouveau wb:
            if (m %% kStep == 0) {
              rr_tmp_divised10 = terDA.rr(wb[sampleidx[m]])
              
              fidx = c(fidx, sampleidx[m])
              fval = c(fval, rr_tmp_divised10)
              l <- buildlp(x_train, y_train, vartype = vartype, fcoefidx = fidx, fcoefval = fval, gamma = gamma, lb = lb, ub = ub)
              wb <- runsolver(l, x_train)
            } else {
              rr_tmp10 = terDA.rr(wb[sampleidx[m]])
              fidx = c(fidx, sampleidx[m])
              fval = c(fval, rr_tmp10)
              wb[fidx] = fval
            }
          } # finis une boucle. Maitenant, tous les coeffs sont deja entiers
          
          L[[count]] = wb
          v_err[count] = error.rate(x_train, y_train, wb) # Ca nous donne une valeur d'erreur chaque cas. 'count' augmente.
          count = count + 1
          fidx = fidx_init
          fval = fval_init
          wb = wb_init10
        } # fini un nRR. Il reste encore (nRR - 1) iterations.
        
        L_min = L[[which.min(v_err)]]
        v.fobj.f[jj] = v.fobj[which.min(v_err)]
        # confusion matrix
        mconf = matrix.conf(x = x_test, y = y_test, wb =  L_min, vIdxCM = c(length = length(y_test)))
        v_err_fold[jj] = error.rate(x_test, y_test, L_min)
        v.AUC[jj] = terDA.AUC(L_min, x_test, y_test)
        print(paste("k = ",jj, "Emp error =", round(v_err_fold[jj]*100, 4),"%",",AUC =", v.AUC[jj]))
      }
      
      v_err = c() # re-initiation v_err
      count = 1   # re-initiation count
      L = list()  # re-initiation L
    } # fini '10-fold'
    
    # v_err_fold <- v_err_fold[!is.na(v_err_fold) & !is.null(v_err_fold)]
    #         result = list(Call = this.call, Error = paste(round((100*sum(v_err_fold)/ length(v_err_fold)), 4), "\u00B1", round(sd(v_err_fold), 4)), Time = round(((proc.time() - temps)[3][[1]]),4))
    # result = list(Call = this.call, fobj = sum(v.fobj.f)/length(v.fobj.f), AUC = round(sum(v.AUC)/length(v.AUC),4), Error = round((100*sum(v_err_fold)/ length(v_err_fold)), 2), Time = round(((proc.time() - temps)[3][[1]]),4))
    # DERNIERE CHANCE
    result = list(Call = this.call,  AUC = round(sum(v.AUC)/length(v.AUC),4), Error = round((100*sum(v_err_fold)/ length(v_err_fold)), 2), Time = round(((proc.time() - temps)[3][[1]]),4))
    
    print("DONE.")
  } else {print("Parameters are wrong!")}
  return(result)
}
###################################### End of the function RRkRand ####################








buildAndSolveReducedLp <- function(x, y, vartype = "integer", k.sparse = NULL, lb = -1.0, ub = 1.0, gamma = 1.0, mytimeout = 120, depthlim = 50, fcoefidx = NULL, fcoefval = NULL) {
  
  if(!is.null(fcoefidx)){
    x.reduced <- x[,-fcoefidx]
    #lhsoffset <- (as.matrix(x[,fcoefidx]) %*% as.matrix(fcoefval))*y
    lhsoffset <- rep(0,nrow(x))
    lr <- simplebuildlp(x.reduced, y=y, vartype = vartype, k.sparse = k.sparse-length(fcoefval!=0), lb = lb, ub = ub, gamma = gamma, mytimeout = mytimeout, depthlim = depthlim, lhsoffset = lhsoffset)
  }else{
    x.reduced <- x
    lhsoffset <- rep(0,nrow(x))
    lr <- simplebuildlp(x.reduced, y=y, vartype = vartype, k.sparse = k.sparse, lb = lb, ub = ub, gamma = gamma, mytimeout = mytimeout, depthlim = depthlim, lhsoffset = lhsoffset)
  }
  
  wbr <- runsolver(lr, x.reduced) #reduced wb
  return(wbr)
}



# Solve and returns a vector (w1,..,wn,b)
runsolver <- function(l,x) {
  m <- nrow(x)
  n <- ncol(x)
  r <- solve(l)
  #0: "optimal solution found"
  #1: "the model is sub-optimal"
  #2: "the model is infeasible"
  #3: "the model is unbounded"
  #4: "the model is degenerate"
  #5: "numerical failure encountered"
  #6: "process aborted"
  #7: "timeout"
  #9: "the model was solved by presolve"
  #10: "the branch and bound routine failed"
  #11: "the branch and bound was stopped because of a break-at-first or break-at-value"
  #12: "a feasible branch and bound solution was found"
  #13: "no feasible branch and bound solution was found"
  if (r >= 2) {
    stop(paste("runsolver exited with error code:",r))
    #return(NULL)
  } else {
    #     if (r == 0) {
    #       print("Optimal solution found!")
    #     } else if (r == 1) {
    #       print("The model is sub-optimal.")
    #     }
    # print(r)
    # print(get.objective(l))
    return( c(get.variables(l)[1:n], get.variables(l)[n + m + 1]))
  }
}







# check.input <- function(x, y, nfolds =1, nRR = 2, k.sparse = c(), gamma = 1, kBest = 10^10, kStep = 10^10, vartype = cplexAPI::CPX_CONTINUOUS) {
#   #
#   v_vartype = c(cplexAPI::CPX_CONTINUOUS,cplexAPI::CPX_INTEGER, cplexAPI::CPX_BINARY)
check.input <- function(x, y, nfolds =1, nRR = 2, k.sparse = c(), gamma = 1, kBest = 10 ^ 10, kStep = 10 ^ 10, vartype = "real") {
  v_vartype = c("real", "integer", "binary")
  if (!is.element(vartype, v_vartype)) {
    stop("vartype : real or binary or integer")
  }
  #
  if ( !is.data.frame(x)) {
    stop("'x' - It's not a dataframe or a matrix.")
  }
  # Check parameter: fold
  if ( !is.wholenumber(nfolds) || (nfolds < 1)) {
    stop("'nfolds' must be a positive integer.")
  }
  # Check nRR
  if ( !is.wholenumber(nRR) || (nRR < 1) ) {
    stop("'nRR' must be a positive integer. Stopped!")
  }
  # First column must be the class (here, there are 2 class.)
  if (length(unique(y)) != 2)
  {
    stop("Just 2 class in the dataset.")
  }
  # Check: value of class: must be -1 or 1?
  if ( length(y[!(y %in% c(-1,1))]) != 0 )
  {
    stop("The class must be -1 or 1.")
  }
  # Check k.sparse
  if (!is.null(k.sparse) &&  (k.sparse < 0 ))
  {
    stop("k.sparse must be a positive or NULL . Stopped!")
  }
  #
  if (is.vector(k.sparse) && (length(k.sparse[which(k.sparse <= 0)]) > 0 )  )
  {
    stop("k.sparse must be a positive. Stopped!")
  }
  #
  if (!is.vector(gamma))
  {
    stop("gamma must be a positive. Stopped!")
  }
  #
  if (is.vector(gamma) && (length(gamma[which(gamma <= 0)]) > 0 )  )
  {
    stop("gamma must be a positive. Stopped!")
  }
  
  # Check kBest
  if ( !is.wholenumber(kBest) || (kBest <= 0) ) {
    stop("'kBest' must be a positive integer. Stopped!")
  }
  # Check kStep
  if ( !is.wholenumber(kStep) || (kStep <= 0) ) {
    stop("'kStep' must be a positive integer. Stopped!")
  }
}
######### End of the function checkInput #############


# Fonction classification
classify <- function(ex,w,b) {
  z <- t(w) %*% t(ex)
  if (z + b >= 0) {
    return(1)
  } else {
    return(-1)
  }
}
######### End of the function classify #############









# Check a value is integer or not?
# http://stackoverflow.com/questions/3476782/how-to-check-if-the-number-is-integer
is.wholenumber <-  function(x, tol = .Machine$double.eps ^ 0.3)  abs(x - round(x)) < tol
# tmp = c(21, 2.3, -2.2, -1, 0, 0.01)
# is.wholenumber(tmp)
# [1]  TRUE FALSE FALSE  TRUE  TRUE FALSE




# Compute the error
error.rate <- function(x,y,wb) {
  w <- wb[-length(wb)]
  b <- wb[length(wb)]
  etot <- 0.0
  for (i in 1:nrow(x)) {
    ex <- x[i,]
    etot <- etot + abs( classify(ex,w,b) - y[i])/2
  }
  return(etot / nrow(x))
}
######### End of the function  error.rate #############




#######################################################################
#' Find the number of weights not yet integer.
#'
#' This method return a maximum number of weights of the model not yet integer.
#' @param x matrix or dataframe of predictors, of dimension n*p; each row is an observation vector.
#' @param y response variable (1 or -1)
#' @param nfolds k-folds cross-validation that we want test. Dedault value is 1.
#' @return An integer is number of coefficients not yet integer.
#' @param lb is the lower bound of coefficients
#' @param ub is the upper bound of coefficients
#' @param k.sparse is the sparsity (non-negative real value). Default value is \code{k.sparse = NULL} - no constraint.
#' @param gamma is the hinge loss parameter.. Defines the margin
#' @param vartype is the type of coefficients : \code{cplexAPI::CPX_INTEGER, cplexAPI::CPX_BINARY, cplexAPI::CPX_CONTINUOUS}. Default \code{vartype = cplexAPI::CPX_INTEGER}
#' @examples
#' library(PredOmics)
#' findk(DATAMETA1[, -1], DATAMETA1[, 1], nfolds = 1)
#' findk <- function(x, y, nfolds = 1, gamma = 1, k.sparse = NULL, vartype =  cplexAPI::CPX_CONTINUOUS, lb = -1.0, ub = 1.0)
findk <- function(x, y, nfolds = 1, gamma = 1, k.sparse = NULL, vartype =  "real", lb = -1.0, ub = 1.0)
{
  check.input(x = x, y = y, gamma = gamma, k.sparse = k.sparse, vartype = vartype)
  #set.seed(42202210)
  this.call = match.call()
  
  # Preparation les donnees
  dim_x2 = dim(x)[2]
  v_k = c()
  count = 1
  if (nfolds == 1) {
    wb <- buildlp(x, y, vartype = vartype, gamma = gamma, k.sparse = k.sparse, lb = lb, ub = ub)
    fidx = which(is.wholenumber(wb[-length(wb)]))
    fidx_res = c(1:dim_x2)[-fidx]
    v_k[count] = length(fidx_res)
  }
  
  if (nfolds > 1)
  {
    lfolds = create.folds(y, k = nfolds, list = TRUE, returnTrain = FALSE) # Index of 10-fold
    for (jj in 1:nfolds)
    {
      # preparation le dataset
      x_train = x[-lfolds[[jj]],]
      row.names(x_train) = NULL
      y_train = y[-lfolds[[jj]]]
      x_test = x[lfolds[[jj]],]
      row.names(x_test) = NULL
      y_test = y[lfolds[[jj]]]
      row.names(y_test) = NULL
      
      # Executer le modele
      wb <- buildlp(x_train,y_train, k.sparse = k.sparse, gamma = gamma, vartype = vartype, lb = lb, ub = ub)
      fidx = which(is.wholenumber(wb[-length(wb)]))
      fidx_res = c(1:dim_x2)[-fidx]
      v_k[count] = length(fidx_res)
      count = count + 1
    }
  }
  print((max(v_k) + 1))  # ca veut dire qu'a partir de (max(v_k) + 1) , alors la valeur de l'erreur ne changera plus.
  result = list(call = this.call, numberk = (max(v_k) + 1))
  return(result)
}
######### End of the function terDA.findk #############






## Make the rule of learning from dataset et the weights 
# x: data
# wb : ouput of solver lp
terDA.LR <- function(x, wb, vartype = "integer"){
  label_res = colnames(x)[which( wb != 0, arr.ind = TRUE)]
  w_no_zero = which(wb[-length(wb)] != 0, arr.ind = TRUE)
  # Show learning rule. Only when vartype = "integer" / "binary"
  if ((vartype == "integer") || (vartype == "binary")) {
    reg = "0"
    for (i in 1:(length(label_res) - 1)) {
      reg  = paste(reg,"+",wb[w_no_zero[i]],"*",label_res[i], sep = " ")
    }
    reg_1 = substr(reg,5,nchar(reg))
    reg_2 = gsub("+ -1 * ", "-", reg_1, ignore.case = FALSE, perl = FALSE, fixed = T)
    reg_3 = gsub("-1 * ", "-", reg_2, ignore.case = FALSE, perl = FALSE, fixed = T)
    reg_4 = gsub("+ 1 * ", "+", reg_3, ignore.case = FALSE, perl = FALSE, fixed = T)
    reg_5 = gsub("1 * ", " ", reg_4, ignore.case = FALSE, perl = FALSE, fixed = T)
    reg_6 = gsub(" ", " ", reg_5, ignore.case = FALSE, perl = FALSE, fixed = T)
  } else {
    reg = "0"
    for (i in 1:(length(label_res) - 1)) {
      reg  = paste(reg,"+",round(wb[w_no_zero[i]],4),"*",label_res[i], sep = " ")
    }
    reg_1 = substr(reg,5,nchar(reg))
    reg_2 = gsub("+ -", "- ", reg_1, ignore.case = FALSE, perl = FALSE, fixed = T)
    reg_6 = reg_2
  }
  reg_7 = reg_6
  # DERNIERE CHANCE 
  #   if ((wb[length(wb)]) <= 0){
  #     reg_7 = paste(reg_6,"\u2265", round(abs(wb[length(wb)]),4))
  #   } else {
  #     reg_7 = paste(reg_6,"+", round(wb[length(wb)],4), "\u2265", toString(0))
  #   }
  # show(paste("Learning rule:"))
  # show(reg_7)
  return(reg_7)
}




#########################################################################
# Find smallest distance between 2 vectors
# Generalisation, cette fonction consomme trop de temps --> BESOIN D'AMELIORER 
#
# vector1 = c(2.1, 2.2,2.3, 2.0009, 0.002, 4, -1.8, -2.001)
# vector2 = c(-2,0,2)
# terDA.mindist(x = vector1,  targets = vector2, u = 3)
# 4 5 8
terDA.mindist <- function(x, targets=c(-1,0,1), u = 1 ) {
  r <- order(sapply(x, function(z) min(abs(z - targets))))
  sort(r[1:u])
}
######### End of the function terDA.mindist #############



# Randomized Rounding method (RR method)
terDA.rr <- function(w) 
{
  wr  <- rep(0.0,length(w))
  z   <- floor(w)
  for (j in 1:length(w)) 
  {
    #  set.seed(clf$params$seed) # Different "seed" give us different result. # want different results
    if(w[j] - z[j] < 0)
    {
      stop("Problem in terDA.rr: non positiv probability for rbinom")
    }
    wr[j] <- rbinom(1, 1, w[j] - z[j]) + z[j]
    if(is.na(wr[j]))
    {
      stop("Problem in terDA.rr: na produced by rbinom")
    }
  }
  return(wr)
}
####### End of the function terDA.rr #######



# Compute AUC for terDA 
# added the 20.9.2015
#' @import caTools
terDA.AUC <- function(wb, x, y){
  w_no_zero = which(wb[-length(wb)] != 0, arr.ind = TRUE)
  tAUC = 0
  for (i in 1:length(w_no_zero)) {
    tAUC = tAUC + wb[[w_no_zero[i]]]*x[w_no_zero[i]]
  }
  auc = caTools::colAUC(tAUC, y, plotROC = F)
  
  # require(pROC)
  # rocobj <- pROC::roc(y, tAUC[,1])
  # resa  = pROC::coords(rocobj, x = "best", input = "threshold", best.method = "youden")
  # print(paste("Optimal threshold:", round(resa[1],4), ".Specificity:", round(resa[2],4), ".Sensitivity:", round(resa[3],4)))
  return(round(auc[[1]],4))
}


##############################
# FUNCTION: confusion matrix
# Ouput: a confusion matrix
matrix.conf <- function(x, y, wb, vIdxCM = c())
{
  for (i in 1:length(y))
  {
    vIdxCM[i] = classify(x[i,], wb[-length(wb)], wb[length(wb)])
  }
  matrix_confusionT = table(data.matrix(y), data.matrix(vIdxCM))
  matrix_confusion  = matrix(data.frame(matrix_confusionT)[,3], ncol = 2)
  # print("Confusion matrix:")
  # show(matrix_confusion)
  return(matrix_confusion)
}
################# End of the function matrix.conf ###################



#  Analyze the model terDA
#
# @description This function is to explore the model terDA
# @param terda.model A object obtained from the function terda()
#  @return
# A list of elements
# analyzeTerDA <- function(data, terda.model){
#   # Best scores
#   bestAUC = apply(terda.model$AUC, 1, max)
#   plot(bestAUC,  ylab = "AUC", xlab = "Number of features", main = "rich_stat_all best models", type = "o")
#   
#   # Best model Index
#   dat = terda.model$AUC
#   coord = which(dat == max(dat), arr.ind = TRUE)
#   bestModelIdx = terda.model$Coef[[coord[1,][1]]][[coord[1,][2]]]
#   
#   # Best Model
#   bestModel = colnames(data[,-1])[which(bestModelIdx != 0, arr.ind = TRUE)]
#   return(list(list_best_score = bestAUC, best_score = max(bestAUC), best_model = bestModel, coeff_best_model = bestModelIdx))
# }

# Added 15/11/2015
# Delete elements of list : all of coefficients zero
delZero <- function(L, len) {
  count = 1
  LnoZero = list()
  for (i in 1:length(L)) {
    if ( max(abs(L[i][[1]][-len])) != 0 ) {
      LnoZero[count] = L[i]
      count = count + 1
    }
  }
  return(LnoZero)
}

# Make the output of terda like as a list of object : out$model_N$....
# added the 3.12.2015
analyze.terda <- function(res, type = "AUC") {
  resterda <- list()
  resterda$res <- res
  v.df = sort(unique(as.vector(res$Df)), decreasing = F)
  for (i in 1:length(v.df)) {
    tmp = which(res$Df == v.df[i], arr.ind = TRUE)
    yy       = as.data.frame(tmp)
    # yy       = yy[order(yy$row),]
    if (length(yy) == 1) {
      extremum = res$best_score
      best_model = res$best_model
      coeff_best_model = res$coeff_best_model 
      time = res$Time
    } else {
      extremum = extremum.type(matrix = res$best_score, type = type, idx = yy)
      tmp2 = as.data.frame(which(res$best_score == extremum, arr.ind = TRUE))
      best_model = res$best_model[[tmp2[1,1]]][[tmp2[1,2]]] 
      coeff_best_model = res$coeff_best_model[[tmp2[1,1]]][[tmp2[1,2]]] 
      time = res$time[tmp2[1,1],tmp2[1,2]] 
    }
    
    # Saving in a list
    resterda[[paste("model_N",v.df[i],sep = "")]] <- list()
    resterda[[paste("model_N",v.df[i],sep = "")]][["best_model"]] <- best_model
    resterda[[paste("model_N",v.df[i],sep = "")]][["coeff_best_model"]] <- coeff_best_model
    resterda[[paste("model_N",v.df[i],sep = "")]][["best_score"]] <- extremum
    resterda[[paste("model_N",v.df[i],sep = "")]][["time"]] <- time
  }
  return(resterda)
} 

# Find the extremum of a list of elements in function of a list idx
# added the 7.12.2015
extremum.type <- function(matrix, type, idx) {
  dim.row = dim(idx)[1]
  extremum = matrix[idx[1,]$row, idx[1,]$col]
  if ( type  == "AUC") {
    for (i in 1:dim.row) {
      if (extremum < matrix[idx[i,]$row, idx[i,]$col] ) {
        extremum = matrix[idx[i,]$row, idx[i,]$col]
      }
    }
  } else if (type == "error") {
    for (i in 1:dim.row) {
      if (extremum > matrix[idx[i,]$row, idx[i,]$col]) {
        extremum = matrix[idx[i,]$row, idx[i,]$col]
      }
    }
  }
  return(extremum)
}


################### END OF FICHIER CODE #################### 
####################### End of function terda() #####################


########################## END OF CODE ###################################################################################################################