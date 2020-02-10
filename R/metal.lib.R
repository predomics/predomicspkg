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
# @script: metal.lib                                      
# @author: Shasha Cui, Edi Prifti
# @date: Avril 2017                                                  
################################################################




#' #' Computes best model of a metal clf
#' #'
#' #' @description Get best metal model
#' #' @param X: dataset to classify
#' #' @param y: variable to predict
#' #' @param clf: an object containing the different parameters of the classifier
#' #' @param clf_res: the result of metal
#' #' @param k_penalty: penalty for k
#' #' @return A list of result of best model for each k, their importance feature of each best model, individuels wrongly classified
#' #' @export
#' getTheBestMetalModel<- function(clf, clf_res, X, k_penalty=0.01, evalToOrder="accuracy_",selected=1)
#' {
#'   if(length(clf_res)==3){
#'     clf_res<-clf_res$classifier 
#'   }
#'   pop<-modelCollectionToPopulation(clf_res$models)
#'   acc <- populationGet_X(evalToOrder)(pop)
#'   k <- populationGet_X("eval.sparsity")(pop)
#'   acc.penalty <- acc-(k*k_penalty)
#'   best.acc <- max(acc.penalty)
#'   epsilon <- sqrt(best.acc*(1-best.acc)/ncol(X))
#'   pop2 <- pop[acc.penalty>(best.acc - epsilon)]
#'   mod <- getMaxMinPrevalenceModel(pop2,X,selected=selected)
#'   return(mod)
#' }


#' Generate a metal list of clfs containing information on the generators and unificators
#'
#' @description Generate a metal list of clfs containing information on the generators and unificators
#' @param mat: a language/learner presence matrix that indicates which algorithms and which languages to explore
#' @param clf: a metal classifier object, with the parameter list.clfs that can be "NULL"
#' @return a list of clfs
#' @export
generator_metal <- function(mat, clf = NULL)
{
  
  # sanity check
  if(!all(colnames(mat) %in% c("terbeam","terda","terga2")))
  {
    stop("generator_metal: unavailable algorithm. Please check configuration matrix")
  }
  
  if(!all(rownames(mat) %in% c("bin","ter","terinter","bininter","ratio")))
  {
    stop("generator_metal: unavailable language. Please check configuration matrix")
  }
  
  # if metal clf exists
  if(is.null(clf))
  {
    clf <- metal()
    print("The metal clf is missing. Creating generator with default parameters")
    
    # create the generator list
    generator <- list()
    for(i in 1:nrow(mat)) # languages
    {
      l <- rownames(mat)[i]
      
      for(j in 1:ncol(mat)) # algorithms
      {
        if(mat[i,j] == 1)
        {
          switch(colnames(mat)[j],
                 terbeam={ 
                   clf.g <- terBeam(language = l)
                 },
                 terda={
                   clf.g <- terda(language = l)
                 },
                 terga2={
                   clf.g <- terga2(language = l, nb_generations = 1)
                 }
          )
          generator[[length(generator)+1]] <- clf.g
        } # end if config yes
      } # end for languages
    } # end for algorithms
    
    # Add the unificator with default
    generator[[length(generator)+1]] <- terga2()
    
    # add the attributes
    generator$attributes  <-  list()
    for(i in 1:(length(generator)-2))
    {
      generator$attributes[[length(generator$attributes)+1]] <- 'G'
    }
    generator$attributes[[length(generator$attributes)+1]] <- 'U'
    
    # give names to the generators
    a = unlist(generator[length(generator)]) 
    b = populationGet_X(element2get = "learner", toVec = TRUE, na.rm = TRUE)(pop = generator[-length(generator)])
    c = c(1: c(length(generator) -2),"")
    names(generator)[-length(generator)] <- paste(paste(a, c, sep=""), b, sep="_")
    
  } else
  {
    # if clf exists fine tune the parameters
    if(clf$params$unificator.method != "terga2")
    {
      stop(paste("generator_metal: the method",clf$params$unificator.method," can not be a unficator. Please provide terga2 as one."))
    }
    
    # create the generator list
    generator <- list()
    for(i in 1:nrow(mat)) # languages
    {
      l <- rownames(mat)[i]
      
      for(j in 1:ncol(mat)) # algorithms
      {
        if(mat[i,j] == 1)
        {
          switch(colnames(mat)[j],
                 terbeam={ 
                   clf.g <- terBeam(language = l, 
                                    maxNbOfModels = 1000,
                                    sparsity = clf$params$sparsity,
                                    nCores = clf$params$nCores,
                                    evalToFit = clf$params$evalToFit,
                                    seed = clf$params$seed)
                 },
                 terda={
                   clf.g <- terda(language = l, 
                                  sparsity = clf$params$sparsity,
                                  nCores = clf$params$nCores,
                                  evalToFit = clf$params$evalToFit,
                                  seed = clf$params$seed)
                 },
                 terga2={
                   clf.g <- terga2(language = l, 
                                   sparsity = clf$params$sparsity,
                                   nCores = clf$params$nCores,
                                   nb_generations = 1, # random generation without evolution
                                   evalToFit = clf$params$evalToFit,
                                   seed = clf$params$seed)
                 }
          )
          generator[[length(generator)+1]] <- clf.g
        } # end if config yes
      } # end for languages
    } # end for algorithms
    
    if(nrow(mat) == 1)
    {
      # Add the unificator using the language space
      l = rownames(mat)
    }else
    {
      # Add the unificator by default in terinter
      l = "terinter"
    }
    
    generator[[length(generator)+1]] <- terga2(sparsity = clf$params$sparsity,
                                               evolver = clf$params$unificator.evolver,
                                               language = l,
                                               nb_generations = 50,
                                               size_pop_random = 0,
                                               size_pop = 500,
                                               evalToFit = clf$params$evalToFit,
                                               seed=clf$params$seed)
    
    # add the attributes
    generator$attributes  <-  list()
    for(i in 1:(length(generator)-2))
    {
      generator$attributes[[length(generator$attributes)+1]] <- 'G'
    }
    generator$attributes[[length(generator$attributes)+1]] <- 'U'
    
    # give names to the generators
    a = unlist(generator[length(generator)]) 
    b = populationGet_X(element2get = "learner", toVec = TRUE, na.rm = TRUE)(pop = generator[-length(generator)])
    c = c(1: c(length(generator) -2),"")
    names(generator)[-length(generator)] <- paste(paste(a,c,sep=""), b, sep="_")
  }
  
  return(generator)
}



