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
# @script: terga1.lib.R                                          
# @author: Edi Prifti
# @author: Lucas Robin
# @date: August 2016                                                    
################################################################


#### SELECTION ####
# this function will select x% of the best individuals of a given population
selectElite <- function(evaluation, percentage) 
{
  perc <- ceiling(percentage*length(evaluation)/100)
  res <- order(evaluation,decreasing = T)[1:perc]
  return(res)
}


# this function will select x% of the individuals of a given population randomly
selectTournoi <- function(selection, percentage, discrete = TRUE) 
{
  if(discrete) 
  {
    perc <- percentage
  }else 
  {
    perc <- ceiling(percentage*length(selection)/100)
  }
  res <- selection[sample(x = 1:length(selection), size  =  perc, replace = FALSE)]
  return(res)
}


# this function will select x% of the best individuals and y% of randomly chosen ones
selectMixed <- function(evaluation, percentage_elite, percentage_tournoi) 
{
  perc <- ceiling((percentage_elite + percentage_tournoi) * length(evaluation)/100)
  perc_elite <- ceiling(percentage_elite * length(evaluation)/100)
  #perc_tournoi <- floor(percentage_tournoi * length(evaluation)/100)
  perc_tournoi <- perc - perc_elite
  
  res <- rep(NA, perc)
  # Elite
  elite <- selectElite(evaluation, percentage_elite)
  res[1:length(elite)] <- elite
  # Tournoi
  ind.pop <- seq(1,length(evaluation),1)
  tournoi <- selectTournoi(selection = ind.pop[!ind.pop %in% res], perc_tournoi, discrete = TRUE)
  
  offset <- length(elite) + 1
  res[offset : min((offset+length(tournoi)),length(res))]  <-  tournoi
  return(res)
}


populationDenseVec <- function(clf, size_world, best_ancestor = NULL)
{
  if(!is.na(clf$params$current_sparsity))
  {
    aimedSpar <- clf$params$current_sparsity
  } else 
  {
    aimedSpar <- min(clf$params$sparsity)
  }
  
  if(is.null(best_ancestor))
  {
    pop2generate <- 1:clf$params$size_pop
  } else 
  {
    if (length(which(best_ancestor != 0)) != (aimedSpar - 1)) 
    { 
      stop("populationDenseVec: The ancestor population does not have the right size!")
    }
    
    number_best_ancestors <- clf$params$size_pop * clf$params$perc_best_ancestor / 100
    best_ancestors <- list()
    world <- which(best_ancestor == 0)
    mininds <- min(size_world - length(which(best_ancestor != 0)), number_best_ancestors)
    missing_gene <- sample(world, mininds, replace = TRUE)
    
    best_ancestors <- foreach(i = 1:mininds) %do% 
    {
      res <- best_ancestor
      res[missing_gene[i]] <- sample(c(-1, 1), 1)
      return(res)
    }
    
    # pop <- best_ancestors
    pop2generate <- 1:(clf$params$size_pop - number_best_ancestors)
  }
    
  prob_1 <- aimedSpar/size_world
  prob <- c(prob_1/2, 1 - prob_1, prob_1/2)
  pop <- foreach(i = pop2generate) %do% 
  {
    res <- sample(-1:1, size_world, replace = TRUE, prob = prob)
    effectiveSpar <- sum(abs(res))
    if(effectiveSpar < aimedSpar)
    {
      notNull <- which(res != 0)
      if(length(notNull) > 0)
      {
        ind2Add <- sample((1:size_world)[-notNull], aimedSpar - effectiveSpar)
      }else
      {
        ind2Add <- sample((1:size_world), aimedSpar - effectiveSpar)
      }
      res[ind2Add] <- sample(c(-1,1), length(ind2Add), replace = TRUE)
    } else 
      
      if(effectiveSpar > aimedSpar) 
      {
        notNull <- which(res != 0)
        ind2Rm <- sample(notNull, effectiveSpar - aimedSpar)
        res[ind2Rm] <- 0
      }
    
    return(res)
  }
  
  if(!is.null(best_ancestors))
  {
    pop[(max(pop2generate)+1):clf$params$size_pop] <- best_ancestors
  }
  
  return(pop)
}


#' Creates a population of index models. 
#'
#' @description This function is used in terga1 and generates a list of index vectors in the variable space. These vectors can be unique or not. NB that if clf$params$unique_vars is set to TRUE it can take a long time to come out of the while loop which ensures the uniqueness of the individuals.
#' @param clf: the classifier parameter object
#' @param size_ind: The sparsity of the models. All the models of this population will have the same number of features.
#' @param size_world: The number of features from which we can choose the indices. This is needed to compute the combinatory space search.
#' @param best_ancestor: We can supply to the popolution an individual (vector with indeces) of a lower sparsity. This will ensure to seed part of the population with at least those genes. We added this feature after an observations that a local optimum of lower sparsity was lost in higher sparsities.
#' @param size_pop: the number of models to produce (default=NULL). This information is stored here clf$params$size_pop, but this parameter allows to override it.
#' @return a population of index models
#' @export
population <- function(clf, size_ind, size_world, best_ancestor = NULL, size_pop = NULL, seed = NULL) 
{
  # if we have a seeding ancestor
  if (!is.null(best_ancestor)) 
  {
    if(is.list(best_ancestor))
    {
      best_ancestor_genes <- unique(unlist(best_ancestor))
    }else
    {
      if(is.vector(best_ancestor))
      {
        best_ancestor_genes <- best_ancestor  
      }else
      {
        stop("population: the best ancestors is in an unknown format. Please provide a list or a vector!")
      }
    }
    
    # sanity checks
    if (max(best_ancestor_genes) > size_world) 
    { 
      stop("population: the best ancestors contain indexes larger then the available space!")
    }

    if (size_ind > size_world)  
    { 
      stop("population: The model can't have more elements than the total number of elements")
    }
    
    if(is.null(size_pop))
    {
      size_pop <- clf$params$size_pop
    }
    
    # if pop.last is not null than we ask for best descendance
    number_best_ancestors <- round(size_pop * clf$params$perc_best_ancestor / 100)
    
    best_ancestors <- list()
    world <- 1:size_world
    world.available <- world[!world %in% best_ancestor_genes]
    mininds <- min(length(world.available), number_best_ancestors)
    
    # Here we can face different configurations
    # a) the available best ancestor genes are equal or more than the sparsity of the current sparsity. In this case we can create "number_best_ancestors" individuals by randomly sampling x individuals without replacement + a pool of other random genes
    # b) the available best ancestor genes are less than the sparsity in which case we can use all of them + sample the missing ones
    
    if(length(best_ancestor_genes) >= size_ind)
    {
      for(i in 1:number_best_ancestors)
      {
        # First sample genes from best_ancestor pool
        if(!is.null(seed)){ set.seed(seed = seed + i) }
        genes_to_keep <- sample(x = best_ancestor_genes, size = size_ind-1, replace = FALSE)
        # Next add one gene from the world
        if(!is.null(seed)){ set.seed(seed = seed + i) }
        missing_gene <- sample(x = world.available, size = 1, replace = FALSE)
        # create the new descendent and save
        best_ancestors[[i]] <- c(genes_to_keep, missing_gene)
      }
    }else # smaller
    {
      for(i in 1:number_best_ancestors)
      {
        # First use the genes from the best_ancestor_genes set
        genes_to_keep <- best_ancestor_genes
        # Next add one gene from the world
        if(!is.null(seed)){ set.seed(seed = seed + i) }
        missing_gene <- sample(x = world.available, size = size_ind-length(best_ancestor_genes), replace = FALSE)
        # create the new descendent and save
        best_ancestors[[i]] <- c(genes_to_keep, missing_gene)
      }
    } # end size condition
  } # end best ancestor
  
  if(is.null(size_pop))
  {
    size_pop <- clf$params$size_pop
  }
  
  combinatory <- choose(size_world, size_ind) 
  if (is.finite(combinatory)) 
  {
    if(size_pop > combinatory)
    {
      print(paste("The chosen population is greater than the combinatory number of solutions. Reducing to the combinatory number",combinatory))
      size_pop <- combinatory
    }
  }else 
  {
    warning("infinite combinatory in population!!")
  }
  pop <- list()

  # if we wish to create a population composed of unique individuals
  if(clf$params$unique_vars)
  { 
    counting <- 1 # needed for the seed
    if(!is.null(seed))
    {
      set.seed(seed = seed + counting) 
    }
    pop[[1]] <- sort(sample(x = size_world, size = size_ind, replace = FALSE))
    
    # We need to create enough unique individuals to fill up the population
    while (length(pop) < size_pop) 
    {
      counting <- counting+1 # needed for the seed
      if(!is.null(seed))
      {
        set.seed(seed = seed + counting) 
      }
      # create an individual
      tmp <- sort(sample(x = size_world, size = size_ind, replace = FALSE))
      tmp.sort <- sort(tmp)
      
      present <- FALSE
      # go over the whole population and test its presence
      for (i in 1:length(pop)) 
      {
        # if found set the marker switch to on
        if (all(tmp.sort == sort(pop[[i]]))) 
        {
          present <- TRUE
        }
      } # end search for
      
      if (!present) # if not found add
      { 
        pop[[length(pop) + 1]] <- tmp
      }
    } # end while
  }else # non unique individuals, this is much faster
  { 
    for(i in 1:size_pop)
    {
      if(!is.null(seed))
      {
        set.seed(seed = seed + i) 
      }
      pop[[i]] <- sort(sample(x = size_world, size = size_ind, replace = FALSE))
    } # end for
  } # end else
  
  # evaluate the population will send a population of model objects
  #pop <- evaluatePopulation(X, y, clf, pop, force.re.evaluation = TRUE, eval.all = FALSE)
  
  #####TODO here try to add the models from the file with the right sparsity 
  if(clf$params$popSourceFile != "NULL")
  {
    popFromFile <- loadPopulation(clf$params$popSourceFile)
    pop2Add <- list()
    j <- 1
    for(i in 1:length(popFromFile))
    {
      if(popFromFile[[i]]$eval.sparsity == clf$params$current_sparsity)
      {
        pop2Add[[j]] <- popFromFile[[i]]$indices_
        j <- j + 1
      }
    }
    if(length(pop2Add) > 0)
    {
      if(length(pop2Add) > length(pop))
      {
        stop("The source population you're trying to use is larger than the population you use")
      } else
      {
        pop[(length(pop) - length(pop2Add) +1):length(pop)] <- pop2Add
      }
    } 
  }
  
  # add the best ancestors of the previous k-sparsity
  if (!is.null(best_ancestor))
  {
    if(length(best_ancestors) < number_best_ancestors) # when few variables this can happen
    {
      pop <- c(best_ancestors, pop)
    }
  }
  
  return(pop)
}


#' Creates new combinations of features based from a parents. 
#'
#' @description This function is used in terga1 will create new combinations of features based of existing ones from the parents.
#' @param clf: the classifier parameter object
#' @param pop: A population (i.e. list) of index vectors
#' @param parents: Indexes of the population pointing to the subset of the population containing the parents (whose genes/features) will be used to create the children.
#' @param seed: For reproductibility purpose to fix the random generator number.
#' @return a population of models, containing parents and children
#' @export
crossing <- function(clf, pop, parents, seed=NULL) 
{
  # save the parents
  old_pop <- pop[parents]
  # select the needed number of couples that are to be created. Each couple will give a child
  number_couples <- length(pop) - length(parents)
 
  # Use the population combinarory algorithm to generate combinations of (couples)
  couples <- population(clf = clf, 
                        size_ind = 2, 
                        size_world = length(parents), 
                        best_ancestor = NULL, 
                        size_pop = number_couples,
                        seed = seed)
  
  children <- list()
  # parallel computing
  if(clf$params$parallel.local)
  {
    children <- foreach(i = 1:length(couples))  %dorng% 
    { # returns a list by default
      parent1 <- old_pop[[couples[[i]][1]]]
      parent2 <- old_pop[[couples[[i]][2]]]
      parents_gene_reservoir <- sort(unique(c(parent1, parent2)))
      # These are not necessarely unique features (depends on the method used to build the )
      
      #### Partie qui doit changer pour les denseVec
      if(!is.null(seed))
      {
        set.seed(seed+i)
      }
  
      sample(x = parents_gene_reservoir, 
             size = clf$params$current_sparsity, 
             replace = FALSE)
    }
  }else
  {
    for (i in 1:length(couples)) 
    {
      parent1 <- old_pop[[couples[[i]][1]]]
      parent2 <- old_pop[[couples[[i]][2]]]
      parents_gene_reservoir <- sort(unique(c(parent1, parent2)))
      # unique vars
      
      #### Partie qui doit changer pour les denseVec
      if(!is.null(seed))
      {
        set.seed(seed+i)
      }

      child <- sample(x = parents_gene_reservoir, 
                      size = clf$params$current_sparsity, 
                      replace = FALSE)
      children[[i]] <- child
    }
  }
  return(c(old_pop, children))
}


#' Changes feature indexes in a given percentage of models. 
#'
#' @description This function is used in terga1 will create new combinations of features based of existing ones from the parents.
#' @param clf: the classifier parameter object
#' @param pop: A population (i.e. list) of index vectors
#' @param selection: Indexes of the population pointing to the subset of the models to be changed
#' @param seed: For reproductibility purpose to fix the random generator number.
#' @return a population of models among which the mutated ones
#' @export
mutate <- function(clf, pop, selection, seed = NULL) 
{
  mutated_pop <- pop
  # degree of mutation (number of genes to be mutated)
  perc <- ceiling(clf$params$current_sparsity * clf$params$mutate_rate/100)
  all_genes <- (1:clf$params$size_world) # all genes
  
  if(clf$params$parallel.local)
  { # if // computing
    res <- foreach(i = selection)  %dorng% 
    { # returns a list by default
      individual_to_be_mutated <- sort(pop[[i]])
      
      # the size of the reservoir
      size_reservoir <- clf$params$size_world - clf$params$current_sparsity
      # if the gene reservoir is smaller than the number of genes to mutate, 
      # lower the mutation rate to the size of the reservoir
      if(perc > size_reservoir){ perc <- size_reservoir }
      
      if(!is.null(seed))
      {
        set.seed(seed+i)
      }
      # identify the indexes of the genes to be mutated
      index_genes_to_mutate <- sample(clf$params$current_sparsity, perc)
      # the unique genes remaining in the reservoir that are not in the untouched part of the genome
      # It is to be noted that the genes to mutate are part of this reservoir and some of them can be 
      # picked again, which lowers the percentage of mutation
      remaining_unique_genes <- all_genes[!all_genes %in% individual_to_be_mutated[-index_genes_to_mutate]]
      
      if(!is.null(seed))
      {
        set.seed(seed+i)
      }
      # draw new genes (some of them might be old)
      new_genes <- sample(x = remaining_unique_genes, size = perc)
      # the mutated individual
      mutated_individual <- individual_to_be_mutated
      mutated_individual[index_genes_to_mutate] <- new_genes; mutated_individual <- sort(mutated_individual)
      # save the individual
    }
    mutated_pop[selection] <- res 
  }else
  {
    for (i in selection) 
    {
      individual_to_be_mutated <- sort(pop[[i]])
      # the size of the reservoir
      size_reservoir <- clf$params$size_world - clf$params$current_sparsity
      # if the gene reservoir is smaller than the number of genes to mutate, 
      # lower the mutation rate to the size of the reservoir
      if(perc > size_reservoir){
        perc <- size_reservoir
      }
      
      if(!is.null(seed))
      {
        set.seed(seed+i)
      }
      # identify the indexes of the genes to be mutated
      index_genes_to_mutate <- sample(clf$params$current_sparsity, perc)
      # the unique genes remaining in the reservoir that are not in the untouched part of the genome
      # It is to be noted that the genes to mutate are part of this reservoir and some of them can be 
      # picked again, which lowers the percentage of mutation
      remaining_unique_genes <- all_genes[!all_genes %in% individual_to_be_mutated[-index_genes_to_mutate]]
      
      if(!is.null(seed))
      {
        set.seed(seed+i)
      }
      # draw new genes (some of them might be old)
      new_genes <- sample(x = remaining_unique_genes, size = perc)
      # the mutated individual
      mutated_individual <- individual_to_be_mutated
      mutated_individual[index_genes_to_mutate] <- new_genes; mutated_individual <- sort(mutated_individual)
      # save the individual
      mutated_pop[[i]] <- mutated_individual
    }
  }
  return(mutated_pop)
}


#' Creates new combinations of features based from a parents. 
#'
#' @description This function is used in terga1 and is the main engine of the algorithm that allows to cross, mutate and select individuals from one generation to the next.
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the response vector
#' @param clf: the classifier parameter object
#' @param pop: A population (i.e. list) of index vectors
#' @param seed: For reproductibility purpose to fix the random generator number.
#' @return a population of models, containing parents and children
#' @export
evolve <- function(X, y, clf, pop, seed = NULL) 
{
  # store information on the evolution
  evolved_pop <- pop
  trace_evolution <- rep(NA, clf$params$nb_generations)
  
  if(length(pop[[1]]) < nrow(X) - 1) # testing zize TODO verify why
  {
    for (i in 1:clf$params$nb_generations) # For all the generations
    {
      # 1. Evaluate the current population
      # transform to a population of model objects
      evolved_pop.mod <- listOfSparseVecToListOfModels(X = X, y = y, clf = clf, v = evolved_pop)
      # evaluate the population
      evolved_pop.eval <- evaluatePopulation(X, y, clf, evolved_pop.mod, eval.all = FALSE)
      # clean the population
      evolved_pop.eval <- cleanPopulation(pop = evolved_pop.eval, clf = clf)
      # get the evaluation
      evaluation <- as.numeric(populationGet_X(element2get = "fit_", toVec = TRUE, na.rm = FALSE)(pop = evolved_pop.eval))
      # clean NA objects
      toclean <- is.na(evaluation)
      evolved_pop.eval <- evolved_pop.eval[!toclean] # reorder to make it easier to check
      evaluation <- evaluation[!toclean]
      # order
      evaluation.ord <- order(abs(evaluation), decreasing = TRUE)
      evolved_pop.eval <- evolved_pop.eval[evaluation.ord] # reorder to make it easier to check
      evaluation <- evaluation[evaluation.ord]
      if(length(evolved_pop.eval)==0)
      {
        break
      }
      # get best
      best_individual_index <- which.max(abs(evaluation))
      
      # transform the models back to list of sparse indexes
      evolved_pop <- listOfModelsToListOfSparseVec(list.models = evolved_pop.eval)
      if(any(is.na(evaluation)))
      {
        warning("evolve: There are individuals with evaluation = NA")
      }
      
      if(clf$params$debug)
      {
        print("Before evolution")
        plot(evaluation, ylim=c(0,1), col="red", pch=19)
        print(paste("length population:",length(evolved_pop)))
        print(paste("length evaluation:",length(evolved_pop.eval)))
        print(paste("best index:",best_individual_index))
        try(cat(paste("gen =",i,"\t", printModel(mod = evolved_pop.eval[[which.max(evaluation)]], method = clf$params$print_ind_method, score = "fit_"),"\n")), silent = TRUE)
      }
      
      # store the best performance for each generation
      trace_evolution[i] <- evaluation[best_individual_index]
      
      # 2 Select the parents of the next generation (indexes in the population). The best should be in the parents
      parents <- NA
      if (clf$params$select_type == "elite")    { parents <- selectElite(evaluation, 50)}
      if (clf$params$select_type == "tournoi")  { parents <- selectTournoi(evaluation, 50, discrete = F)}
      if (clf$params$select_type == "mixed")    { parents <- selectMixed(evaluation, 
                                                                         clf$params$select_perc1, 
                                                                         clf$params$select_perc2)}
      if(clf$params$debug)
      {
        if(best_individual_index %in% parents)
        {
          print(paste("Best in parents:",TRUE))
        }else 
        {
          print(paste("Best in parents:",FALSE))
        }
      }
      
      # 3. Cross parents genetic material to give the children and put back to the population
      evolved_pop <- crossing(clf, evolved_pop, parents, seed = seed)
      
      # Fill up the population size by adding random if some are missing
      if(length(evolved_pop) < clf$params$size_pop)
      {
        nb.tomake <- clf$params$size_pop - length(evolved_pop)
        if(clf$params$debug & nb.tomake>0) 
        {
          print(paste("Adding", nb.tomake, "new individuals,"))
        }
        
        evolved_pop <- c(evolved_pop,
                         population(clf=clf, 
                                    size_ind = clf$params$current_sparsity, 
                                    size_world=nrow(X), 
                                    best_ancestor = NULL, 
                                    size_pop = nb.tomake,
                                    seed=clf$params$current_seed)
                         )
      }
      
      # transform to a population of model objects
      evolved_pop.mod <- listOfSparseVecToListOfModels(X = X, y = y, clf = clf, v = evolved_pop)
      # evaluate the population
      evolved_pop.eval <- evaluatePopulation(X, y, clf, evolved_pop.mod, force.re.evaluation = TRUE, eval.all = FALSE)
      # clean the population
      evolved_pop.eval <- cleanPopulation(pop = evolved_pop.eval, clf = clf)
      # get the evaluation vector
      evaluation <- as.numeric(populationGet_X(element2get = "fit_", toVec = TRUE, na.rm = FALSE)(pop = evolved_pop.eval))
      # clean NA objects
      toclean <- is.na(evaluation)
      evolved_pop.eval <- evolved_pop.eval[!toclean] # reorder to make it easier to check
      evaluation <- evaluation[!toclean]
      # order
      evaluation.ord <- order(abs(evaluation), decreasing = TRUE)
      evolved_pop.eval <- evolved_pop.eval[evaluation.ord] # reorder to make it easier to check
      evaluation <- evaluation[evaluation.ord]
      if(length(evolved_pop.eval)==0)
      {
        break
      }
      # get best
      best_individual_index <- which.max(abs(evaluation))
      # take into account the missing ones
      # convert back to list of sparse vectors
      evolved_pop <- listOfModelsToListOfSparseVec(list.models = evolved_pop.eval)
      if(any(is.na(evaluation)))
      {
        warning("evolve: There are individuals with evaluation = NA")
      }
      if(clf$params$debug)
      {
        print("After crossing")
        plot(evaluation, ylim=c(0,1), col="red", pch=19)
        print(paste("length population:",length(evolved_pop)))
        print(paste("length evaluation:",length(evolved_pop.eval)))
        print(paste("best index:",best_individual_index))
        try(cat(paste("gen =",i,"\t", printModel(mod = evolved_pop.eval[[which.max(evaluation)]], method = clf$params$print_ind_method, score = "fit_"),"\n")), silent = TRUE)
      }
      
      # 4. Mutate the new population
      # Select the samples to mutate (and protect the best one/s)
      protect.ind <- order(abs(evaluation), decreasing = TRUE)[1:5] # keep the 5 best samples untouched
      selection <- !(1:length(evaluation)) %in% protect.ind
      if(!is.null(seed))
      {
        set.seed(seed = seed)  
      }
      selection.ind <- sample(x = which(selection), 
                              size = ceiling(sum(selection) * clf$params$mutate_size/100),
                              replace = FALSE)
      # MUTATE
      evolved_pop <- predomics::mutate(clf, evolved_pop, selection.ind, seed = seed)
      # transform to a population of model objects
      evolved_pop.mod <- listOfSparseVecToListOfModels(X = X, y = y, clf = clf, v = evolved_pop)
      # evaluate the population
      evolved_pop.eval <- evaluatePopulation(X, y, clf, evolved_pop.mod, force.re.evaluation = TRUE, eval.all = FALSE)
      # clean the population
      evolved_pop.eval <- cleanPopulation(pop = evolved_pop.eval, clf = clf)
      # get the evaluation
      evaluation <- as.numeric(populationGet_X(element2get = "fit_", toVec = TRUE, na.rm = FALSE)(pop = evolved_pop.eval))
      # clean NA objects
      toclean <- is.na(evaluation)
      evolved_pop.eval <- evolved_pop.eval[!toclean] # reorder to make it easier to check
      evaluation <- evaluation[!toclean]
      # order
      evaluation.ord <- order(abs(evaluation), decreasing = TRUE)
      evolved_pop.eval <- evolved_pop.eval[evaluation.ord] # reorder to make it easier to check
      evaluation <- evaluation[evaluation.ord]
      if(length(evolved_pop.eval)==0)
      {
        break
      }
      # get best
      best_individual_index <- which.max(abs(evaluation))
      # transform back to a list of sparse models
      evolved_pop <- listOfModelsToListOfSparseVec(list.models = evolved_pop.eval)
      if(any(is.na(evaluation)))
      {
        warning("evolve: There are individuals with evaluation = NA")
      }
      if(clf$params$debug)
      {
        print("After mutation")
        plot(evaluation, ylim=c(0,1), col="red", pch=19)
        print(paste("length population:",length(evolved_pop)))
        print(paste("length evaluation:",length(evolved_pop.eval)))
        print(paste("best index:",best_individual_index))
        try(cat(paste("gen =",i,"\t", printModel(mod = evolved_pop.eval[[which.max(evaluation)]], method = clf$params$print_ind_method, score = "fit_"),"\n")), silent = TRUE)
      }
      
      # print best
      if (clf$params$debug)
      {
        # Print out the best individual
        best_individual <- evolved_pop.eval[[best_individual_index]]
        if(isModel(best_individual))
        {
          try(cat(paste("gen =",i,"\t", printModel(mod = best_individual, method = clf$params$print_ind_method, score = "fit_"),"\n")), silent = TRUE)
        }
      }
      
      # 5. Test convergence and exit if atteint
      if (clf$params$convergence) # if we are testing convergence
      {  
        trace.nona <- trace_evolution[!is.na(trace_evolution)]
        if (length(trace.nona) > clf$params$convergence_steps) # run at least enough steps to test the convergence
        {  
          if (convergence.test(x = trace.nona, steps = clf$params$convergence_steps)) # if convergence atteint
          {  
            break("")
          } # end if convergence
        }
      }# if convergence
    } # end of generation

    if (clf$params$plot)  
    { 
      if(!all(is.na(trace_evolution))) # For k = 1 in ter language there is no valid models so it won't be able to plot
      {
        plot(trace_evolution[!is.na(trace_evolution)], ylab = "AUC score",  
             pch = 20,col = "red",type = "l",ylim = c(0,1), 
             main = paste("Best score evolution: k=",length(pop[[1]]),sep = ""))   
      }
    }
  }
  return(evolved_pop)
}



# function needed for testing the convergence of a generation
convergence.test <- function(x, steps = 10) {
  test <- FALSE
  x <- x[!is.na(x)]
  if (length(x) < steps) { stop("Convergence steps number greater than the generations")}
  count0 <- 0
  for (k in 1:(length(x) - 1)) {
    #print(x[k])
    if ((x[k + 1]  - x[k]) == 0) {
      count0 <- count0 + 1;
      #print(paste("up", count0))
    }
    else  {count0 <- 0} # reset the counter
    if (count0 == steps)  {
      test <- TRUE
      break("converged")
    }
  }
  return(test)
}


#######################################################
###  CLASSE INDIVIDUAL #####
#######################################################
# this function will estimate the coefficients of multiple logistic regression model
estimateCoefficientsIndividual <- function(X, y, ind, type  =  "logit") {
  trait.bin <- c(0,1)[as.factor(y)]
  res <- list()
  if (type == "logit") 
  {
    # testing AUC
    if (length(ind) == 1) 
    {
      dat <- as.data.frame(X[ind,])
      m <- glm(trait.bin ~ ., data = dat, family = binomial(logit))
      #summary(m)
      coefs <- m$coefficients[2:length(m$coefficients)]
      aic <- m$aic
    }else
    {
      dat <- as.data.frame(t(X[ind,]))
      m <- glm(trait.bin ~ ., data = dat, family = binomial(logit))
      #summary(m)
      coefs <- m$coefficients[2:length(m$coefficients)]
      aic <- m$aic
    }
    names(coefs) <- ind
    res$coefs <- coefs
    res$aic <- aic
    res$glm <- m
  }else {
    print("estimateCoefficientsIndividual: Method unknown !!!")
  }
  return(res)
}



#######################################################
###  CLASSE PRINCIPALE #####
#######################################################

# (c) authors, Edi Prifti and Emmanuelle le Chatelier.

#######################################################
# General functions for result exploration
#######################################################
# modelSampling : function that select different training/test datasets randomly or by a sliding window
modelSampling <- function(vect, sampling.type =  c("window","random"), training.pc =  90, iteration =  10)  {
  samples <- sample(colnames(vect),size =  ncol(vect),replace =  FALSE) # random permutations of the cohort samples
  sampling <- matrix(0,nrow =  ncol(vect), ncol =  iteration)
  rownames(sampling) <- samples
  colnames(sampling) <- 1:iteration
  size.train <- round(ncol(vect)*training.pc*0.01)
  size.valid <- nrow(sampling) - size.train
  
  if (sampling.type   ==    "window")  {
    for (i in 1:iteration)  {
      beg <- ((i - 1)*size.valid) + 1
      if (beg > nrow(sampling))   {
        beg <- beg %% nrow(sampling)
      }
      end <- beg - 1 + size.train
      if (end <=  nrow(sampling))  {
        sampling[beg:end,i] <- 1
      } else   {
        sampling[beg:nrow(sampling),i] <- 1
        sampling[0:(end %% nrow(sampling)),i] <- 1
      }
    }
  }
  if (sampling.type   ==    "random")  {
    for (i in 1:iteration)  {
      temp <- sample(x =  c(1:nrow(sampling)), size =  size.train, replace =  FALSE)
      sampling[temp, i] <- 1 # train cohort
    }
  }
  return(sampling)
}

tableBestModels <- function(model.sim)  {                          ## corr 03/2014
  N <- c(rep(NA,min(20,length(model.sim))))
  AUC <- c(rep(NA,min(20,length(model.sim))))
  MGS <- c(rep(NA,min(20,length(model.sim))))
  tableau <- data.frame(cbind(N, AUC, MGS))
  for (i in 1:length(model.sim))  {
    tableau[i,1] <- paste("N",i,sep =  "")
    ana <- paste("model_N",i,sep =  "")
    tableau[i,2] <- round(model.sim[[ana]][[grep("score",names(model.sim[[ana]]))]],4)
    tableau[i,3] <- paste(c(sort(model.sim[[ana]][[grep("model",names(model.sim[[ana]]))]])), collapse  =   " & ")
  }
  return(tableau)
}

