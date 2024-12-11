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
#' Select Elite Individuals Based on Evaluation Scores
#'
#' This function selects the top-performing individuals from a population based
#' on their evaluation scores. It returns the indices of the highest-scoring
#' individuals, making up a specified percentage of the total population.
#'
#' @param evaluation A numeric vector of evaluation scores for each individual
#'   in the population. Higher scores indicate better performance.
#' @param percentage A numeric value representing the percentage of the top
#'   individuals to select, specified as a value between 0 and 100.
#'
#' @details The function calculates the number of individuals to select based on
#' the given `percentage` and the length of `evaluation`. It then sorts
#' `evaluation` in descending order and returns the indices of the
#' top-performing individuals.
#'
#' This selection method is often used to maintain high-quality solutions
#' (elites) in evolutionary algorithms, ensuring that the best individuals are
#' retained across generations.
#'
#' @return A vector of indices representing the top-performing individuals in
#'   the population.
#'
#' @examples
#' \dontrun{
#' evaluation <- c(0.8, 0.6, 0.9, 0.7, 0.5)
#' elite_indices <- selectElite(evaluation, percentage = 40)
#' print(elite_indices)  # Returns indices of the top 40% based on evaluation scores
#' }
#'
#' @export
selectElite <- function(evaluation, percentage) 
{
  perc <- ceiling(percentage*length(evaluation)/100)
  res <- order(evaluation,decreasing = T)[1:perc]
  return(res)
}


#' Select Individuals Using Tournament Selection
#'
#' This function performs tournament selection on a population, randomly
#' selecting a specified percentage of individuals for the next generation. The
#' selection can be either a fixed number (`discrete = TRUE`) or a percentage of
#' the population size.
#'
#' @param selection A vector representing the indices or identifiers of the
#'   individuals from which to select.
#' @param percentage A numeric value indicating either the exact number (if
#'   `discrete = TRUE`) or percentage (if `discrete = FALSE`) of individuals to
#'   select.
#' @param discrete A logical value indicating whether `percentage` should be
#'   treated as a discrete number of individuals (`TRUE`) or a percentage of the
#'   total population (`FALSE`).
#'
#' @details If `discrete` is set to `TRUE`, the function selects exactly
#' `percentage` individuals from `selection`. If `discrete` is `FALSE`, the
#' function interprets `percentage` as a proportion of the population size, and
#' the number of individuals selected is calculated accordingly.
#'
#' This selection method introduces randomness while maintaining a proportionate
#' selection of individuals, which can help maintain diversity within the
#' population in evolutionary algorithms.
#'
#' @return A vector of selected individuals from `selection`.
#'
#' @examples
#' \dontrun{
#' selection <- 1:10
#' # Select 3 individuals directly
#' selected_indices <- selectTournoi(selection, percentage = 3, discrete = TRUE)
#' print(selected_indices)
#'
#' # Select 30% of individuals from the population
#' selected_indices <- selectTournoi(selection, percentage = 30, discrete = FALSE)
#' print(selected_indices)
#' }
#'
#' @export
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


#' Select Individuals Using Mixed Elite and Tournament Selection
#'
#' This function combines elite selection and tournament selection to select
#' individuals from a population. It first selects a specified percentage of
#' top-performing individuals (elite) and then fills the remaining percentage
#' with randomly selected individuals (tournament).
#'
#' @param evaluation A numeric vector of evaluation scores for each individual
#'   in the population, with higher values indicating better performance.
#' @param percentage_elite A numeric value representing the percentage of
#'   individuals to select based on elite selection (top performers).
#' @param percentage_tournoi A numeric value representing the percentage of
#'   individuals to select based on tournament selection (random selection).
#'
#' @details The function calculates the total number of individuals to select
#' (`perc`) based on `percentage_elite` and `percentage_tournoi`. It then uses
#' `selectElite` to select the top-performing individuals as elites, filling
#' `perc_elite` slots. Afterward, it uses `selectTournoi` to fill the remaining
#' `perc_tournoi` slots by selecting randomly from the remaining individuals.
#'
#' This mixed selection approach ensures that high-quality solutions are
#' retained while also introducing diversity through random selection.
#'
#' @return A vector of indices representing the selected individuals, combining
#'   elite and tournament selections.
#'
#' @examples
#' \dontrun{
#' evaluation <- c(0.9, 0.85, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05)
#' selected_indices <- selectMixed(evaluation, percentage_elite = 30, percentage_tournoi = 20)
#' print(selected_indices)  # Returns indices selected via elite and tournament strategies
#' }
#'
#' @export
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


#' Generate a Population with Controlled Sparsity Using Dense Vectors
#'
#' This function generates a population of dense vectors for a given size of the
#' "world." It controls the sparsity of each vector based on the parameters
#' defined in `clf`, either by creating random dense vectors or by building on a
#' given `best_ancestor` to produce similar solutions.
#'
#' @param clf A classifier object containing parameters for population
#'   generation, including:
#'   - `size_pop`: Size of the population to generate.
#'   - `current_sparsity`: Target sparsity level for the population (if specified).
#'   - `sparsity`: Vector of sparsity levels, from which the minimum is used if `current_sparsity` is not provided.
#'   - `perc_best_ancestor`: Percentage of the population to generate based on the best ancestor.
#' @param size_world An integer representing the size of the dense vector for
#'   each individual.
#' @param best_ancestor An optional vector representing the best ancestor to use
#'   as a basis for generating a portion of the population. If provided, this
#'   ancestor must have a non-zero gene count equal to `aimedSpar - 1`.
#'
#' @details The function: 1. Determines the target sparsity (`aimedSpar`) based
#' on `current_sparsity` or the minimum value in `sparsity`. 2. Generates a set
#' of dense vectors based on `aimedSpar`:
#'    - **Random Vectors**: For the general population, randomly generated with a specified probability for non-zero elements.
#'    - **Ancestor-based Vectors**: If `best_ancestor` is provided, some population members are generated by modifying the ancestor with additional random genes.
#'
#' **Sparsity Control**: Each generated vector's sparsity is adjusted to exactly match `aimedSpar` by adding or removing non-zero elements as necessary.
#'
#' @return A list of dense vectors representing the generated population, each
#'   with controlled sparsity.
#'
#' @examples
#' \dontrun{
#' clf <- list(
#'   params = list(size_pop = 10, current_sparsity = 5, sparsity = c(3, 5, 7), perc_best_ancestor = 50)
#' )
#' size_world <- 20
#' best_ancestor <- sample(c(-1, 0, 1), size_world, replace = TRUE)
#' population <- populationDenseVec(clf, size_world, best_ancestor)
#' print(population)
#' }
#'
#' @export
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


#' Generate a Population with Controlled Sparsity and Unique Individuals
#'
#' This function generates a population of individuals (feature subsets) with
#' specified sparsity, optional use of a best ancestor, and controls for unique
#' individuals. It creates feature subsets according to specified parameters,
#' optionally seeding some individuals based on a provided best ancestor.
#'
#' @param clf A classifier object containing parameters for population
#'   generation, including:
#'   - `size_pop`: Size of the population to generate.
#'   - `perc_best_ancestor`: Percentage of the population to generate based on the best ancestor.
#'   - `unique_vars`: Boolean indicating whether each individual in the population should be unique.
#'   - `popSourceFile`: Path to a file containing saved populations to import.
#'   - `current_sparsity`: Target sparsity for each individual in the population.
#' @param size_ind Integer specifying the number of features (sparsity) for each
#'   individual.
#' @param size_world Integer representing the size of the "world" (total number
#'   of possible features).
#' @param best_ancestor Optional vector or list representing the best ancestor
#'   to use as a template for generating a portion of the population.
#' @param size_pop Optional integer specifying the number of individuals in the
#'   population. Defaults to `clf$params$size_pop`.
#' @param seed Optional integer seed for random number generation, ensuring
#'   reproducibility.
#'
#' @details The function generates a population by: 1. **Best Ancestor**: If
#' `best_ancestor` is provided, a portion of the population is generated by
#' slightly modifying the ancestor. 2. **Population Size**: Determines the final
#' population size based on combinatory limits and `size_pop`. 3.
#' **Uniqueness**: If `unique_vars` is `TRUE`, each individual is checked for
#' uniqueness before adding it to the population.
#'
#' **Additional Options**:
#' - The function can import previously saved populations from `popSourceFile` and add them to the current population.
#' - Ensures that population size and sparsity constraints are respected.
#'
#' @return A list of individuals (feature subsets), each represented as a sorted
#'   vector of selected feature indices.
#'
#' @examples
#' \dontrun{
#' clf <- list(
#'   params = list(size_pop = 10, perc_best_ancestor = 20, unique_vars = TRUE, popSourceFile = "NULL", current_sparsity = 5)
#' )
#' size_ind <- 5
#' size_world <- 20
#' best_ancestor <- sample(1:size_world, size_ind - 1, replace = FALSE)
#' population <- population(clf, size_ind, size_world, best_ancestor)
#' print(population)
#' }
#'
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


#' Generate Children by Crossover from Parent Population
#'
#' This function performs crossover on a population by selecting pairs of parent
#' individuals to generate children. Each child is created by combining unique
#' genes from two parents, up to a specified sparsity level.
#'
#' @param clf A classifier object containing parameters for crossover,
#'   including:
#'   - `current_sparsity`: The target sparsity level for each child.
#'   - `parallel.local`: Boolean indicating whether parallel processing should be used.
#' @param pop A list representing the population of individuals (feature
#'   subsets), where each individual is represented as a vector of feature
#'   indices.
#' @param parents A vector of indices representing the parent individuals within
#'   `pop` from which to generate children.
#' @param seed Optional integer seed for random number generation, ensuring
#'   reproducibility.
#'
#' @details The function: 1. Selects a specified number of couples from the
#' parents to generate children. 2. Combines unique genes from each pair of
#' parents to create a child with `current_sparsity` number of features. 3.
#' Optionally uses parallel processing to accelerate child generation.
#'
#' For each child, the genes are selected from the combined unique genes of both
#' parents (`parents_gene_reservoir`). The final population consists of both the
#' original parents and the generated children.
#'
#' @return A list representing the updated population, including both the
#'   original parents and the newly generated children.
#'
#' @examples
#' \dontrun{
#' clf <- list(
#'   params = list(current_sparsity = 5, parallel.local = FALSE)
#' )
#' pop <- list(1:5, 2:6, 3:7, 4:8)  # Example parent population
#' parents <- c(1, 2, 3, 4)
#' new_population <- crossing(clf, pop, parents, seed = 42)
#' print(new_population)
#' }
#'
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


#' Mutate Selected Individuals in a Population
#'
#' This function applies mutation to a selected set of individuals in a
#' population by randomly modifying a specified percentage of genes. The
#' mutation changes some genes to new, unique values from the gene pool.
#'
#' @param clf A classifier object containing parameters for mutation, including:
#'   - `current_sparsity`: The target sparsity level for each individual.
#'   - `mutate_rate`: The mutation rate as a percentage of the total genes in each individual.
#'   - `size_world`: The total number of genes available (size of the "world").
#'   - `parallel.local`: Boolean indicating whether parallel processing should be used.
#' @param pop A list representing the population, where each individual is a
#'   vector of gene indices.
#' @param selection A vector of indices indicating which individuals in `pop`
#'   should undergo mutation.
#' @param seed Optional integer seed for random number generation to ensure
#'   reproducibility.
#'
#' @details The function: 1. Calculates the number of genes to mutate based on
#' `current_sparsity` and `mutate_rate`. 2. For each selected individual,
#' identifies genes to mutate and replaces them with new, unique genes from the
#' remaining gene pool. 3. Uses parallel processing if specified by
#' `parallel.local`.
#'
#' **Mutation Process**:
#' - Genes are randomly selected for mutation from each individual’s current genes.
#' - New genes are randomly chosen from a "reservoir" of available genes not currently in the individual, ensuring they are unique and do not exceed the specified mutation percentage.
#'
#' @return A list representing the population after mutation, with selected
#'   individuals modified.
#'
#' @examples
#' \dontrun{
#' clf <- list(
#'   params = list(current_sparsity = 5, mutate_rate = 20, size_world = 10, parallel.local = FALSE)
#' )
#' pop <- list(1:5, 2:6, 3:7, 4:8)  # Example population
#' selection <- c(1, 3)  # Indices of individuals to mutate
#' mutated_population <- mutate(clf, pop, selection, seed = 42)
#' print(mutated_population)
#' }
#'
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


#' Evolve a Population of Models Over Generations
#'
#' This function evolves a population of models across a specified number of generations, using selection, crossover, and mutation operations. It tracks the best-performing individuals in each generation to encourage convergence toward optimal solutions.
#'
#' @param X A matrix or data frame representing the features of the training dataset.
#' @param y A vector representing the target variable for the training dataset.
#' @param clf A classifier object containing parameters for the evolutionary algorithm, including:
#'   - `nb_generations`: The number of generations to evolve the population.
#'   - `size_pop`: The target population size.
#'   - `current_sparsity`: The desired sparsity level for individuals.
#'   - `select_type`: Selection method, which can be "elite," "tournoi," or "mixed."
#'   - `select_perc1` and `select_perc2`: Percentages used in mixed selection.
#'   - `mutate_size`: Percentage of the population to mutate in each generation.
#'   - `plot`: Boolean indicating whether to plot the evolution of the best score.
#'   - `convergence`: Boolean indicating if convergence testing should be applied.
#'   - `convergence_steps`: The number of steps used for convergence testing.
#' @param pop A list representing the initial population of models, with each individual represented as a vector of selected features.
#' @param seed Optional integer seed for random number generation to ensure reproducibility.
#'
#' @details
#' The function evolves a population of models over multiple generations using the following steps:
#' 1. **Evaluation**: Evaluates each individual in the population.
#' 2. **Selection**: Selects parents based on specified selection criteria (`select_type`).
#' 3. **Crossover**: Combines genetic material from selected parents to produce children.
#' 4. **Mutation**: Mutates a portion of the population to introduce genetic diversity.
#' 5. **Convergence Testing**: Optionally tests for convergence, stopping the evolution early if stability is detected in performance.

#' **Population Structure**: The population is a list of vectors, where each vector represents a sparse model (a subset of features). After each generation, the function retains the best-performing individuals, tracking the evolution of the best score across generations.

#' @return A list representing the final evolved population of models.
#'
#' @examples
#' \dontrun{
#' clf <- list(
#'   params = list(nb_generations = 10, size_pop = 20, current_sparsity = 5, 
#'                 select_type = "elite", mutate_size = 10, plot = TRUE, convergence = FALSE)
#' )
#' X <- matrix(runif(200), nrow = 20)
#' y <- sample(0:1, 20, replace = TRUE)
#' initial_pop <- list(1:5, 2:6, 3:7)  # Example initial population
#' final_pop <- evolve(X, y, clf, initial_pop, seed = 42)
#' print(final_pop)
#' }
#'
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



#' Test for Convergence in Evolutionary Algorithm Scores
#'
#' This function checks for convergence in a vector of scores, identifying if
#' there has been a stable sequence of scores over a specified number of steps.
#' Convergence is defined as no change in score for `steps` consecutive
#' generations.
#'
#' @param x A numeric vector representing scores from successive generations in
#'   an evolutionary algorithm.
#' @param steps An integer specifying the number of consecutive steps required
#'   with no change in score to consider the sequence as converged.
#'
#' @details The function iterates through the scores in `x` and counts
#' consecutive instances where there is no change in score. If the count reaches
#' `steps`, the function returns `TRUE`, indicating convergence. If there are
#' fewer scores in `x` than `steps`, the function stops with an error.
#'
#' **Usage in Evolutionary Algorithms**: This convergence test is useful in evolutionary algorithms to detect when further evolution is unlikely to yield improved scores, allowing the algorithm to terminate early.
#'
#' @return A logical value: `TRUE` if convergence is detected, `FALSE`
#'   otherwise.
#'
#' @examples
#' \dontrun{
#' scores <- c(0.8, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85)
#' has_converged <- convergence.test(scores, steps = 5)
#' print(has_converged)  # Returns TRUE, as the score did not change for 5 steps
#' }
#'
#' @export
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


#' Estimate Model Coefficients for a Subset of Features
#'
#' This function estimates coefficients for a logistic regression model
#' (`logit`) on a subset of features specified by `ind`. It evaluates the
#' individual model’s performance by calculating the Akaike Information
#' Criterion (AIC) and returns the coefficients.
#'
#' @param X A matrix or data frame representing the features of the dataset.
#' @param y A vector representing the binary target variable for the dataset.
#' @param ind A vector of feature indices specifying the subset of features to
#'   include in the model.
#' @param type A character string specifying the type of model to fit.
#'   Currently, only `"logit"` (logistic regression) is supported.
#'
#' @details The function fits a logistic regression model using the features
#' specified in `ind`. If `ind` contains a single feature, the function formats
#' `X` accordingly. It then calculates the coefficients for each feature and the
#' AIC of the model, providing insights into model performance.
#'
#' **Note**: Currently, only `"logit"` models are supported. Other types will result in an error message.
#'
#' @return A list containing:
#'   - `coefs`: A named vector of estimated coefficients for each feature in `ind`.
#'   - `aic`: The AIC value for the fitted logistic model, indicating model fit.
#'   - `glm`: The fitted `glm` object.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100), nrow = 10)
#' y <- sample(0:1, 10, replace = TRUE)
#' ind <- c(1, 2, 3)  # Example subset of feature indices
#' result <- estimateCoefficientsIndividual(X, y, ind, type = "logit")
#' print(result$coefs)  # Estimated coefficients
#' print(result$aic)    # AIC of the model
#' }
#'
#' @export
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

#' Generate Training and Validation Sampling Matrices
#'
#' This function creates a sampling matrix for cross-validation, specifying
#' which samples should be used for training and validation across multiple
#' iterations. It supports two sampling types: "window" (systematic shifting)
#' and "random" (random selection).
#'
#' @param vect A matrix or data frame with samples as columns, used to determine
#'   the sampling order and dimensions.
#' @param sampling.type A character string indicating the sampling method,
#'   either `"window"` or `"random"`. Default is `c("window", "random")`.
#'   - `"window"`: Sequentially shifts the training window across samples.
#'   - `"random"`: Randomly selects training samples in each iteration.
#' @param training.pc An integer representing the percentage of samples to be
#'   allocated for training in each iteration. Default is 90.
#' @param iteration An integer specifying the number of sampling iterations.
#'   Default is 10.
#'
#' @details The function creates a binary sampling matrix with rows as sample
#' names and columns as iterations. Cells with `1` indicate the sample is in the
#' training set, while `0` indicates the sample is in the validation set. For
#' `"window"` sampling, a moving window of training samples shifts across each
#' iteration. For `"random"` sampling, training samples are selected randomly in
#' each iteration.
#'
#' **Sampling Process**:
#' - For `"window"`: A sequential window of `size.train` is assigned for training, wrapping around if needed.
#' - For `"random"`: A random subset of size `size.train` is selected as the training set in each iteration.
#'
#' @return A binary matrix with dimensions `(number of samples) x (iteration)`,
#'   indicating training (`1`) and validation (`0`) samples for each iteration.
#'
#' @examples
#' \dontrun{
#' vect <- matrix(rnorm(100), nrow = 10)  # Example data
#' sampling_matrix <- modelSampling(vect, sampling.type = "random", training.pc = 80, iteration = 5)
#' print(sampling_matrix)
#' }
#'
#' @export
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


#' Generate a Summary Table of Best Models
#'
#' This function creates a data frame summarizing the best models from a list,
#' displaying model names, AUC scores, and model gene sets.
#'
#' @param model.sim A list of models, where each model is represented by a
#'   sublist containing elements with names that include `"score"` for the
#'   model's AUC score and `"model"` for the selected features.
#'
#' @details The function iterates through the provided models, extracting and
#' formatting the AUC score and sorted feature list for each. The resulting
#' table includes:
#'   - **N**: Model identifier (e.g., "N1" for the first model).
#'   - **AUC**: The AUC score for the model.
#'   - **MGS**: The model gene set, showing selected features in sorted order.
#'
#' **Structure**:
#' - The function processes up to 20 models, or the number of available models if fewer than 20.
#' - For each model, it retrieves the AUC score (by searching for `"score"` in the model name) and the model gene set (by searching for `"model"`).
#'
#' @return A data frame with columns:
#'   - `N`: Identifier for each model (e.g., "N1", "N2").
#'   - `AUC`: The AUC score for each model, rounded to 4 decimal places.
#'   - `MGS`: A string representing the sorted model gene set for each model.
#'
#' @examples
#' \dontrun{
#' model.sim <- list(
#'   model_N1 = list(score = 0.85, model = c("geneA", "geneB", "geneC")),
#'   model_N2 = list(score = 0.90, model = c("geneC", "geneD", "geneE"))
#' )
#' best_models_table <- tableBestModels(model.sim)
#' print(best_models_table)
#' }
#'
#' @export
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

