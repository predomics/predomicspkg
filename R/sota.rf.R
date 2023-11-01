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


#' sota.rf: launching Random Forest classifier
#'
#' @title sota.rf
#' @importFrom randomForest randomForest
#' @description sota.svm is a wrapper that executes svm using the same framework as for the predomics package.
#' @param sparsity: number of features in a given model. This is a vector with multiple lengths.
#' @param objective: prediction mode (default: auc)
#' @param max.nb.features: create the glmnet object using only the top most significant features (default:1000)
#' @param language is the language that is used by the different algorithms {bin, bininter, ter, terinter, ratio}, (default:"sota")
#' @param intercept: (Interceot for the a given model) (default:NULL)
#' @param evalToFit: Which model property will be used to select the best model among different k_sparsities (default: auc_)
#' @param k_penalty: Penalization of the fit by the k_sparsity (default: 0)
#' @param ntree: ??
#' @param mtry: Number of variables randomly sampled as candidates at each split. Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3)
#' @param replace: Should sampling of cases be done with or without replacement?
#' @param classwt: Priors of the classes. Need not add up to one. Ignored for regression.
#' @param sampsize: Size(s) of sample to draw. For classification, if sampsize is a vector of the length the number of strata, then sampling is stratified by strata, and the elements of sampsize indicate the numbers to be drawn from the strata.
#' @param nodesize: Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that the default values are different for classification (1) and regression (5).
#' @param maxnodes: Maximum number of terminal nodes trees in the forest can have. If not given, trees are grown to the maximum possible (subject to limits by nodesize). If set larger than maximum possible, a warning is issued.
#' @param importance: ??
#' @param localImp: ??
#' @param nPerm: ??
#' @param norm.votes: (??)
#' @param do.trace: ??
#' @param keep.forest: ??
#' @param cor.bias: ??
#' @param keep.inbag: ??
#' @param popSaveFile: (??)
#' @param seed: the seed to be used for reproductibility. If seed=NULL than it is not taken into account (default:NULL).
#' @param nCores: the number of CPUs to run the programm in parallel
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

# initialize the functions
sota.rf <- function(sparsity = c(1:30), # when sparsityis null it means that we can not fix it and the learner will use all the features.
                    objective = "auc",
                    max.nb.features = 1000,
                    intercept = "NULL",
                    language = "rf",
                    evalToFit = "auc_",
                    k_penalty=0,
                    ntree=500,
                    mtry=NULL,
                    replace=TRUE, 
                    classwt=NULL, 
                    #cutoff, strata,
                    sampsize = NULL,
                    nodesize = NULL,
                    maxnodes = NULL,
                    importance = FALSE, 
                    localImp = FALSE, 
                    nPerm = 1,
                    #proximity, 
                    #oob.prox=proximity,
                    norm.votes = TRUE, 
                    do.trace = FALSE,
                    keep.forest = TRUE, 
                    corr.bias = FALSE,
                    keep.inbag = FALSE,
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
  clf$learner                   <- "sota.rf" # name of the method
  clf$params                    <- list() # parameter list
  clf$params$sparsity           <- sparsity
  clf$params$objective          <- objective
  clf$params$max.nb.features    <- max.nb.features
  clf$params$language           <- language
  clf$params$intercept          <- intercept
  clf$params$parallel           <- FALSE
  clf$params$ntree              <- ntree  # Number of trees to grow.
  clf$params$mtry               <- mtry  # Number of variables randomly sampled as candidates at each split.
  clf$params$replace            <- replace # Should sampling of cases be done with or without replacement?
  clf$params$classwt            <- classwt # Priors of the classes
  #clf$params$cutoff            <- cutoff # A vector of length equal to number of classes. The ‘winning’ class for an observation is the one with the maximum ratio of proportion of votes to cutoff
  #clf$params$strata            <- strata # A (factor) variable that is used for stratified sampling
  clf$params$sampsize           <- sampsize # Size(s) of sample to draw
  clf$params$nodesize           <- nodesize # Minimum size of terminal nodes.
  clf$params$maxnodes           <- maxnodes # Maximum number of terminal nodes trees in the forest can have.
  clf$params$importance         <- importance # Should importance of predictors be assessed?
  clf$params$localImp           <- localImp # Should casewise importance measure be computed?
  clf$params$nPerm              <- nPerm # Number of times the OOB data are permuted per tree for assessing variable importance.
  #clf$params$proximity         <- proximity # Should proximity measure among the rows be calculated?
  #clf$params$oob.prox          <- oob.prox # Should proximity be calculated only on “out-of-bag” data?
  clf$params$norm.votes         <- norm.votes # If TRUE (default), the final result of votes are expressed as fractions. If FALSE, raw vote counts are returned
  clf$params$do.trace           <- do.trace # If set to TRUE, give a more verbose output as randomForest is run
  
  clf$params$plot               <- plot
  clf$params$verbose            <- verbose
  clf$params$warnings           <- warnings       # print out warnings
  clf$params$debug              <- debug          # print out logs.
  clf$params$print_ind_method   <- print_ind_method # method to print individual
  
  clf$params$evalToFit          <- evalToFit
  clf$params$k_penalty          <- k_penalty
  
  clf$params$keep.forest        <- keep.forest # If set to FALSE, the forest will not be retained in the output object.
  clf$params$corr.bias          <- corr.bias # perform bias correction for regression?
  clf$params$keep.inbag         <- keep.inbag # Should an n by ntree matrix be returned that keeps track of which samples are  “in-bag” in which trees
  clf$params$popSaveFile        <- popSaveFile
  
  clf$params$seed               <- seed
  clf$params$nCores             <- nCores         # parallel computing
  clf$params$parallel           <- nCores > 1     # parallel computing
  
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


sota.rf_fit <- function(X, y, clf) {
  # transpose (in the predomics package we have observations in the columns)
  x <- as.matrix(t(X))
  # compute the feature correlation for feature selection
  #feature.cor <- filterfeaturesK(data = t(x), trait = y, k = nrow(X), sort = TRUE) # to avoid having to recompute this all the time
  
  # initialize some parameters
  if(is.null(clf$params$mtry))
  {
    clf$params$mtry = if(!is.null(y) && !is.factor(y))
      max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x)))
  }  
  
  if(is.null(clf$params$sampsize)) 
  { 
    clf$params$sampsize = if(clf$params$replace) nrow(x) else ceiling(.632*nrow(x)) 
  }
  
  if(is.null(clf$params$nodesize)) 
  { 
    clf$params$nodesize = if(!is.null(y) && !is.factor(y)) 5 else 1 
  }
  
  if(is.null(clf$params$sparsity)) 
  { 
    clf$params$sparsity <- nrow(X)
  }
  
  model_collection <- list()
  
  for(i in 1:length(clf$params$sparsity))  # sparsity is = k, i.e. the number of features in a model
  {
    k <- clf$params$sparsity[i]
    if(clf$params$verbose) cat("... ... Resolving problem with\t", k, "\tvariables ...\n")
    
    selected.features <- rownames(clf$feature.cor)[order(clf$feature.cor$p)][1:k]
    x.reduced <- x[,selected.features]
    
    ### Launch randomForest with the parameters from the clf
    set.seed(clf$params$seed)
    rf <- randomForest(x = x.reduced, y=y,
                       ntree = clf$params$ntree,
                       mtry = clf$params$mtry, 
                       replace = clf$params$replace, 
                       classwt = clf$params$classwt,  
                       maxnodes = clf$params$maxnodes,
                       importance = clf$params$importance, 
                       localImp = clf$params$localImp,  
                       nPerm = nPerm,
                       norm.votes = clf$params$norm.votes, 
                       do.trace = clf$params$do.trace, 
                       keep.forest = clf$params$keep.forest, 
                       corr.bias = clf$params$corr.bias, 
                       keep.inbag = clf$params$keep.inbag
    )
    
    
    
    # build the model for this given sparsity
    mod.res                  <- list() # to send the results
    # Fill the object up with the rest of the values
    mod.res$learner          <- clf$learner
    mod.res$language         <- clf$params$language
    mod.res$names_           <- selected.features
    # match the index
    mod.res$indices_         <- match(selected.features, rownames(X))
    mod.res$eval.sparsity    <- length(unique(mod.res$names_))
    # add the objective in the model, needed for visualization
    mod.res$objective        <- clf$params$objective
    
    # Add the rf object
    mod.res$obj              <- rf
    mod.res$auc_             <- NA
    mod.res$cor_             <- NA
    mod.res$aic_             <- NA
    mod.res$score_           <- NA
    mod.res$fit_             <- NA
    mod.res$unpenalized_fit_ <- NA
    mod.res$intercept_       <- NA
    mod.res$sign_            <- NA
    # evaluate all
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
    
  }# end of sparsity loop
  
  # names(model_collection) <- paste("k",clf$params$sparsity, sep="_")
  # list(model_collection=model_collection, best.model.id=best.model)
  names(model_collection)    <- paste("k",clf$params$sparsity, sep="_")
  return(model_collection)
}


# NOTE to print a random forest object we used the code written by Rafael Zambrano here in stacks exchange https://stats.stackexchange.com/questions/41443/how-to-actually-plot-a-sample-tree-from-randomforestgettree

#**************************
#return the rules of a tree
#**************************
getConds <- function(tree)
{
  #store all conditions into a list
  conds <- list()
  #start by the terminal nodes and find previous conditions
  id.leafs <- which(tree$status==-1)
  j<-0
  for(i in id.leafs)
  {
    j <- j+1
    prevConds <- prevCond(tree,i)
    conds[[j]] <- prevConds$cond
    while(prevConds$id > 1)
    {
      prevConds <- prevCond(tree,prevConds$id)
      conds[[j]] <- paste(conds[[j]]," & ",prevConds$cond)
    }
    
    if(prevConds$id==1)
    {
      conds[[j]]<-paste(conds[[j]]," => ",tree$prediction[i])
    }
  } # end for
  return(conds)
}

#**************************
#find the previous conditions in the tree
#**************************
prevCond <- function(tree,i)
{
  if(i %in% tree$right_daughter)
  {
    id <- which(tree$right_daughter==i)
    cond <- paste(tree$split_var[id],">",tree$split_point[id])
  }
  
  if(i %in% tree$left_daughter)
  {
    id <- which(tree$left_daughter==i)
    cond <- paste(tree$split_var[id],"<",tree$split_point[id])
  }
  return(list(cond=cond,id=id))
}

#remove spaces in a word
collapse <- function(x)
{
  x <- sub(" ","_",x)
  return(x)
}

# written by Edi Prifti
printRFObject <-  function(obj)
{
  tree <- getTree(obj, k=1, labelVar=TRUE)

  # tree <- getTree(obj, k=2, labelVar = TRUE)
  # reprtree:::as.tree(gTree = tree, rforest = obj,max.depth = 2)
  #

  #rename the name of the column
  colnames(tree) <- sapply(colnames(tree), collapse)
  rules <- getConds(tree)
  print(unlist(rules))

}


tree_func <- function(final_model, 
                      tree_num, 
                      col.class = c("deepskyblue1", "firebrick1"), 
                      main = "",
                      node.text.color = "black") 
{
  # this function comes from Shirin Glander with some minor modifications on the visu.
  # https://shiring.github.io/machine_learning/2017/03/16/rf_plot_ggraph
  
  require(dplyr)
  require(ggraph)
  require(igraph)
  
  # get tree by index
  tree <- randomForest::getTree(final_model, 
                                k = tree_num, 
                                labelVar = TRUE) %>%
    tibble::rownames_to_column() %>%
    # make leaf split points to NA, so the 0s won't get plotted
    dplyr::mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))
  
  # prepare data frame for graph
  graph_frame <- data.frame(from = rep(tree$rowname, 2),
                            to = c(tree$`left daughter`, tree$`right daughter`))
  
  # convert to graph and delete the last node that we don't want to plot
  graph <- graph_from_data_frame(graph_frame) %>%
    delete_vertices("0")
  
  # set node labels
  V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
  V(graph)$leaf_label <- as.character(tree$prediction)
  V(graph)$split <- as.character(round(tree$`split point`, digits = 2))
  
  # plot
  plot <- ggraph(graph, 'dendrogram') + 
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    scale_fill_manual("", values = col.class) +
    geom_node_text(aes(label = node_label), na.rm = TRUE, repel = TRUE, colour = node.text.color) +
    geom_node_label(aes(label = split), vjust = 2.5, na.rm = TRUE, fill = "white") +
    geom_node_label(aes(label = leaf_label, fill = leaf_label), na.rm = TRUE, 
                    repel = TRUE, colour = "white", fontface = "bold", show.legend = FALSE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18))
  
  return(plot)
}

#tree_func(final_model = obj, 100)

