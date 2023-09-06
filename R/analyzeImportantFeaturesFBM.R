#' Visualize a summary of an experiment/set of experiments
#'
#' @description Visualization of 4 panels corresponding to feature prevalence in FBM, feature importance, feature prevalence in groups, effect sizes of feature abundances vs y-variable (cliff's delta for binary y; spearman rho for continuous y)
#' Can be applied to single classification task or to multiple classification tasks carried out on the same X-y dataset
#' @param clf_res The result of an experiment or multiple experiments (list of experiments)
#' @param X The feature table used as input of fit function behind experiments in clf_res
#' @param y The target class (binary/continuous)
#' @param makeplot  make a pdf file with the resulting plots (default:TRUE)
#' @param saveplotobj make a .Rda file with a list of the individual plots (default:TRUE)
#' @param name the suffix of the pdf file (default:"")
#' @param verbose print out informaiton
#' @param pdf.dims dimensions of the pdf object (default: c(w = 25, h = 20))
#' @param filter.cv.prev keep only features found in at least (default: 0.25, i.e 25 percent) of the cross validation experiments 
#' @param nb.top.features the maximum number (default: 100) of most important features to be shown.
#' If the number of features in FBM < nb.top.features, the number of features in FBM will be shown instead
#' @param scaled.importance the scaled importance is the importance multipied by the prevalence in the folds. If (default = TRUE) this will be used, the mean mda 
#' will be scaled by the prevalence of the feature in the folds and ordered subsequently 
#' @param k_penalty the sparsity penalty needed to select the best models of the population (default:0.75/100).
#' @param k_max select the best population below a given threshold. If (default:0) no selection is performed.
#'
#' @return plots if makeplot is FALSE; plot.list list object saved locally with individual plots (including source data) if saveplotobj
#' @export
#'
#' @examples
analyzeImportanceFeaturesFBM <- function(clf_res, X, y, makeplot = TRUE, saveplotobj=TRUE, name = "", 
                                         verbose = TRUE, pdf.dims = c(width = 25, height = 20), 
                                         filter.cv.prev = 0.25, 
                                         nb.top.features = 100, 
                                         scaled.importance = FALSE, 
                                         k_penalty = 0.75/100,
                                         k_max = 0)
{
  multiple.experiments <- FALSE
  mode <- NULL
  
  if(!isExperiment(clf_res))
  {
    if(any(!unlist(lapply(clf_res, isExperiment))))
    {
      stop("analyzeLearningFeatures: please provide a valid experiment results!")  
    }
    multiple.experiments <- TRUE
    
    if(clf_res[[1]]$classifier$params$objective == "cor")
    {
      mode <- "regression"
    }else
    {
      mode <- "classification"
    }
    
  }else # if not multiple experiments
  {
    if(clf_res$classifier$params$objective == "cor")
    {
      mode <- "regression"
    }else
    {
      mode <- "classification"
    }
  }
  
  if(!is.null(mode)) 
  {
    cat(paste("... Estimating mode: ", mode,"\n"))
  }else
  {
    stop("analyzeImportanceFeatures: mode not founding stopping ...")
  }
  
  if(!multiple.experiments)
  {
    # get the final population
    pop <- modelCollectionToPopulation(clf_res$classifier$models)
    if(verbose) print(paste("There are",length(pop), "models in this population"))
    # select the best population
    pop <- selectBestPopulation(pop = pop, score = clf_res$classifier$params$evalToFit, p = 0.05, k_penalty = k_penalty, k_max = k_max)
    if(verbose) print(paste("There are",length(pop), "models in this population after selection of the best"))
    
    if(length(pop) == 1)
    {
      print("analyzeImportanceFeatures: only one model after filtering. Plot can not be built... returing empty handed.")
      return(NULL)
    }
    
    # get the population information into a dataframe
    pop.df <- populationToDataFrame(pop = pop)
    # get the feature to model dataframe
    pop.noz <- listOfModelsToDenseCoefMatrix(clf = clf_res$classifier, X = X, y = y, list.models = pop)
    if(verbose) print(paste("Pop noz object is created with", nrow(pop.noz), "features and", ncol(pop.noz), "models"))
    
    # make the feature annots
    fa <- makeFeatureAnnot(pop = pop, X = X, y = y, clf = clf_res$classifier)
    pop.noz <- fa$pop.noz
    
    # get the feature importance information if it exists; do it for all features in X (subsequent filtering from FBM)
    lr <- list(clf_res)
    names(lr) <- paste(clf_res$classifier$learner, clf_res$classifier$params$language, sep=".")
    
    feat.import <- mergeMeltImportanceCV(list.results = lr, 
                                         filter.cv.prev = filter.cv.prev, 
                                         min.kfold.nb = FALSE, 
                                         learner.grep.pattern = "*", 
                                         # nb.top.features = nb.top.features, 
                                         nb.top.features = nrow(X),
                                         feature.selection = NULL,
                                         scaled.importance = scaled.importance,
                                         make.plot = TRUE)
    if(is.null(feat.import))
    {
      print("analyzeImportanceFeatures: no feature importance data found... returning empty handed.")
      return(NULL)
    }
    
    # Process the feat.import$summary to wide format (will be used for filtering features to plot)
    feat.import.dcast <- data.table::dcast(data = feat.import$summary, formula = feature~method, value.var="value")
    rownames(feat.import.dcast) <- feat.import.dcast$feature ; feat.import.dcast <- feat.import.dcast[,-1,drop=FALSE]
    feat.import.dcast$meanval <- rowMeans(feat.import.dcast, na.rm = TRUE)
    feat.import.dcast$feature <- factor(rownames(feat.import.dcast), levels=rownames(feat.import.dcast)[order(feat.import.dcast$meanval, decreasing = TRUE)])
    
    #Do the selection of features from feat.import.dcast
    features.import <- feat.import.dcast[rownames(pop.noz),]
    features.import$feature <- droplevels(features.import$feature)
    if(nrow(features.import)>nb.top.features) #if the number of features in FBM higher than the nb.top.features, select the top nb.top.features in FBM based on feature importance
    {
      print(paste(nrow(features.import), " features in FBM vs ", nb.top.features, " nb.top.features; select the ", nb.top.features, " top features with higher feature importance", sep = ""))
      features.import.vec <- levels(features.import$feature)[1:nb.top.features]
    } else #else, keep all features in FBM
    {
      print(paste(nrow(features.import), " features in FBM vs ", nb.top.features, " nb.top.features; reporting all features in FBM", sep = ""))
      features.import.vec <- levels(features.import$feature)
    }
    
    #Build the feature prevalence plot
    g1.df <- data.frame(pop.noz)
    g1.df$feature <- rownames(g1.df) ; g1.df <- melt(g1.df)
    g1.df$learner <- unlist(lapply(strsplit(as.character(g1.df$variable), split="_"), function(x){x[1]}))
    #add the language from clf object
    g1.df$learner <- paste(g1.df$learner, clf_res$classifier$params$language, sep = ".")
    g1.df$model <- as.character(unlist(lapply(strsplit(as.character(g1.df$variable), split="_"), function(x){x[2]})))
    #Subset g1.df to target features
    g1.df <- g1.df[g1.df$feature %in% features.import.vec,]
    g1.df$feature <- factor(g1.df$feature, levels=rev(features.import.vec))
    g1.df$value <- factor(g1.df$value, levels = c(-1,0,1))
    g1.df$value <- droplevels(g1.df$value)
    g1.df.cols <- c("deepskyblue1","white","firebrick1")
    names(g1.df.cols) <- c("-1","0","1")
    g1.df.cols <- g1.df.cols[levels(g1.df$value)]
    g1 <- ggplot(g1.df, aes(x=model, y=feature, fill=value)) + 
      geom_tile(colour=NA) + 
      facet_grid(.~learner, scales = "free_x", space = "free_x") + 
      scale_fill_manual(values = g1.df.cols) + 
      xlab("FBM") + 
      theme_classic() + 
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    
    # get the importance graphic
    g6.df <- feat.import$summary
    g6.df <- g6.df[g6.df$feature %in% features.import.vec,]
    g6.df$feature <- factor(g6.df$feature, levels=rev(features.import.vec))
    g6.df$sign <- factor(g6.df$sign, levels=c("-1","0","1"))
    g6.df$sign <- droplevels(g6.df$sign)
    g6.col <- c("deepskyblue1","white","firebrick1")
    names(g6.col) <- c("-1","0","1")
    g6.col <- g6.col[levels(g6.df$sign)]
    g6 <- ggplot(g6.df, aes(x=feature, y=value)) + 
      # the prevalence data on the bottom
      # geom_bar(data = fprev.melt, stat="identity", position="identity", aes(fill = "1")) +
      geom_hline(yintercept = min(0, g6.df$value, na.rm = TRUE), col="gray") +
      # ylab("Feature importance & prevalence (CV)") +
      ylab("Feature importance") +
      xlab("") +
      facet_grid(.~method, scales = "free") + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme_bw() +
      coord_flip() +
      # scale_fill_manual("Dataset", values = c("gray90","gray90")) +
      scale_color_manual("Dataset", values = g6.col) +
      # the importance data 
      geom_errorbar(aes(ymin = value - se, ymax = value + se, color = sign), width=.1, position=position_dodge(0.3)) +
      #geom_line(aes(group = feature, color = sign), position=pd) +
      geom_point(position = position_dodge(0.3), size=2, shape=19, aes(color = sign)) + # 21 is filled circle
      guides(colour = "none", fill = "none")
    g6 <- g6 + theme(axis.text.y = element_blank(),
                     strip.background = element_rect(fill = NA))
    
    # get the prevalence graphic
    if(mode == "regression")
    {
      g7 <- plotPrevalence(features = features.import.vec, X, y = NULL)
    }else
    {
      g7 <- plotPrevalence(features = features.import.vec, X, y)
    }
    g7 <- g7 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    # get the effect size graphic
    g8.df <- computeEffectSizes(X = X, y = y, mode = mode)
    #subset g8 to features.import
    g8.df <- g8.df[g8.df$feature %in% features.import.vec,]
    g8.df$feature <- factor(g8.df$feature, levels=rev(features.import.vec))
    #add FDR adjustment + labels
    if(mode=="classification")
    {
      g8.df$fdr <- stats::p.adjust(g8.df$pval.wilcox, method = "BH")
      g8.df$label <- ifelse(g8.df$fdr<0.05,"**", ifelse(g8.df$pval.wilcox<0.05,"*",""))
      g8 <- ggplot(g8.df, aes(x=feature, y=cdelta, fill=cdelta)) + 
        geom_bar(stat = "identity") + 
        geom_text(aes(label=label)) + 
        scale_fill_gradient(low = "deepskyblue1", high = "firebrick1", limits=c(-1,1)) + 
        ylab("Cliff's delta 1 vs. -1") + 
        ggtitle("") + 
        coord_flip() + 
        theme_bw() + 
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")
    } else if(mode=="regression")
    {
      g8.df$fdr <- stats::p.adjust(g8.df$pval, method = "BH")
      g8.df$label <- ifelse(g8.df$fdr<0.05,"**", ifelse(g8.df$pval<0.05,"*",""))
      g8 <- ggplot(g8.df, aes(x=feature, y=rho, fill=rho)) + 
        geom_bar(stat = "identity") + 
        geom_text(aes(label=label)) + 
        scale_fill_gradient(low = "deepskyblue1", high = "firebrick1", limits=c(-1,1)) + 
        ylab("Spearman Rho vs. y var") + 
        ggtitle("") + 
        coord_flip() + 
        theme_bw() + 
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")
    }
    
    ##Build a a list with plots to visualize
    plot.list <- list("featprevFBM"=g1,
                      "featImp"=g6,
                      "effectSizes"=g8,
                      "featPrevGroups"=g7)
    #Save plotlist object if specified
    if(saveplotobj)
    {
      save(plot.list, file=paste("population features",name,".Rda", sep=""))
    }
    #plotting statements
    if(makeplot)
    {
      if(verbose) print(paste("Making plots in a dedicated pdf"))
      #Set the layout for patchwork plotting
      layout <- "
      AAABBCD
      "
      #Set the file name
      fname <- paste("population features",name,".pdf", sep="")
      #Build the patchwork plot
      fname.plot <- patchwork::wrap_plots(plot.list, design = layout)
      #Save the plot with cowplot::ggsave2
      cowplot::ggsave2(filename = fname, 
                       plot = fname.plot, width = as.numeric(pdf.dims[1]), height = as.numeric(pdf.dims[2]))
      #Return the plot
      return(patchwork::wrap_plots(plot.list, design = layout))
    }else
    {
      layout <- "
      AAABBCD
      "
      return(patchwork::wrap_plots(plot.list, design = layout))
      # return(plot.list)
    }
  }else # multiple experiments
  {
    lr <- clf_res
    feat.import <- mergeMeltImportanceCV(list.results = lr, 
                                         filter.cv.prev = filter.cv.prev, 
                                         min.kfold.nb = FALSE, 
                                         learner.grep.pattern = "*", 
                                         # nb.top.features = nb.top.features, 
                                         nb.top.features = nrow(X), 
                                         make.plot = TRUE)
    
    if(is.null(feat.import))
    {
      warning("analyzeImportantFeatures: These learners have no feature importance as implemented in the BTR ones. Sending an empty plot.")
      return(ggplot() + theme_void())
    }
    
    # Process the feat.import$summary to wide format (will be used for filtering features to plot)
    feat.import.dcast <- data.table::dcast(data = feat.import$summary, formula = feature~method, value.var="value")
    rownames(feat.import.dcast) <- feat.import.dcast$feature ; feat.import.dcast <- feat.import.dcast[,-1]
    feat.import.dcast$meanval <- rowMeans(feat.import.dcast, na.rm = TRUE)
    feat.import.dcast$feature <- factor(rownames(feat.import.dcast), levels=rownames(feat.import.dcast)[order(feat.import.dcast$meanval, decreasing = TRUE)])
    
    # # the most important features along with the order
    # features.import <- rev(levels(feat.import$summary$feature))
    
    ## get the FBM feature prevalence graphic
    #Get the population of best models
    pop.list <- lapply(clf_res, function(x){modelCollectionToPopulation(x[["classifier"]][["models"]])})
    #Select the best populations (FBM)
    pop.list.fbm <- list()
    for(i in names(pop.list))
    {
      pop.list.fbm[[i]] <- selectBestPopulation(pop = pop.list[[i]], 
                                                score = clf_res[[i]][["classifier"]][["params"]][["evalToFit"]],
                                                p = 0.05, 
                                                k_penalty = k_penalty, 
                                                k_max = k_max)
    }
    # get the population information into a dataframe
    pop.list.fbm.df <- lapply(pop.list.fbm, function(x){populationToDataFrame(pop = x)})
    # get the feature to model dataframe
    pop.list.fbm.noz <- list()
    for(i in names(pop.list.fbm))
    {
      pop.list.fbm.noz[[i]] <- listOfModelsToDenseCoefMatrix(clf = clf_res[[i]][["classifier"]], X = X, y = y, list.models = pop.list.fbm[[i]])
    }
    
    # make the feature annots
    pop.list.fbm.fa <- list()
    for(i in names(pop.list.fbm))
    {
      pop.list.fbm.fa[[i]] <- makeFeatureAnnot(pop = pop.list.fbm[[i]], X = X, y = y, clf = clf_res[[i]][["classifier"]])
    }
    
    #merge both pop.nozdfs
    pop.list.fbm.fa.popnoz.df <- list()
    for(i in names(pop.list.fbm.fa))
    {
      idf <- data.frame(pop.list.fbm.fa[[i]][["pop.noz"]])
      idf$feature <- rownames(idf) ; idf <- melt(idf)
      idf$learner <- unlist(lapply(strsplit(as.character(idf$variable), split="_"), function(x){x[1]}))
      idf$model <- as.character(unlist(lapply(strsplit(as.character(idf$variable), split="_"), function(x){x[2]})))
      idf$source <- paste(idf$learner,i, sep = ".")
      pop.list.fbm.fa.popnoz.df[[i]] <- idf
    }
    pop.list.fbm.fa.popnoz.df <- do.call("rbind", pop.list.fbm.fa.popnoz.df)
    pop.list.fbm.fa.popnoz.df$value <- factor(pop.list.fbm.fa.popnoz.df$value, levels=c(-1,0,1))
    pop.list.fbm.fa.popnoz.df$value <- droplevels(pop.list.fbm.fa.popnoz.df$value)
    
    #Do the selection of features from feat.import.dcast
    features.import <- feat.import.dcast[unique(pop.list.fbm.fa.popnoz.df$feature),]
    features.import$feature <- droplevels(features.import$feature)
    if(nrow(features.import)>nb.top.features)
    {
      print(paste(nrow(features.import), " features in FBM from ", length(clf_res), " experiments vs ", nb.top.features, " nb.top.features; select the ", nb.top.features, " top features with higher mean feature importance across experiments", sep = ""))
      features.import.vec <- levels(features.import$feature)[1:nb.top.features]
    } else
    {
      print(paste(nrow(features.import), " features in FBM from ", length(clf_res), "experiments vs ", nb.top.features, " nb.top.features; reporting all features in FBM", sep = ""))
      features.import.vec <- levels(features.import$feature)
    }
    
    #Build the feature prevalence plot
    g1.df <- pop.list.fbm.fa.popnoz.df[pop.list.fbm.fa.popnoz.df$feature %in% features.import.vec,]
    g1.df$feature <- factor(g1.df$feature, levels=rev(features.import.vec))
    g1.col <- c("deepskyblue1","white","firebrick1")
    names(g1.col) <- c("-1","0","1")
    g1.col <- g1.col[levels(g1.df$value)]
    g1 <- ggplot(g1.df, aes(x=model, y=feature, fill=value)) + 
      geom_tile(colour=NA) + 
      facet_grid(.~source, scales = "free_x", space = "free_x") + 
      scale_fill_manual(values = g1.col) + 
      xlab("FBM") + 
      theme_classic() + 
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    
    # get the importance graphic
    g6.df <- feat.import$summary
    g6.df <- g6.df[g6.df$feature %in% features.import.vec,]
    g6.df$feature <- factor(g6.df$feature, levels=rev(features.import.vec))
    g6.df$sign <- factor(g6.df$sign, levels=c("-1","0","1"))
    g6.df$sign <- droplevels(g6.df$sign)
    #Add the learner from the clf list
    g6.df.learner <- data.frame(sapply(clf_res, function(x){x[["classifier"]][["learner"]]}))
    colnames(g6.df.learner) <- "learner"
    g6.df <- merge(g6.df, g6.df.learner, by.x="method", by.y=0, all.x=TRUE)
    g6.df$learner <- paste(g6.df$learner, g6.df$method, sep = ".")
    g6.col <- c("deepskyblue1","white","firebrick1")
    names(g6.col) <- c("-1","0","1")
    g6.col <- g6.col[levels(g6.df$sign)]
    g6 <- ggplot(g6.df, aes(x=feature, y=value)) + 
      # the prevalence data on the bottom
      # geom_bar(data = fprev.melt, stat="identity", position="identity", aes(fill = "1")) +
      geom_hline(yintercept = min(0, g6.df$value, na.rm = TRUE), col="gray") +
      # ylab("Feature importance & prevalence (CV)") +
      ylab("Feature importance") +
      xlab("") +
      facet_grid(.~learner, scales = "free") + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme_bw() +
      coord_flip() +
      # scale_fill_manual("Dataset", values = c("gray90","gray90")) +
      scale_color_manual("Dataset", values = g6.col) +
      # the importance data 
      geom_errorbar(aes(ymin = value - se, ymax = value + se, color = sign), width=.1, position=position_dodge(0.3)) +
      #geom_line(aes(group = feature, color = sign), position=pd) +
      geom_point(position = position_dodge(0.3), size=2, shape=19, aes(color = sign)) + # 21 is filled circle
      guides(colour = "none", fill = "none")
    g6 <- g6 + theme(axis.text.y = element_blank(), 
                     strip.background = element_rect(fill = NA))
    
    # get the prevalence graphic
    if(mode == "regression")
    {
      g7 <- plotPrevalence(features = features.import.vec, X, y = NULL)
    }else
    {
      g7 <- plotPrevalence(features = features.import.vec, X, y)
    }
    g7 <- g7 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    # get the effect size graphic
    g8.df <- computeEffectSizes(X = X, y = y, mode = mode)
    #subset g8 to features.import
    g8.df <- g8.df[g8.df$feature %in% features.import.vec,]
    g8.df$feature <- factor(g8.df$feature, levels=rev(features.import.vec))
    #add FDR adjustment + labels
    if(mode=="classification")
    {
      g8.df$fdr <- stats::p.adjust(g8.df$pval.wilcox, method = "BH")
      g8.df$label <- ifelse(g8.df$fdr<0.05,"**", ifelse(g8.df$pval.wilcox<0.05,"*",""))
      g8 <- ggplot(g8.df, aes(x=feature, y=cdelta, fill=cdelta)) + 
        geom_bar(stat = "identity") + 
        geom_text(aes(label=label)) + 
        scale_fill_gradient(low = "deepskyblue1", high = "firebrick1", limits=c(-1,1)) + 
        ylab("Cliff's delta 1 vs. -1") + 
        ggtitle("") + 
        coord_flip() + 
        theme_bw() + 
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")
    } else if(mode=="regression")
    {
      g8.df$fdr <- stats::p.adjust(g8.df$pval, method = "BH")
      g8.df$label <- ifelse(g8.df$fdr<0.05,"**", ifelse(g8.df$pval<0.05,"*",""))
      g8 <- ggplot(g8.df, aes(x=feature, y=rho, fill=rho)) + 
        geom_bar(stat = "identity") + 
        geom_text(aes(label=label)) + 
        scale_fill_gradient(low = "deepskyblue1", high = "firebrick1", limits=c(-1,1)) + 
        ylab("Spearman Rho vs. y var") + 
        ggtitle("") + 
        coord_flip() + 
        theme_bw() + 
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")
    }
    
    ##Build a a list with plots to visualize
    plot.list <- list("featprevFBM"=g1,
                      "featImp"=g6,
                      "effectSizes"=g8,
                      "featPrevGroups"=g7)
    #Save the plotobject if specified
    if(saveplotobj)
    {
      save(plot.list, file=paste("population features",name,".Rda", sep=""))
    }
    #Plotting
    if(makeplot)
    {
      if(verbose) print(paste("Making plots in a dedicated pdf"))
      #Set the layout for patchwork plotting
      layout <- "
      AAABBCD
      "
      #Set the file name
      fname <- paste("population features",name,".pdf", sep="")
      #Build the patchwork plot
      fname.plot <- patchwork::wrap_plots(plot.list, design = layout)
      #Save the plot with cowplot::ggsave2
      cowplot::ggsave2(filename = fname, 
                       plot = fname.plot, width = as.numeric(pdf.dims[1]), height = as.numeric(pdf.dims[2]))
      #Return the plot
      return(patchwork::wrap_plots(plot.list, design = layout))
    } else
    {
      layout <- "
      AAABBCD
      "
      return(patchwork::wrap_plots(plot.list, design = layout))
    }
  } # end multiple experiments
}


#' Compute effect sizes for features in binary classification/regression tasks
#'
#' @param X The X matrix (rows=features; columns=samples)
#' @param y The y vector of sample class (-1,1 in binary classification; continuous variable in regression)
#' @param mode classification or regression
#'
#' @description In binary classification tasks, compute the cliff's delta effect sizes btw groups (1 vs. -1) + pvalues from wilcoxon rank-sum tests; in regression tasks, compute spearman correlations (rho + pvalue) vs. continuous y variable
#' 
#' @return data frame of features, effect sizes (cliff's delta for binary classification; spearman rho for regression), and pvalues (wicoxon rank-sum test for binary classification task; spearman correlation for regression)
#' @export
#'
#' @examples
computeEffectSizes <- function(X,y,mode)
{
  if(mode=="classification")
  {
    #Transform y to factor
    y.factor <- factor(y, levels=c(1,-1))
    ##Check that y is binary (1, -1)
    if(length(which(is.na(y.factor)))>0)
    {
      stop("y vector contains other levels than 1,-1; please check if classification task")
    }
    #Build vector of sample classes
    names(y.factor) <- colnames(X)
    y.factor <- data.frame(y.factor)
    #melt the X
    X.melt <- as.data.frame(X)
    X.melt$feature <- rownames(X) ; X.melt <- melt(X.melt)
    X.melt <- merge(X.melt, y.factor, by.x="variable", by.y=0, all.x=TRUE)
    #Do the univariate tests + cliff delta calculations
    X.melt.tests <- list()
    for(v in unique(X.melt$feature))
    {
      v.cdelta <- effsize::cliff.delta(X.melt$value[X.melt$feature %in% v]~X.melt$y.factor[X.melt$feature %in% v])
      v.wilcox <- stats::wilcox.test(X.melt$value[X.melt$feature %in% v]~X.melt$y.factor[X.melt$feature %in% v])
      X.melt.tests[[v]] <- data.frame(feature=v,
                                      cdelta=v.cdelta$estimate,
                                      pval.wilcox=v.wilcox$p.value)
    }
    X.melt.tests <- do.call("rbind", X.melt.tests)
    rownames(X.melt.tests) <- seq(1,nrow(X.melt.tests))
    return(X.melt.tests)
  } else if(mode=="regression")
  {
    y.cont <- y
    names(y.cont) <- colnames(X)
    y.cont <- data.frame(y.cont)
    #melt the X
    X.melt <- as.data.frame(X)
    X.melt$feature <- rownames(X) ; X.melt <- melt(X.melt)
    X.melt <- merge(X.melt, y.cont, by.x="variable", by.y=0, all.x=TRUE)
    #Do the univariate tests + cliff delta calculations
    X.melt.tests <- list()
    for(v in unique(X.melt$feature))
    {
      v.spearman <- stats::cor.test(X.melt$value[X.melt$feature %in% v], X.melt$y.cont[X.melt$feature %in% v], method = "spearman")
      X.melt.tests[[v]] <- data.frame(feature=v,
                                      rho=v.spearman$estimate,
                                      pval=v.spearman$p.value)
    }
    X.melt.tests <- do.call("rbind", X.melt.tests)
    rownames(X.melt.tests) <- seq(1,nrow(X.melt.tests))
    return(X.melt.tests)
  } else
  {
    stop("mode not regression nor classification")
  }
}

#' Get objects needed for a merged visualization task combining different experiments from different datasets (different X and y)
#'
#' @description Here we get the 4 datasets from a given prediction experiment (clf object + X + y) needed for subsequent combination with other 
#' predition experiments for combined visualization (feature prevalence in FBM + feature importance + featureEffSizes + feature prevalence in groups)
#' @param clf_res The result of a single experiment
#' @param X The feature table used as input of fit function behind experiments in clf_res
#' @param y The target class (binary/continuous)
#' @param verbose print out informaiton
#' @param filter.cv.prev keep only features found in at least (default: 0.25, i.e 25 percent) of the cross validation experiments 
#' @param scaled.importance the scaled importance is the importance multipied by the prevalence in the folds. If (default = TRUE) this will be used, the mean mda 
#' will be scaled by the prevalence of the feature in the folds and ordered subsequently 
#' @param k_penalty the sparsity penalty needed to select the best models of the population (default:0.75/100).
#' @param k_max select the best population below a given threshold. If (default:0) no selection is performed.
#'
#' @return list of objects for subsequent combination
#' @export
#'
#' @examples
getImportanceFeaturesFBMobjects <- function(clf_res, 
                                            X, 
                                            y, 
                                            verbose = TRUE, 
                                            filter.cv.prev = 0.25, 
                                            scaled.importance = FALSE, 
                                            k_penalty = 0.75/100,
                                            k_max = 0)
{
  mode <- NULL
  if(!isExperiment(clf_res))
  {
    stop("analyzeLearningFeatures: please provide a valid experiment results!")  
  }
  if(clf_res$classifier$params$objective == "cor")
  {
    mode <- "regression"
  }else
  {
    mode <- "classification"
  }
  
  if(!is.null(mode)) 
  {
    cat(paste("... Estimating mode: ", mode,"\n"))
  }else
  {
    stop("analyzeImportanceFeatures: mode not founding stopping ...")
  }
  #############
  # get the final population
  #############
  pop <- modelCollectionToPopulation(clf_res$classifier$models)
  if(verbose) print(paste("There are",length(pop), "models in this population"))
  # select the best population
  pop <- selectBestPopulation(pop = pop, score = clf_res$classifier$params$evalToFit, p = 0.05, k_penalty = k_penalty, k_max = k_max)
  if(verbose) print(paste("There are",length(pop), "models in this population after selection of the best"))
  
  if(length(pop) == 1)
  {
    stop("analyzeImportanceFeatures: only one model after filtering. Plot can not be built... returing empty handed.")
  }
  
  # get the population information into a dataframe
  pop.df <- populationToDataFrame(pop = pop)
  # get the feature to model dataframe
  pop.noz <- listOfModelsToDenseCoefMatrix(clf = clf_res$classifier, X = X, y = y, list.models = pop)
  if(verbose) print(paste("Pop noz object is created with", nrow(pop.noz), "features and", ncol(pop.noz), "models"))
  # make the feature annots + get the data in plotting format
  fa <- makeFeatureAnnot(pop = pop, X = X, y = y, clf = clf_res$classifier)
  pop.noz <- fa$pop.noz
  pop.noz <- data.frame(pop.noz)
  pop.noz$feature <- rownames(pop.noz) ; pop.noz <- melt(pop.noz)
  pop.noz$learner <- unlist(lapply(strsplit(as.character(pop.noz$variable), split="_"), function(x){x[1]}))
  #add the language from clf object
  pop.noz$learner <- paste(pop.noz$learner, clf_res$classifier$params$language, sep = ".")
  pop.noz$model <- as.character(unlist(lapply(strsplit(as.character(pop.noz$variable), split="_"), function(x){x[2]})))
  pop.noz$value <- factor(pop.noz$value, levels = c(-1,0,1))
  pop.noz$value <- droplevels(pop.noz$value)
  
  ########
  # get the feature importance information if it exists; do it for all features in X (subsequent filtering from FBM)
  ########
  lr <- list(clf_res)
  names(lr) <- paste(clf_res$classifier$learner, clf_res$classifier$params$language, sep=".")
  feat.import <- mergeMeltImportanceCV(list.results = lr, 
                                       filter.cv.prev = filter.cv.prev, 
                                       min.kfold.nb = FALSE, 
                                       learner.grep.pattern = "*", 
                                       # nb.top.features = nb.top.features, 
                                       nb.top.features = nrow(X),
                                       feature.selection = NULL,
                                       scaled.importance = scaled.importance,
                                       make.plot = TRUE)
  if(is.null(feat.import))
  {
    stop("analyzeImportanceFeatures: no feature importance data found... returning empty handed.")
  }
  
  # get the prevalence graphic object
  if(mode == "regression")
  {
    featPrevPlot <- plotPrevalence(features = rownames(X), X, y = NULL)
  }else
  {
    featPrevPlot <- plotPrevalence(features = rownames(X), X, y)
  }
  # get the effect size graphic
  effSizes.df <- computeEffectSizes(X = X, y = y, mode = mode)
  
  ##Build a a list with objects to save
  outlist <- list("featprevFBM"=pop.noz,
                  "featImp"=feat.import$summary,
                  "effectSizes"=effSizes.df,
                  "featPrevGroups"=featPrevPlot$data)
  ## Add specific class label to outlist for subsequent checking in plotImportanceFeaturesFBMobjects
  class(outlist) <- "listFBMfeatures"
  return(outlist)
}


#' Visualize a list containing outouts of getImportanceFeaturesFBMobjects
#'
#' @description Here we combine the 4 datasets generated by getImportanceFeaturesFBMobjects function from different prediction experiments (clf object + X + y) ; designed to combine predomics results with different X,y source data
#' for unified visualization (feature prevalence in FBM, feature importance, feature effect size across groups, feature prevalence across groups)
#' @param FBMobjList List of outputs of getImportanceFeaturesFBMobjects function (1 list per experiment to combine)
#' @param verbose print out informaiton
#' @param nb.top.features features to retain for visualization (top features with highest mean feature importance across datasets)
#' @param makeplot  make a pdf file with the resulting plots (default:TRUE)
#'
#' @return Combined visualization of feature prevalence in FBM + feature importance + feature effect size across groups + feature prevalence across groups in different predomics prediction tasks
#' @export
#'
#' @examples
#' 

plotImportanceFeaturesFBMobjects <- function(FBMobjList, 
                                             verbose = TRUE, 
                                             nb.top.features=100,
                                             makeplot=TRUE)
{
  #check that elements in FBMobjList are of class listFBMfeatures
  FBMobjList.check <- sapply(FBMobjList, function(x){class(x)=="listFBMfeatures"})
  if(sum(FBMobjList.check)!=length(FBMobjList))
  {
    FBMobjList.check.wrong <- names(FBMobjList.check)[!FBMobjList.check]
    stop(print(paste(paste(FBMobjList.check.wrong, collapse = ","), "elements in FBMobjList not produced by getImportanceFeaturesFBMobjects (not of class listFBMfeatures)"))) 
  }
  #get the feature importance data + add dataset source
  FBMobjList.fimp <- lapply(FBMobjList, function(x){x[["featImp"]]})
  for(i in names(FBMobjList.fimp))
  {
    FBMobjList.fimp[[i]][,"dataset"] <- i
  }
  #Get the FBM dataframes + the feature importances of the corresponding features
  FBMobjList.fbm <- lapply(FBMobjList, function(x){x[["featprevFBM"]]})
  for(i in names(FBMobjList.fbm))
  {
    FBMobjList.fbm[[i]][,"dataset"] <- i
  }
  #get the unique features in FBMobjList.fbm
  FBMobjList.fbm_features <- unique(unlist(lapply(FBMobjList.fbm, function(x){x[,"feature"]})))
  
  #get the features to plot
  FBMobjList.fimp.df <- do.call("rbind", FBMobjList.fimp)
  FBMobjList.fimp.df.dcast <- dcast(data = FBMobjList.fimp.df, formula = feature~dataset, value.var="value")
  rownames(FBMobjList.fimp.df.dcast) <- FBMobjList.fimp.df.dcast[,1] ; FBMobjList.fimp.df.dcast <- FBMobjList.fimp.df.dcast[,-1,drop=FALSE]
  #compute the average of feat importance across datasets
  FBMobjList.fimp.df.dcast$avg <- rowSums(FBMobjList.fimp.df.dcast, na.rm = TRUE)/ncol(FBMobjList.fimp.df.dcast)
  tfeats <- rownames(FBMobjList.fimp.df.dcast)
  tfeats.fbm <- intersect(tfeats, FBMobjList.fbm_features)
  
  if(length(tfeats.fbm)>nb.top.features)
  {
    if(verbose)
    {
      print(paste(length(tfeats.fbm), " unique features in FBM of combined datasets with feature importance vs. ",nb.top.features, " nb.top.features; selecting ", nb.top.features, " features with highest mean feature importance across datasets"))
    }
    tfeats.fbm <- FBMobjList.fimp.df.dcast[tfeats.fbm,]
    tfeats.fbm <- rownames(tfeats.fbm)[order(tfeats.fbm$avg, decreasing = TRUE)]
    tfeats.fbm <- tfeats.fbm[1:nb.top.features]
  } else
  {
    if(verbose)
    {
      print(paste(length(tfeats.fbm), " unique features in FBM of combined datasets with fefature importance vs. ",nb.top.features, " nb.top.features; keeping all features in FBM of combined datasets"))
    }
    tfeats.fbm <- FBMobjList.fimp.df.dcast[tfeats.fbm,]
    tfeats.fbm <- rownames(tfeats.fbm)[order(tfeats.fbm$avg, decreasing = TRUE)]
  }
  ### Visu prevalence in FBM
  FBMobjList.fbm_plotdf <- lapply(FBMobjList.fbm, function(x){x[x[,"feature"] %in% tfeats.fbm,]})
  FBMobjList.fbm_plotdf <- do.call("rbind", FBMobjList.fbm_plotdf)
  FBMobjList.fbm_plotdf$feature <- factor(FBMobjList.fbm_plotdf$feature, levels = rev(tfeats.fbm))
  plot1.cols <- c("deepskyblue1","white","firebrick1")
  names(plot1.cols) <- c("-1","0","1")
  plot1.cols <- plot1.cols[levels(FBMobjList.fbm_plotdf$value)]
  plot1 <- ggplot(FBMobjList.fbm_plotdf, aes(x=model, y=feature, fill=value)) + 
    geom_tile(colour=NA) + 
    facet_grid(.~dataset, scales = "free_x", space = "free_x") + 
    scale_fill_manual(values = plot1.cols) + 
    xlab("FBM") + 
    theme_classic() + 
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  ### Visu the feature importance
  FBMobjList.fimp_plot <- lapply(FBMobjList.fimp, function(x){x[x[,"feature"] %in% tfeats.fbm,]})
  FBMobjList.fimp_plot <- do.call("rbind", FBMobjList.fimp_plot)
  FBMobjList.fimp_plot$feature <- factor(FBMobjList.fimp_plot$feature, levels=rev(tfeats.fbm))
  FBMobjList.fimp_plot$sign <- factor(FBMobjList.fimp_plot$sign, levels=c("-1","1"))
  plot2.col <- c("deepskyblue1","white","firebrick1")
  names(plot2.col) <- c("-1","0","1")
  plot2.col <- plot2.col[levels(FBMobjList.fimp_plot$sign)]
  plot2 <- ggplot(FBMobjList.fimp_plot, aes(x=feature, y=value)) + 
    # the prevalence data on the bottom
    # geom_bar(data = fprev.melt, stat="identity", position="identity", aes(fill = "1")) +
    geom_hline(yintercept = min(0, FBMobjList.fimp_plot$value, na.rm = TRUE), col="gray") +
    # ylab("Feature importance & prevalence (CV)") +
    ylab("Feature importance") +
    xlab("") +
    facet_grid(.~dataset) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme_bw() +
    coord_flip() +
    # scale_fill_manual("Dataset", values = c("gray90","gray90")) +
    scale_color_manual("Dataset", values = plot2.col) +
    # the importance data 
    geom_errorbar(aes(ymin = value - se, ymax = value + se, color = sign), width=.1, position=position_dodge(0.3)) +
    #geom_line(aes(group = feature, color = sign), position=pd) +
    geom_point(position = position_dodge(0.3), size=2, shape=19, aes(color = sign)) + # 21 is filled circle
    guides(colour = "none", fill = "none")
  plot2 <- plot2 + theme(axis.text.y = element_blank(),
                         strip.background = element_rect(fill = NA))
  
  ### Visu the feature effectsizes
  FBMobjList.effsize <- lapply(FBMobjList, function(x){x[["effectSizes"]]})
  #subset to tfeats
  FBMobjList.effsize <- lapply(FBMobjList.effsize, function(x){x[x[,"feature"] %in% tfeats.fbm,]})
  #fdr adjustment of pvalues
  for(i in names(FBMobjList.effsize))
  {
    idf <- FBMobjList.effsize[[i]]
    idf[,"fdr"] <- p.adjust(idf[,3], method = "BH")
    idf$dataset <- i
    FBMobjList.effsize[[i]] <- idf
  }
  FBMobjList.effsize <- do.call("rbind", FBMobjList.effsize)
  FBMobjList.effsize$label <- ifelse(FBMobjList.effsize[,4]<0.05,"**", ifelse(FBMobjList.effsize[,3]<0.05,"*",""))
  FBMobjList.effsize$feature <- factor(FBMobjList.effsize$feature, levels=rev(tfeats.fbm))
  plot3 <- ggplot(FBMobjList.effsize, aes(x=feature, y=cdelta, fill=cdelta)) + 
    geom_bar(stat = "identity") + 
    geom_text(aes(label=label)) + 
    scale_fill_gradient(low = "deepskyblue1", high = "firebrick1", limits=c(-1,1)) + 
    ylab("Cliff's delta 1 vs. -1") + 
    ggtitle("") + 
    facet_grid(.~dataset) + 
    coord_flip() + 
    theme_bw() + 
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill = NA),
          legend.position = "none")
  #Visu the feature prevalence across groups
  FBMobjList.fprev <- lapply(FBMobjList, function(x){x[["featPrevGroups"]]})
  FBMobjList.fprev <- lapply(FBMobjList.fprev, function(x){x[x[,"feature"] %in% tfeats.fbm,]})
  for(i in names(FBMobjList.fprev))
  {
    FBMobjList.fprev[[i]][,"dataset"] <- i
  }
  FBMobjList.fprev <- do.call("rbind", FBMobjList.fprev)
  FBMobjList.fprev$feature <- factor(FBMobjList.fprev$feature, levels = rev(tfeats.fbm))
  col.pt = c("deepskyblue4", "firebrick4") 
  col.bg = c("deepskyblue1", "firebrick1")
  plot4 <- ggplot(FBMobjList.fprev, aes(x=feature, y=prevalence, fill=group)) + 
    geom_bar(data=subset(FBMobjList.fprev, group %in% c("all")), stat="identity", position="identity") + 
    facet_grid(.~dataset) + 
    coord_flip() + 
    geom_point(data = subset(FBMobjList.fprev, group %in% c("-1", "1")), aes(x=feature, y=prevalence, color=group, shape=group)) + 
    scale_color_manual("Dataset", values = c("all"="gray90", "-1"=col.pt[1], "1"=col.pt[2])) +
    scale_fill_manual("Dataset", values = c("all"="gray90", "-1"=col.bg[1], "1"=col.bg[2])) +
    scale_shape_manual(values=c(25,24)) + 
    theme_bw() + 
    theme(legend.position="none", axis.text=element_text(size=9),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_rect(fill = NA)) + 
    ggtitle("") 
  ##Build a a list with plots to visualize
  plot.list <- list("featprevFBM"=plot1,
                    "featImp"=plot2,
                    "effectSizes"=plot3,
                    "featPrevGroups"=plot4)
  if(makeplot)
  {
    if(verbose) print(paste("Making plots in a dedicated pdf"))
    #Set the layout for patchwork plotting
    layout <- "
      AAABBCD
      "
    #Set the file name
    fname <- paste("population features multipleExperiments",name,".pdf", sep="")
    #Build the patchwork plot
    fname.plot <- patchwork::wrap_plots(plot.list, design = layout)
    #Save the plot with cowplot::ggsave2
    cowplot::ggsave2(filename = fname, 
                     plot = fname.plot, width = as.numeric(pdf.dims[1]), height = as.numeric(pdf.dims[2]))
    #Return the plot
    return(patchwork::wrap_plots(plot.list, design = layout))
  }else
  {
    layout <- "
      AAABBCD
      "
    return(patchwork::wrap_plots(plot.list, design = layout))
  }
}
