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
# @script: global.visu.R                                          
# @author: Edi Prifti
# @author: Lucas Robin
# @author: Yann Chevaleyre
# @author: Jean-Daniel Zucker
# @date: August 2016                                         
################################################################
# CONTENTS
# ========= PLOT EVALUATION RESULTS
# ========= MODEL VISUALISATION
# ========= POPULATION VISUALISATION
# ========= PRINT DIFFERENT OBJECTS
################################################################


################################################################
# PLOT EVALUATION RESULTS
################################################################

#' Plots a graph for a given score
#' 
#' @import RColorBrewer
#' @import ggplot2
#' @description plotComparativeEmpiricalScore plots a digested.results data object for a given score.
#' @param digested.results: a list of data.frames containing performance results from a lists of learners. This data object is returned by the function merge_digestScores()
#' @param ylim: y-axis zoom in the plot
#' @param score: default (auc_) score
#' @param main: name of the graphic
#' @return A ggplot graphs
#' @export
plotComparativeEmpiricalScore <- function(digested.results, ylim = c(0.5,1), score = "auc_", main = "")
{
  nbe <- length(table(digested.results$empirical[[score]]$method))# number of elements to compare

  if(!is.null(digested.results$colors))
  {
    colors <- digested.results$colors
  }else
  {
    colors <- c(brewer.pal(9,"Set1"),brewer.pal(9, "Set3"), brewer.pal(9, "Spectral"),brewer.pal(8, "Dark2"))[1:nbe]  
  }
  
  if(!is.null(digested.results$pch))
  {
    pch <- digested.results$pch
  }else
  {
    pch <- rep(19, nbe)  
  }
  
  pd <- position_dodge(0.3) # move them .05 to the left and right
  
  # EMPIRICAL SCORES 
  data <- digested.results$empirical[[score]]
  # g <- ggplot(na.omit(data), aes(x=variable, y=value, colour=method)) + 
  #   #geom_errorbar(aes(ymin=value, ymax=value), width=.1, position=pd) +
  #   geom_line(aes(group=method), position=pd) +
  #   geom_point(position=pd, size=2, shape=19) + # 21 is filled circle
  #   xlab("k sparse") +
  #   ylab(score) +
  #   ggtitle(paste("WD","empirical",score, sep="; ")) +
  #   expand_limits(y=ylim) +                        # Expand y range
  #   theme_bw() + scale_colour_manual(values=colors) + scale_fill_manual(values=colors) +
  #   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  #   theme(legend.justification=c(1,0), legend.position=c(1,0), legend.direction="vertical") # Position legend in bottom right
  
  maj.class <- NA
  
  # Majoritary class
  if(score == "auc_")
  {
    maj.class <- 0.5
  }
  
  if(score == "accuracy_")
  {
    maj.class <- unique(digested.results$maj.class)
  }
  
  
  # if no data exists
  if(is.null(data))
  {
    g <- ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(ylim) + 
      theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
      ylab(score) +
      xlab("k sparse") +
      ggtitle(paste("WD","empirical",score, main, sep="; ")) +
      geom_hline(aes(yintercept=1), lty=2, col="lightgray") + 
      geom_hline(aes(yintercept=0.5), lty=2, col="lightgray") + 
      theme(legend.position="bottom", legend.direction="horizontal") +
      guides(col = guide_legend(ncol = 3, byrow = TRUE))
    return(g)
  }
  
  # if not specified set it as the range
  if(is.null(ylim))
  {
    ylim <- range(data$value)
  }
  
  # if specified range is not outside the real values
  if(min(data$value, na.rm = TRUE) < min(ylim, na.rm = TRUE))
  {
    ylim[1] <- min(data$value, na.rm = TRUE)
  }
  
  if(max(data$value, na.rm = TRUE) > max(ylim, na.rm = TRUE))
  {
    ylim[2] <- max(data$value, na.rm = TRUE)
  }
  
  g <- ggplot(na.omit(data), aes(x=variable, y=value, colour=method, shape = method)) + 
    #geom_errorbar(aes(ymin=value, ymax=value), width=.1, position=pd) +
    geom_line(aes(group=method), position=pd) +
    geom_point(position=pd, size=2) + # 21 is filled circle
    xlab("k sparse") +
    ggtitle(paste("WD","empirical",score, main, sep="; ")) +
    geom_hline(aes(yintercept=1), lty=2, col="lightgray") + 
    geom_hline(aes(yintercept=0.5), lty=2, col="lightgray") + 
    ylab(score) +
    coord_cartesian(ylim = ylim) +
    #scale_y_continuous(limits = ylim) +
    scale_colour_manual(values=colors, name = "") + 
    scale_fill_manual(values=colors) + 
    scale_shape_manual(values=pch, name = "") +
    theme_bw() +
    #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.text.x=element_text(angle=90,hjust=1, vjust=0.5)) +
    theme(legend.position="bottom", legend.direction="horizontal") +
    guides(col = guide_legend(ncol = 3, byrow = TRUE))
  
  if(!is.na(maj.class))
  {
    g <-  g +  geom_hline(aes(yintercept=maj.class), lty=2, col="black")
  }
  
  return(g)
}


#' Plots a graph for a given score
#' 
#' @import RColorBrewer
#' @import ggplot2
#' @description plotComparativeCV plots a digested.results data object for a given score.
#' @param digested.results: a list of data.frames containing performance results from a lists of learners. This data object is returned by the function merge_digestScores()
#' @param ylim: y-axis zoom in the plot
#' @param generalization: when (default:TRUE) then the generalization score will be used
#' @param score: default (auc_) score for the cross-validation representation
#' @param ci: should the confidence intereval be plotted (default:TRUE)
#' @param main: name of the graphic
#' @return A ggplot graphs
#' @export
plotComparativeCV <- function(digested.results, ylim = c(0.5,1), generalization = TRUE, score="auc_", ci = TRUE, main = "")
{
  nbe <- length(table(digested.results$empirical[[score]]$method))# number of elements to compare
  
  if(!is.null(digested.results$colors))
  {
    colors <- digested.results$colors
  }else
  {
    colors <- c(brewer.pal(9,"Set1"),brewer.pal(9, "Set3"), brewer.pal(9, "Spectral"),brewer.pal(8, "Dark2"))[1:nbe]  
  }
  
  if(!is.null(digested.results$pch))
  {
    pch <- digested.results$pch
  }else
  {
    pch <- rep(19, nbe)  
  }
  
  pd <- position_dodge(0.3) # move them .05 to the left and right

  # EMPIRICAL SCORES 
  # TODO ORDER BY K not alphabetically
  
  if(generalization) 
  { # Generalization
    if(score=="auc_")
    {
      data <- digested.results$cv.generalization.auc$cv.se
    }else if(score=="accuracy_")
    {
      data <- digested.results$cv.generalization.acc$cv.se
    }else if(score=="recall_")
    {
      data <- digested.results$cv.generalization.rec$cv.se
    }else if(score=="precision_")
    {
      data <- digested.results$cv.generalization.prc$cv.se
    }else if(score=="f1_")
    {
      data <- digested.results$cv.generalization.f1s$cv.se
    }else if(score=="cor_")
    {
      data <- digested.results$cv.generalization.cor$cv.se
    } else
    {
      stop("plotComparativeCV: please enter an existing score!")
    }
    type <- "generalization"
  }else 
  { # empirical
    if(score=="auc_")
    {
      data <- digested.results$cv.empirical.auc$cv.se
    }else if(score=="accuracy_")
    {
      data <- digested.results$cv.empirical.acc$cv.se
    }else if(score=="recall_")
    {
      data <- digested.results$cv.empirical.rec$cv.se
    }else if(score=="precision_")
    {
      data <- digested.results$cv.empirical.prc$cv.se
    }else if(score=="f1_")
    {
      data <- digested.results$cv.empirical.f1s$cv.se
    }else if(score=="cor_")
    {
      data <- digested.results$cv.empirical.cor$cv.se
    }else
    {
      stop("plotComparativeCV: please enter an existing score!")
    }
    type <- "empirical"
  }
  
  maj.class <- NA
  
  # Majoritary class
  if(score == "auc_")
  {
    maj.class <- 0.5
  }
  
  if(score == "accuracy_")
  {
    maj.class <- unique(digested.results$maj.class)
  }
  
  # if no data exists
  if(is.null(data))
  {
    g <- ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(ylim) + 
      theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1,vjust=0.5)) +
      ylab(score) +
      xlab("k sparse") +
      ggtitle(paste("WD","empirical",score, main, sep="; ")) +
      geom_hline(aes(yintercept=1), lty=2, col="lightgray") + 
      geom_hline(aes(yintercept=0.5), lty=2, col="lightgray") + 
      theme(legend.position="bottom", legend.direction="horizontal") +
      guides(col = guide_legend(ncol = 3, byrow = TRUE))
  }else
  {
    data$shape = pch[data$method]
    
    if(ci)
    {
      g <- ggplot(na.omit(data), aes(x=variable, y=value, colour=method, shape = method)) + 
        geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.1, position=pd) +
        geom_line(aes(group=method), position=pd) +
        geom_point(position=pd, size=2) + 
        ggtitle(paste("CV", type, score, main, sep="; ")) +
        geom_hline(aes(yintercept=1), lty=2, col="lightgray") + 
        geom_hline(aes(yintercept=0.5), lty=2, col="lightgray") + 
        #geom_hline(aes(yintercept=maj.class), lty=2, col="black") + 
        ylab(score) +
        coord_cartesian(ylim = ylim) +
        #scale_y_continuous(limits = ylim) +
        scale_colour_manual(values=colors, name = "") + 
        scale_fill_manual(values=colors) + 
        theme_bw() +
        #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
        theme(legend.position="bottom", legend.direction="horizontal") +
        guides(col = guide_legend(ncol = 3, byrow = TRUE)) + 
        scale_shape_manual(values=pch, name = "")
    }else
    {
      g <- ggplot(na.omit(data), aes(x=variable, y=value, colour=method, shape = method)) + 
        #geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.1, position=pd) +
        geom_line(aes(group=method),position=pd) +
        geom_point(position=pd, size=2) + # 21 is filled circle
        ggtitle(paste("CV", type, score, main, sep="; ")) +
        geom_hline(aes(yintercept=1), lty=2, col="lightgray") + 
        geom_hline(aes(yintercept=0.5), lty=2, col="lightgray") + 
        #geom_hline(aes(yintercept=maj.class), lty=2, col="black") + 
        ylab(score) +
        coord_cartesian(ylim = ylim) +
        #scale_y_continuous(limits = ylim) +
        scale_colour_manual(values=colors, name = "") + 
        scale_fill_manual(values=colors) + 
        theme_bw() +
        #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
        theme(legend.position="bottom", legend.direction="horizontal") +
        guides(col = guide_legend(ncol = 3, byrow = TRUE)) +
        scale_shape_manual(values=pch, name = "")
    }
    
    if(!is.na(maj.class))
    {
      g <-  g +  geom_hline(aes(yintercept=maj.class), lty=2, col="black")
    }
  }
  
  return(g)
}


#' Plots a graph for a given score
#' 
#' @import RColorBrewer
#' @import ggplot2
#' @description plotComparativeCV plots a digested.results data object for a given score.
#' @param digested.results: a list of data.frames containing performance results from a lists of learners. This data object is returned by the function merge_digestScores()
#' @param ylim: y-axis zoom in the plot
#' @param score: default (auc_) score for the cross-validation representation
#' @param main: name of the graphic
#' @return A ggplot graphs
#' @export
plotComparativeBestCV <- function(digested.results, ylim = c(0.5,1), generalization = TRUE, score = "auc_", ci = TRUE, main = "")
{
  
  # global information for the plots
  nbe <- length(table(digested.results$empirical[[score]]$method))# number of elements to compare
  if(!is.null(digested.results$colors))
  {
    colors <- digested.results$colors
  }else
  {
    colors <- c(brewer.pal(9,"Set1"),brewer.pal(9, "Set3"), brewer.pal(9, "Spectral"),brewer.pal(8, "Dark2"))[1:nbe]  
  }
  
  if(!is.null(digested.results$pch))
  {
    pch <- digested.results$pch
  }else
  {
    pch <- rep(19, nbe)  
  }
  
  pd <- position_dodge(0.3) # move them .05 to the left and right
  
  if(generalization) 
  { # Generalization
    if(score=="auc_")
    {
      data <- digested.results$best.auc$summary.gen
      best <- digested.results$best.auc$split
    }else if(score=="accuracy_")
    {
      data <- digested.results$best.acc$summary.gen
      best <- digested.results$best.acc$split
    }else if(score=="recall_")
    {
      data <- digested.results$best.rec$summary.gen
      best <- digested.results$best.rec$split
    }else if(score=="precision_")
    {
      data <- digested.results$best.prc$summary.gen
      best <- digested.results$best.prc$split
    }else if(score=="f1_")
    {
      data <- digested.results$best.f1s$summary.gen
      best <- digested.results$best.f1s$split
    }else if(score=="cor_")
    {
      data <- digested.results$best.cor$summary.gen
      best <- digested.results$best.cor$split
    } else
    {
      stop("plotComparativeBestCV: please enter an existing score!")
    }
    type <- "generalization"
  }else 
  { # empirical
    if(score=="auc_")
    {
      data <- digested.results$best.auc$summary.emp
      best <- digested.results$best.auc$split
    }else if(score=="accuracy_")
    {
      data <- digested.results$best.acc$summary.emp
      best <- digested.results$best.acc$split
    }else if(score=="recall_")
    {
      data <- digested.results$best.rec$summary.emp
      best <- digested.results$best.rec$split
    }else if(score=="precision_")
    {
      data <- digested.results$best.prc$summary.emp
      best <- digested.results$best.prc$split
    }else if(score=="f1_")
    {
      data <- digested.results$best.f1s$summary.emp
      best <- digested.results$best.f1s$split
    }else if(score=="cor_")
    {
      data <- digested.results$best.cor$summary.emp
      best <- digested.results$best.cor$split
    }else
    {
      stop("plotComparativeBestCV: please enter an existing score!")
    }
    type <- "empirical"
  }
  
  if(is.null(data))
  {
    warning(paste("plotComparativeBestCV: no data for score",score))
    return(NULL)
  }
  
  data$colors <- colors
  
  # best k_sparsity (for the legend)
  k_sparsity <- unlist(lapply(best, function(x){return(unique(as.character(x[,"k_sparse"])))}))
  data$Lmethod <- paste(data$Lmethod, " (",k_sparsity[data$Lmethod],")", sep="")
  maj.class <- NA
  
  # Majoritary class
  if(score == "auc_")
  {
    maj.class <- 0.5
  }
  
  if(score == "accuracy_")
  {
    maj.class <- unique(digested.results$maj.class)
  }
  
  if(ci)
  {
    g <- ggplot(na.omit(data), aes(x=Lmethod, y=value, colour=Lmethod, shape=Lmethod)) +
      geom_point(position=pd, size=2) +
      geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.1, position=pd) +
      ggtitle(paste("CV", type, score, main, sep="; ")) +
      ylab(score) +
      coord_cartesian(ylim = ylim) +
      #scale_y_continuous(limits = ylim) +
      geom_hline(aes(yintercept=1), lty=2, col="lightgray") + 
      geom_hline(aes(yintercept=0.5), lty=2, col="lightgray") +
      scale_colour_manual(values = as.character(colors), name = "") + 
      scale_fill_manual(values = as.character(colors)) + 
      theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(legend.position="bottom", legend.direction="horizontal") +
      guides(col = guide_legend(ncol = 3, byrow = TRUE)) +
      scale_shape_manual(values=as.numeric(pch), name = "")
    
  }else
  {
    g <- ggplot(na.omit(data), aes(x=Lmethod, y=value, colour=Lmethod, shape = Lmethod)) +
      geom_point(position=pd, size=2) +
      #geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.1, position=pd) +
      ggtitle(paste("CV", type, score, main, sep="; ")) +
      ylab(score) +
      coord_cartesian(ylim = ylim) +
      #scale_y_continuous(limits = ylim) +
      geom_hline(aes(yintercept=1), lty=2, col="lightgray") + 
      geom_hline(aes(yintercept=0.5), lty=2, col="lightgray") + 
      #geom_hline(aes(yintercept=maj.class), lty=2, col="black") + 
      scale_colour_manual(values = as.character(colors), name = "") + 
      scale_fill_manual(values = as.character(colors)) + 
      theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(legend.position="bottom", legend.direction="horizontal") +
      guides(col = guide_legend(ncol = 3, byrow = TRUE)) +
      scale_shape_manual(values=as.numeric(pch), name = "")
  }
  
  if(!is.na(maj.class))
  {
    g <-  g +  geom_hline(aes(yintercept=maj.class), lty=2, col="black")
  }
  
  return(g)
}


#' Plot performance scores for multiple learners
#' 
#' @import RColorBrewer
#' @import ggplot2
#' @description plotComparativeResults plots a digested.results data object to compare performance results between different learners.
#' @param digested.results: a list of data.frames containing performance results from a lists of learners. This data object is returned by the function merge_digestScores()
#' @param ylim: y-axis zoom in the plot
#' @param best: a swith to plot the best values instead of declining by k_sparsity
#' @param main: name of the graphic
#' @param mode: either classification or regression (default:classification)
#' @return A list of ggplot graphs if plot is set to FALSE and a pannel organized graph otherwise.
#' @export
plotComparativeResults <- function(digested.results, plot = TRUE, ylim = c(0.5, 1), best = FALSE, ci = FALSE, main = "", mode = "classification")
{
  # estimate mode
  mode <- NULL
  if(!is.null(digested.results$list.results.digest[[1]]))
  {
    if(digested.results$list.results.digest[[1]]$best$model$objective == "cor")
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
    stop("plotComparativeResults: mode not founding stopping ...")
  }
  
  # EMPIRICAL SCORES 
  g.auc                   <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="auc_")
  g.accuracy              <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="accuracy_")
  g.recall                <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="recall_")
  g.precision             <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="precision_")
  g.f1                    <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="f1_")
  g.fit                   <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="fit_")
  g.cor                   <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="cor_")
  
  # CROSS-VALIDATION SCORES 
  # k_sparse declined
  # auc
  g.cv.empirical.auc      <- plotComparativeCV(digested.results, ylim = ylim, generalization = FALSE, score = "auc_", ci = ci, main = main)
  g.cv.generalization.auc <- plotComparativeCV(digested.results, ylim = ylim, generalization = TRUE, score = "auc_", ci = ci, main = main)
  # accuracy
  g.cv.empirical.acc      <- plotComparativeCV(digested.results, ylim = ylim, generalization = FALSE, score = "accuracy_", ci = ci, main = main)
  g.cv.generalization.acc <- plotComparativeCV(digested.results, ylim = ylim, generalization = TRUE, score = "accuracy_", ci = ci, main = main)
  # cor
  g.cv.empirical.cor      <- plotComparativeCV(digested.results, ylim = ylim, generalization = FALSE, score = "cor_", ci = ci, main = main)
  g.cv.generalization.cor <- plotComparativeCV(digested.results, ylim = ylim, generalization = TRUE, score = "cor_", ci = ci, main = main)
  
  
  # best regardless of k_sparse
  if(best)
  {
    # auc
    g.best.empirical.auc      <- plotComparativeBestCV(digested.results, ylim = ylim, generalization = FALSE, score = "auc_", ci = ci, main = main)
    g.best.generalization.auc <- plotComparativeBestCV(digested.results, ylim = ylim, generalization = TRUE, score = "auc_", ci = ci, main = main)
    
    # accuracy
    g.best.empirical.acc      <- plotComparativeBestCV(digested.results, ylim = ylim, generalization = FALSE, score = "accuracy_", ci = ci, main = main)
    g.best.generalization.acc <- plotComparativeBestCV(digested.results, ylim = ylim, generalization = TRUE, score = "accuracy_", ci = ci, main = main)
    
    # corelation
    g.best.empirical.cor      <- plotComparativeBestCV(digested.results, ylim = ylim, generalization = FALSE, score = "cor_", ci = ci, main = main)
    g.best.generalization.cor <- plotComparativeBestCV(digested.results, ylim = ylim, generalization = TRUE, score = "cor_", ci = ci, main = main)
  }
  
  
  if(plot)
  {
    # PUT PLOTS TOGETHER IN THE SAME CANVAS 
    if(mode == "classification") # if this is a classification
    {
      if(!best) # if by k_spase
      {
        multiplot(g.auc, 
                  g.cv.empirical.auc, 
                  g.cv.empirical.acc, 
                  g.accuracy,
                  g.cv.generalization.auc, 
                  g.cv.generalization.acc, 
                  cols = 2
        )
      }
      else
      {
        multiplot(g.auc, 
                  g.best.empirical.auc, 
                  g.best.empirical.acc, 
                  g.accuracy,
                  g.best.generalization.auc, 
                  g.best.generalization.acc, 
                  cols = 2
        )
      }
    }else  # if correlation
    {
      if(!best) # if by k_spase
      {
        multiplot(g.cor, 
                  g.cv.empirical.cor, 
                  g.fit,
                  g.cv.generalization.cor, 
                  cols = 2
        )
      }
      else
      {
        multiplot(g.cor, 
                  g.best.empirical.cor, 
                  g.fit,
                  g.best.generalization.cor, 
                  cols = 2
        )
      }
    }
  }else
  {
    res <- list(# empirical scores
      g.auc = g.auc, 
      g.accuracy = g.accuracy, 
      g.recall = g.recall, 
      g.precision = g.precision, 
      g.f1 = g.f1,
      g.fit = g.fit,
      g.cor = g.cor,
      # cross validation by k_sparse
      g.cv.empirical.auc = g.cv.empirical.auc, 
      g.cv.generalization.auc = g.cv.generalization.auc,
      g.cv.empirical.acc = g.cv.empirical.acc, 
      g.cv.generalization.acc = g.cv.generalization.acc,
      g.cv.empirical.cor = g.cv.empirical.cor, 
      g.cv.generalization.cor = g.cv.generalization.cor,
      # cross validation regardless of k_sparse
      g.best.empirical.auc = g.best.empirical.auc,
      g.best.generalization.auc = g.best.generalization.auc,
      g.best.empirical.acc = g.best.empirical.acc,
      g.best.generalization.acc = g.best.generalization.acc,
      g.best.empirical.cor = g.best.empirical.cor,
      g.best.generalization.cor = g.best.generalization.cor
    )
    return(res)
  }
}


#' Plot performance scores for multiple learners
#' 
#' @import RColorBrewer
#' @import ggplot2
#' @description plotComparativeResultsBest plots a digested.results data object to compare performance results between different learners focusing at the best model.
#' @param digested.results: a list of data.frames containing performance results from a lists of learners. This data object is returned by the function merge_digestScores()
#' @param ylim: y-axis zoom in the plot
#' @return A list of ggplot graphs if plot is set to FALSE and a pannel organized graph otherwise.
#' @export
plotComparativeResultsBest <- function(digested.results, plot = TRUE, ylim = c(0.5,1))
{
  
  # EMPIRICAL SCORES 
  g.auc                   <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="auc_")
  g.accuracy              <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="accuracy_")
  g.recall                <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="recall_")
  g.precision             <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="precision_")
  g.f1                    <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="f1_")
  g.fit                   <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="fit_")
  g.cor                   <- plotComparativeEmpiricalScore(digested.results, ylim = ylim, score="cor_")
  
  # CROSS-VALIDATION SCORES 
  
  # auc
  g.cv.empirical.auc      <- plotComparativeCV(digested.results, ylim = ylim, generalization = FALSE, score = "auc_")
  g.cv.generalization.auc <- plotComparativeCV(digested.results, ylim = ylim, generalization = TRUE, score = "auc_")
  
  # accuracy
  g.cv.empirical.acc      <- plotComparativeCV(digested.results, ylim = ylim, generalization = FALSE, score = "accuracy_")
  g.cv.generalization.acc <- plotComparativeCV(digested.results, ylim = ylim, generalization = TRUE, score = "accuracy_")
  
  # cor
  g.cv.empirical.cor      <- plotComparativeCV(digested.results, ylim = ylim, generalization = FALSE, score = "cor_")
  g.cv.generalization.cor <- plotComparativeCV(digested.results, ylim = ylim, generalization = TRUE, score = "cor_")
  
  if(plot)
  {
    # PUT PLOTS TOGETHER IN THE SAME CANVAS 
    if(!is.null(digested.results$empirical$auc_)) # if this is a classification
    {
      multiplot(g.auc, 
                g.cv.empirical.auc, 
                g.cv.empirical.acc, 
                g.accuracy,
                g.cv.generalization.auc, 
                g.cv.generalization.acc, 
                cols = 2
      )
    }else
    {
      multiplot(g.cor, 
                g.cv.empirical.cor, 
                g.fit,
                g.cv.generalization.cor, 
                cols = 2
      )
    }
  }else
  {
    res <- list(g.auc=g.auc, 
                g.accuracy=g.accuracy, 
                g.recall=g.recall, 
                g.precision=g.precision, 
                g.f1=g.f1,
                g.fit=g.fit,
                g.cor=g.cor,
                g.cv.empirical.auc=g.cv.empirical.auc, 
                g.cv.generalization.auc=g.cv.generalization.auc,
                g.cv.empirical.acc=g.cv.empirical.acc, 
                g.cv.generalization.acc=g.cv.generalization.acc)
    return(res)
  }
}



#load("~/Research/predomics/package/predomics/test/2016-03-09_analyse_comparative_db3/terga;db=lgc_db3_;sparsity=1_to_30;population=100;convergence=10.rda")
#load("../test/2016-03-09_analyse_comparative_db3/terga_db=lgc_db3__sparsity=1_to_30_population=100_convergence=10.rda")
# load this for the plotCors function
#source("~/Research/workspace_r/momr_1.1_corrections.R")


#### Do we need docuentation for this ?
#' @export
plotMGSQuality <- function (dat, main = "mgs", return.scores = TRUE) {
  par(mfcol = c(6, 2), xaxs = "i", yaxs = "i", mar = c(2, 2, 2, 1))
  colpan <- c("black", gplots::colorpanel(n = 20, low = "darkblue", mid = "darkorchid", high = "red"))
  size <- 50
  if (nrow(dat) < size) {
    size <- nrow(dat)
    warning("need more than 50 genes")
  }
  dat50 <- dat[1:size, ]
  plotBarcode(dat50, main = paste(main, nrow(dat50), "genes"))
  scores50 <- computeSignalMetrics(dat50)
  for (i in 1:ncol(scores50)) {
    if (all(scores50[, i] == 0)) {
      cols = "black"
    }
    else {
      cols <- colpan[as.numeric(cut(scores50[, i], breaks = 20))]
    }
    plot(scores50[, i], pch = 19, cex = 0.4, main = colnames(scores50)[i], 
         ylab = "", xlab = "individual index", col = cols)
  }
  plotBarcode(dat, main = paste(main, nrow(dat), "genes"))
  scores <- computeSignalMetrics(dat)
  for (i in 1:ncol(scores)) {
    if (all(scores[, i] == 0)) {
      cols = "black"
    }
    else {
      cols <- colpan[as.numeric(cut(scores[, i], breaks = 20))]
    }
    plot(scores[, i], pch = 19, cex = 0.4, main = colnames(scores)[i], 
         ylab = "", xlab = "individual index", col = cols)
  }
  if (return.scores) 
    return(scores)
}



# Multiple plot function
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#' @export
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) 
  {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, 
                     nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) 
  {
    print(plots[[1]])
  } else 
  {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), 
                                               ncol(layout)))
    )
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) 
    {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# # plot a horizontal barplot
# plotBarplot <- function(v, rev=TRUE, xlim=range(v), main="", col="darkgray"){
#   if(rev) v <- rev(v)
#   #barplot(v, las=2, horiz=TRUE, col="black", main=main, xlim=xlim)
#   plot.barchart <- barchart(x=v, col=col, pch=22, cex=1, xlab="", main=main)
#   return(plot.barchart)
# }


#' Plots the prevalence of a list of features in the whole dataset and per each class
#'
#' @description Plots the prevalence of a given number of features
#' @param features: a list of features or features indexes for which we wish to compute prevalence
#' @param X: dataset where to compute the prevalence
#' @param y: if provided it will also compute hte prevalence per each class (default:NULL)
#' @param topdown: showing features from top-down or the other way around (default:TRUE)
#' @param main: main title (default:none)
#' @param plot: if TRUE this provides a plot, otherwise will return different metrics such as prevalence and enrichment statistics
#' @param col.pt: colors for the point border (-1:deepskyblue4, 1:firebrick4)
#' @param col.bg: colors for the point fill (-1:deepskyblue1, 1:firebrick1)
#' @param zero.value: the value that specifies what is zero. This can be a different than 0 in log transformed data for instance (default = 0)
#' @return a ggplot object
#' @export
plotPrevalence <- function(features, X, y, topdown = TRUE, main = "", plot = TRUE, 
                           col.pt = c("deepskyblue4", "firebrick4"), 
                           col.bg = c("deepskyblue1", "firebrick1"),
                           zero.value = 0)
{
  # build object
  v.prop <- getFeaturePrevalence(features = features, X = X, y = y, prop = TRUE, zero.value = zero.value)
  v.card <- getFeaturePrevalence(features = features, X = X, y = y, prop = FALSE, zero.value = zero.value)
  v.prop.mat <- do.call(rbind, v.prop)
  v.card.mat <- do.call(rbind, v.card)
  # build melted version
  v.prop.melt <- meltScoreList(v = v.prop, prepare.for.graph = FALSE, topdown = topdown)
  colnames(v.prop.melt) <- c("feature","prevalence", "group")
  
  # convert to percentage
  v.prop.melt$prevalence <- v.prop.melt$prevalence *100
  
  # get enrichment information
  prev.enrichment <- computeCardEnrichment(v.card.mat = v.card.mat, y = y)
  if(!is.null(prev.enrichment$chisq.q))
  {
    qvals <- rep("",length(prev.enrichment$chisq.q))
    qvals[prev.enrichment$chisq.q<0.05] <- "*"
    
    # plot object
    p <- ggplot(v.prop.melt, aes(x=feature, y=prevalence, fill=group)) + 
      geom_bar(data=subset(v.prop.melt, group %in% c("all")), stat="identity", position="identity") + 
      coord_flip() + 
      geom_point(data = subset(v.prop.melt, group %in% c("-1", "1")), aes(x=feature, y=prevalence, color=group, shape=group)) + 
      scale_color_manual("Dataset", values = c("all"="gray90", "-1"=col.pt[1], "1"=col.pt[2])) +
      scale_fill_manual("Dataset", values = c("all"="gray90", "-1"=col.bg[1], "1"=col.bg[2])) +
      scale_shape_manual(values=c(25,24)) + 
      theme_bw() + 
      theme(legend.position="none", axis.text=element_text(size=9)) + 
      ggtitle(main) 
    
    if(topdown){
      p <- p +  annotate("text", y = rep(101,length(qvals)), x = seq(1,length(qvals),1)-0.3, label = rev(qvals), color="gray", size=7)
    }else{
      p <- p +  annotate("text", y = rep(101,length(qvals)), x = seq(1,length(qvals),1)-0.3, label = qvals, color="gray", size=7)
    }
    
  }else
  {
    # plot object
    p <- ggplot(v.prop.melt, aes(x=feature, y=prevalence, fill=group)) + 
      geom_bar(data=subset(v.prop.melt, group %in% c("all")), stat="identity", position="identity") + 
      coord_flip() + 
      scale_color_manual("Dataset", values = c("all"="gray90")) +
      scale_fill_manual("Dataset", values = c("all"="gray90")) +
      theme_bw() + 
      theme(legend.position="none", axis.text=element_text(size=9)) + 
      ggtitle(main)
  }
  
  if(!plot)
  {
    return(prev.enrichment)
  }else
  {
    return(p)
  }
}



#' Plots the prevalence of a list of features in the whole dataset and per each class
#'
#' @description Plots the abundance of a given number of features for each class and tests significance
#' @import reshape2
#' @param features: a list of features or features indexes for which we wish to compute prevalence
#' @param X: dataset where to compute the prevalence
#' @param y: if provided it will also compute hte prevalence per each class (default:NULL)
#' @param topdown: showing features from top-down or the other way around (default:TRUE)
#' @param main: main title (default:none)
#' @param plot: if TRUE this provides a plot, otherwise will return different metrics such as prevalence and enrichment statistics
#' @param col.pt: colors for the point border (-1:deepskyblue4, 1:firebrick4)
#' @param col.bg: colors for the point fill (-1:deepskyblue1, 1:firebrick1)
#' @return a ggplot object
#' @export
plotAbundanceByClass <- function(features, X, y, topdown = TRUE, main = "", plot = TRUE, col.pt = c("deepskyblue4", "firebrick4"), col.bg = c("deepskyblue1", "firebrick1"))
{
  check.X_y_w(X,y)
  
  if(any(is.na(match(features, rownames(X)))))
  {
    stop(paste("plotAbundanceByClass: These features are not found in the dataset",
               features[is.na(match(features, rownames(X)))]))
  }
  
  if(!is.matrix(X))
  {
    X <- as.matrix(X)
  }
  
  mode <- "classification"
  if(class(y) == "numeric" & length(table(y)) > 2)
  {
    #cat("... plotAbundanceByClass will not work for a continous y - probably in regression mode. Adapting as a uniclass\n")
    mode <- "regression"
  }
  
  if(mode == "classification")
  {
    # get levels
    lev <- names(table(y))
    
    # compute p-value of the non parametric abundance test
    if(length(features) == 1)
    {
      dat <- t(as.matrix(X[features, ]))
      rownames(dat) <- features
      datl1 <- t(as.matrix(X[features, y == lev[1]]))
      rownames(datl1) <- features
      datl2 <- t(as.matrix(X[features, y == lev[2]]))
      rownames(datl2) <- features
    }else
    {
      dat <- X[features, ]
      datl1 <- X[features, y == lev[1]]
      datl2 <- X[features, y == lev[2]]
      if(ncol(X) == 1)
      {
        dat <- as.matrix(dat)
        datl1 <- as.matrix(datl1)
        datl2 <- as.matrix(datl2)
      }
    }
    
    dat.test <- filterfeaturesK(dat, y, k = nrow(dat), sort = FALSE)
    
    if(plot)
    {
      pvals <- dat.test$p
      qvals <- rep("",nrow(dat.test))
      qvals[dat.test$q<0.05] <- "*"
      
      datl1.reshape <- melt(datl1)
      colnames(datl1.reshape) <- c("feature","observation","abundance")
      datl1.reshape$class <- rep(lev[1], nrow(datl1.reshape))
      
      datl2.reshape <- melt(datl2)
      colnames(datl2.reshape) <- c("feature","observation","abundance")
      datl2.reshape$class <- rep(lev[2], nrow(datl2.reshape))
      
      dat.reshape <- as.data.frame(t(data.frame(t(datl1.reshape), t(datl2.reshape))))
      dat.reshape$abundance <- as.numeric(as.character(dat.reshape$abundance))
      
      # fix factor level order
      if(topdown) 
      {
        # use the same factor levels as features
        dat.reshape$feature <- factor(dat.reshape$feature, levels=rev(features))
      }else
      {
        # use the same factor levels as features
        dat.reshape$feature <- factor(dat.reshape$feature, levels=features)
      }
      
      # plot object
      p <- ggplot(dat.reshape, aes(x=feature, y = abundance, fill=class, color=class)) +
        geom_boxplot() + 
        #scale_x_continuous(limits = range(dat.reshape$abundance)) +
        coord_flip() +
        #facet_grid(. ~ class) +
        theme_bw() +
        scale_color_manual(values = col.pt) +
        scale_fill_manual(values = col.bg) +
        theme(legend.position="none") +
        ggtitle(main)
      
      pad <- max(dat.reshape$abundance) + max(dat.reshape$abundance)*0.1
      
      if(topdown){
        p <- p +  annotate("text", y = rep(pad,length(qvals)), x = seq(1,length(qvals),1) - 0.3, label = rev(qvals), color="gray", size=7)
      }else{
        p <- p +  annotate("text", y = rep(pad,length(qvals)), x = seq(1,length(qvals),1) - 0.3, label = qvals, color="gray", size=7)
      }
      return(p)
    }else
    {
      return(dat.test)
    }
  }else # mode regression
  {
    # get levels
    
    # compute p-value of the non parametric abundance test
    if(length(features) == 1)
    {
      dat <- t(as.matrix(X[features, ]))
      rownames(dat) <- features
    }else
    {
      dat <- X[features, ]
      if(ncol(X) == 1)
      {
        dat <- as.matrix(dat)
      }
    }
    
    # we can still correlate and compute p-values
    dat.test <- filterfeaturesK(dat, y, k = nrow(dat), sort = FALSE)
    
    if(plot)
    {
      pvals <- dat.test$p
      qvals <- rep("",nrow(dat.test))
      qvals[dat.test$q<0.05] <- "*"
      
      dat.reshape <- melt(dat)
      colnames(dat.reshape) <- c("feature","observation","abundance")
      dat.reshape$class <- rep("all", nrow(dat.reshape))
      
      # fix factor level order
      if(topdown) 
      {
        # use the same factor levels as features
        dat.reshape$feature <- factor(dat.reshape$feature, levels=rev(features))
      }else
      {
        # use the same factor levels as features
        dat.reshape$feature <- factor(dat.reshape$feature, levels=features)
      }
      
      # plot object
      p <- ggplot(dat.reshape, aes(x=feature, y = abundance, fill=class, color=class)) +
        geom_boxplot() + 
        #scale_x_continuous(limits = range(dat.reshape$abundance)) +
        coord_flip() +
        #facet_grid(. ~ class) +
        theme_bw() +
        scale_color_manual(values = "gray40") +
        scale_fill_manual(values = "gray80") +
        theme(legend.position="none") +
        ggtitle(main)
      
      pad <- max(dat.reshape$abundance) + max(dat.reshape$abundance)*0.1
      
      if(topdown){
        p <- p +  annotate("text", y = rep(pad,length(qvals)), x = seq(1,length(qvals),1) - 0.3, label = rev(qvals), color="gray", size=7)
      }else{
        p <- p +  annotate("text", y = rep(pad,length(qvals)), x = seq(1,length(qvals),1) - 0.3, label = qvals, color="gray", size=7)
      }
      return(p)
    }else
    {
      return(dat.test)
    }
  }
}



#' Plots the prevalence of a list of features in the whole dataset and per each class
#'
#' @description Plots the coefficients of subset of features in the models where they are found
#' @importFrom reshape2 melt
#' @param feat.model.coeffs: feature vs. model coeffient table
#' @param topdown: showing features from top-down or the other way around (default:TRUE)
#' @param main: main title (default:none)
#' @param col: colors to be used for the coeffients (default: -1 = deepskyblue1, 0 = white, 1 = firebrick1)
#' @param vertical.label: wether the x-axis labels should be vertical or not (default:TRUE)
#' @return a ggplot object
#' @export
plotFeatureModelCoeffs <- function(feat.model.coeffs, topdown = TRUE, main="", col = c("deepskyblue1","white","firebrick1"), vertical.label = TRUE)
{
  
  # get the corresponding data
  data <- feat.model.coeffs 
  
  # prepare the data for the 
  if(topdown) 
  {
    data <- t(data[nrow(data):1, ])
  }else
  {
    data <- t(data)
  }
  #rownames(data) <- c(1:nrow(data))
  prev <- signif(colSums(data!=0)/nrow(data),1)
  # melt the data
  data.m <- melt(data)
  colnames(data.m) <- c("models","feature","value")
  
  col.n <- c("-1","0","1")
  tab.v <- table(data.m$value)
  if(length(tab.v)<3)
  {
    col <- col[col.n %in% names(tab.v)]
  }
  
  p <- ggplot(data.m, aes(models, feature)) + 
    geom_tile(aes(fill = value), colour = "darkgray") + 
    theme_bw() + 
    scale_fill_gradientn(colours = col) +
    ggtitle(main) 
  
  if(vertical.label)
  {
    p <- p + theme(legend.position="none", axis.text=element_text(size=9), axis.text.x = element_text(angle = 90, hjust = 1))
  }else
  {
    p <- p + theme(legend.position="none", axis.text=element_text(size=9))
  }
  
  return(p)
}


# # usees gridExtra
# plotUsedFeatures <- function(feature.data=NULL, y = y, feature.presence=NULL, topdown = TRUE){
#   p1 <- plotCoeffsPop(data = feature.presence, main = "Model appearence", topdown = topdown)
#   p2 <- plotPrevalence(features = feature.data, y = y, main="X/y prevalence", topdown = topdown)
#   return(grid.arrange(p1, p2, ncol=2, widths=c(3,1)))
# }


################################################################
# POPULATION PLOTS
################################################################
#' @title Plots a population of models (or a single model) objects as barplots of scaled coefficients.
#'
#' @description Plots an model or a population of models as a barplots, representing each feature, the length being the coefficient
#' @import gridExtra
#' @param pop: a population of models to plot
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the class vector
#' @param sort.features: wether the features need to be sorted by correlation with 'y' or not
#' @param sort.ind: computing sorting can take time if computed for every model and can be computed outside the function and passed as a parameter
#' @param col.sign: the colors of the cofficients based on the sign of the coefficients (default: -1=deepskyblue1,1:firebrick1)
#' @param ncol: number of graphics for each line (default: 10)
#' @param slim: plot without axis information (default:FALSE)
#' @param importance: the importance (mda) of the features in crossval
#' @export
plotPopulation <- function(pop, X, y, 
                           sort.features = FALSE, sort.ind = NULL, 
                           col.sign = c("deepskyblue1", "firebrick1"), ncol = 10, 
                           slim = FALSE, importance = FALSE)
{
  
  # test population validity
  if(!isPopulation(obj = pop) & !isModel(obj = pop))
  {
    stop("plotPopulation: a population of models shouls be provided.")
  }
  
  if(length(col.sign) != 2) 
  {
    print("plotPopulation: please provide 2 colors for the ternary coefficients excluding zero")
    return(NULL)
  }
  
  # if this is a model plot the model
  if(isModel(obj = pop))
  {
    g <- plotModel(pop, X, y, sort.features=sort.features, sort.ind = sort.ind, col.sign = col.sign, slim = slim, importance = importance)
    return(g)
  }
  
  # else prepare the data to print a population
  if(sort.features)
  {
    if(is.null(sort.ind))
    {
      # get the order of the features in terms of discriminance compared to the class.
      ind       <- order(filterfeaturesK(X, y, k=nrow(X), sort = FALSE)$p, decreasing = FALSE)
    }else
    {
      ind       <- sort.ind
    }
  }else # no order (use the default X ordering)
  {
    ind         <- c(1:nrow(X))
  }
  
  # put all the graphs here
  list.graphs <- list()
  for(i in 1:length(pop))
  {
    list.graphs[[i]] <- plotModel(mod = pop[[i]], X, y, 
                                  sort.features = sort.features, sort.ind = ind, 
                                  col.sign = col.sign, 
                                  slim = slim, importance = importance)
  }
  
  do.call("grid.arrange", c(list.graphs, ncol=ncol))
}



################################################################
# SINGLE MODEL PLOTS
################################################################

#' @title Plots a model or a population of model objectsas barplots of scaled coefficients.
#'
#' @description Plots a model or a population of models as a barplots, representing each feature, the length being the coefficient
#' @param mod: a model to plot
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the class vector
#' @param sort.features: wether the features need to be sorted by correlation with 'y' or not (default: TRUE)
#' @param sort.ind: computing sorting can take time if computed for every model and can be computed outside the function and passed as a parameter
#' @param feature.name: show the name of the features (default:FALSE)
#' @param col.sign: the colors of the cofficients based on the sign of the coefficients (default: -1=deepskyblue1, 1:firebrick1)
#' @param main: possibility to change the title of the function (default:"")
#' @param slim: plot without axis information (default:FALSE)
#' @param importance: the importance (mda) of the features in crossval
#' @param res_clf: the result of the learning process (default:NULL). If provided information on MDA will be extracted for the importance graphic.
#' @export
plotModel <- function(mod, X, y, 
                      sort.features = FALSE, 
                      sort.ind = NULL, 
                      feature.name = FALSE,
                      col.sign = c("deepskyblue1", "firebrick1"), 
                      main = "", 
                      slim = FALSE, 
                      importance = FALSE,
                      res_clf = NULL)
{
  
  # test model validity
  if(!isModel(obj = mod))
  {
    print("plotModel: The model object is not valid!")
    return(NULL)
  }
  
  if(length(col.sign) != 2) 
  {
    print("plotModel: please provide 2 colors for the ternary coefficients excluding zero")
    return(NULL)
  }
  
  # disable importance for SOTA
  if(isModelSota(mod) & importance)
  {
    importance <- FALSE
    print("plotModel: importance graphic is disabled for SOTA models")
  }
  
  if(!isModelSotaRF(mod))
  {
    
    # Fix display names
    if(isModelSotaSVM(mod))
    {
      mod$learner <- "sota"
    }
    
    # reset model attributes for glmnet for proper viewing
    if(mod$learner == "terda" & mod$language == "logreg")
    {
      mod$learner <- "sota"
      mod$language <- "glmnet"
    }
    
    if(sort.features)
    {
      if(is.null(sort.ind))
      {
        # get the order of the features in terms of discriminance compared to the class.
        ind       <- order(filterfeaturesK(data = X, trait = y, k=nrow(X), sort = FALSE)$p, decreasing = FALSE)
      }else
      {
        ind       <- sort.ind
      }
    }else # no order (use the default X ordering)
    {
      ind         <- c(1:nrow(X))
    }
    
    # get the normalized coefficients
    coeffsl <- normModelCoeffs(mod = mod, X = X, y = y, sort.features = sort.features, sort.ind = ind)
    # and the sparsity
    k_sparsity <- mod$eval.sparsity
    #k_sparsity <- k_sparsity[!unlist(lapply(coeffsl,is.null))]
    #tmp <- unlist(lapply(coeffsl,length))
    coeffs.name <- c()
    coeffs.nr   <- c()
    
    coeffs.name  <- rep(mod$learner,length(coeffsl))
    coeffs.nr    <- c(1:length(coeffsl))  
    
    coeffs.data         <- data.frame(coeffsl, coeffs.name, coeffs.nr)
    colnames(coeffs.data) <- c("coeff","classifier","feature")
    coeffs.data$coeff   <- as.numeric(coeffs.data$coeff)
    coeffs.data$sign    <- sign(coeffs.data$coeff)
    coeffs.data$sign[coeffs.data$sign == 0] <- NA
    coeffs.data$col     <- col.sign[factor(coeffs.data$sign, levels = c(-1,1))]
    
    # add information on importance
    rownames(coeffs.data) <- rownames(X)[ind]
    coeffs.data$importance <- 0
    coeffs.data$importance.col <- NA
    
    if(!is.null(mod$mda.cv_))
    {
      coeffs.data[mod$names_,]$importance <- mod$mda.cv_
      coeffs.data[mod$names_,]$importance.col <- "black"
    }
    
    # get the features
    features <- rownames(coeffs.data)[which(!is.na(coeffs.data$sign))]
    
    names(col.sign) <- c("-1","1")
    col.sign <- as.character(col.sign[names(table(sign(coeffs.data$coeff[coeffs.data$coeff != 0])))])
    
    # make the main plot
    g1 <- ggplot(coeffs.data, aes(feature, coeff, fill=col)) + 
      geom_bar(stat="identity", position="dodge") + ylim(-1, 1) + 
      xlab("") +
      theme(legend.position="none", axis.text=element_text(size=9)) + 
      scale_fill_manual("Sign", values = col.sign) +
      geom_hline(yintercept = 0, col="gray") +
      theme_bw() + guides(fill = "none") +
      coord_flip() + 
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(colour = "white", fill = "white", size = 0.3),
        strip.text.x = element_text(colour = "darkred", size = 10))
    
    if(feature.name)
    {
      g1 <- g1 + scale_x_continuous(breaks=which(!is.na(coeffs.data$sign)), labels = features) + theme(aspect.ratio = 2) 
    }else
    {
      g1 <- g1 + scale_x_continuous(breaks=which(!is.na(coeffs.data$sign)))
    }
    
    if(!slim)
    {
      if(main == "")
      {
        main <- paste("alg:", mod$learner, " | lang:",mod$language," | k:",mod$eval.sparsity, sep="")
      }
      g1 <- g1 + 
        ggtitle(main)
    }else
    {
      if(main == "")
      {
        main <- paste(mod$learner, "|",mod$language,"|",mod$eval.sparsity, sep="")
      }
      
      g1 <- g1 + 
        ggtitle(main) +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank())
    }
    
    # make the plot on importance
    if(!is.null(mod$mda.cv_) & importance)
    {
      importance.se <- FALSE
      if(!is.null(res_clf))
      {
        if(!is.null(res_clf$crossVal$fip))
        {
          importance.se <- TRUE
        }
      }
      
      if(importance.se)
      {
        
        mdacv <- mergeMeltImportanceCV(list.results = list(exp = res_clf), 
                                       filter.cv.prev = 0, 
                                       min.kfold.nb = FALSE, 
                                       learner.grep.pattern = "*", 
                                       nb.top.features = NULL,
                                       feature.selection = features,
                                       scaled.importance = FALSE,
                                       make.plot = TRUE,
                                       main = TRUE,
                                       cv.prevalence = FALSE)
        g2 <- mdacv$g + theme(aspect.ratio = 2) + ggtitle("mda|cv")
        
      }else
      {
        g2 <- ggplot(coeffs.data, aes(feature, importance, fill=importance.col)) + 
          geom_bar(stat="identity", position="dodge") + 
          theme(legend.position="none", axis.text=element_text(size=9)) + 
          scale_fill_manual("Sign", values = "black") +
          xlab("") +
          geom_hline(yintercept = 0, col="gray") +
          theme_bw() + guides(fill = "none") +
          coord_flip() + 
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.background = element_rect(colour = "white", fill = "white", size = 0.3),
            strip.text.x = element_text(colour = "darkred", size = 10)
          )
        
        if(feature.name )
        {
          g2 <- g2 + scale_x_continuous(breaks=which(!is.na(coeffs.data$sign)), labels = features) + theme(aspect.ratio = 2) +
            ggtitle("mda|cv")
          
        }else
        {
          g2 <- g2 + scale_x_continuous(breaks=which(!is.na(coeffs.data$sign)))
          
          if(slim)
          {
            g2 <- g2 + 
              ggtitle("mda|cv") + 
              ylim(0, 1) + 
              theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                strip.background = element_rect(colour = "white", fill = "white", size = 0.3),
                strip.text.x = element_text(colour = "darkred", size = 10),
                axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.title.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.text.y=element_blank()
              )
          }else
          {
            g2 <- g2 + 
              ggtitle("mda|cv")
          }
          
        } # end zoom in features
      } # end CV mda plot type
      
      #return(grid.arrange(g, g2, ncol=2, widths = c(2,1)))
      return(g2)
    }else # end else plot importance
    {
      warning("plotModel: importance is not available in this model, returning model plot")
      return(g1)
    }
  }else
  {
    # RANDOM FOREST MODEL
    mod$learner <- "sota"
    
    if(main == "")
    {
      main <- paste(mod$learner, "|",mod$language,"|",mod$eval.sparsity, sep="")
    }
    
    if(!is.null(mod$obj))
    {
      g <- tree_func(final_model = mod$obj, tree_num = 1, main = main, node.text.color = "black")  
    }else
    {
      # return an empty graph
      df <- data.frame(
        x = rep(1, 1),
        y = rep(1, 1),
        z = rep(1, 1)
      )
      
      g <- ggplot(df, aes(x, y, fill = "black")) +
        geom_tile() +
        ggtitle(main) +
        #theme(legend.position="none", axis.text=element_text(size=9)) + 
        scale_fill_manual("Sign", values = "lightgray") +
        theme_bw() + guides(fill = "none") +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              legend.position="none",
              panel.background=element_blank(),
              #panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()#,
              #plot.background=element_blank()
        ) +
        coord_flip()
    }
  }
  
  return(g)
}


#' @title Plots a model or a population of model objectsas barplots of scaled coefficients.
#' @import ggplot2
#' @description Plots a model score or a population of models as a barplots, representing each feature, the length being the coefficient
#' @param mod: a model to plot
#' @param y: the class to predict
#' @param col.sign: the colors of the cofficients based on the sign of the coefficients (default: -1=deepskyblue1, 1:firebrick1)
#' @param main: possibility to change the title of the function (default:"")
#' @export
plotModelScore <- function(mod = NULL,
                           y = NULL, # the class to predict
                           col.sign = c("deepskyblue1", "firebrick1"), 
                           main = ""
)
{
  
  # test model validity
  if(!isModel(obj = mod))
  {
    print("plotModelScore: The model object is not valid!")
    return(NULL)
  }
  
  if(length(col.sign) != 2) 
  {
    print("plotModelScore: please provide 2 colors for the ternary coefficients excluding zero")
    return(NULL)
  }
  
  # disable this plot for SOTA
  if(isModelSota(mod))
  {
    print("plotScore: This method is not available for SOTA models")
    return(NULL)
  }
  
  # check the existence of the models score
  if(!myAssertNotNullNorNa(mod$score_))
  {
    print("plotModelScore: The score_ attribute is missing. Returning empty handed, please provide it to the model.")
    return(NULL)
  }
  
  # check the existence of the models score
  if(!myAssertNotNullNorNa(y))
  {
    print("plotModelScore: The parameter y is missing. Returning empty handed, please provide it.")
    return(NULL)
  }
  
  # check the existence of the models score
  if(length(mod$score_) != length(y))
  {
    print("plotModelScore: The score_ attribute does not match y in length. Returning empty handed, please verify the input data.")
    return(NULL)
  }
  
  # reset model attributes for glmnet for proper viewing
  if(mod$learner == "terda" & mod$language == "logreg")
  {
    mod$learner <- "sota"
    mod$language <- "glmnet"
  }
  
  if(mod$objective == "auc")
  {
    y <- as.factor(y)
    df <- data.frame(yhat = mod$score_,
                     y = y)
    
    main <- paste("Acc =", signif(mod$accuracy_,2), 
                  " | AUC =", signif(mod$auc_,2), 
                  " | Size =", signif(mod$eval.sparsity,2), 
                  " | Lang =", mod$language)
    g1 <- ggplot(data = na.omit(df), aes(x = y, y = yhat, color = y)) + 
      geom_boxplot(aes(alpha = 0.1, color = y, fill = y)) +
      geom_jitter(width = 0.2) +
      geom_hline(yintercept = mod$intercept_, linetype = "dashed", col = "darkgray") +
      scale_color_manual(values = col.sign) +
      scale_fill_manual(values = col.sign) +
      ggtitle(main) +
      ylab("y^") +
      theme_bw() +
      theme(aspect.ratio = 1) + 
      theme(legend.position="none")
    
  }else
  {
    df <- data.frame(yhat = mod$score_,
                     y = y)
    main <- paste("Rho =", signif(mod$cor_,2), 
                  " | R2 =", signif(mod$rsq_,2), 
                  " | SER =", signif(mod$ser_,2))
    g1 <- ggplot(data = na.omit(df), aes(x = y, y = yhat, alpha = 0.9)) + 
      geom_smooth(method = 'lm') + 
      geom_point(aes(color = "darkred", size = 4)) +
      ggtitle(main) +
      ylab("y^") +
      theme_bw() +
      theme(aspect.ratio = 1) + 
      theme(legend.position="none")
  }
  
 return(g1)
  
}


#' @title Normalize the model coefficients needed for the plot
#'
#' @description Normalize the model coefficients needed for the plot
#' @param mod: a model to plot
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the class vector
#' @param sort.features: wether the features need to be sorted by correlation with 'y' or not (default:FALSE)
#' @param sort.ind: computing sorting can take time if computed for every model and can be computed outside the function and passed as a parameter
#' @return the normalized coefficients
#' @export
normModelCoeffs <- function(mod, X, y, sort.features = FALSE, sort.ind = NULL)
{
  
  # test model validity
  if(!isModel(obj = mod))
  {
    stop("normModelCoeffs: Please make sure the model object is valid.")
  }
  
  if(is.null(mod$coeffs_))
  {
    warning("normModelCoeffs: This model does not have any coefficients and can not be plotted")
    return(NULL)
  }
  
  if(any(is.na(mod$coeffs_)))
  {
    warning("normModelCoeffs: This model does not have any coefficients and can not be plotted")
    return(NULL)
  }
  
  # transform the model onto the dense format
  bm.sparse     <- modelToDenseVec(natts = nrow(X), mod = mod)
  
  if(sort.features)
  {
    if(is.null(sort.ind))
    {
      # get the order of the features in terms of discriminance compared to the class.
      ind       <- order(filterfeaturesK(X, y, k=nrow(X), sort = FALSE)$p, decreasing = FALSE)
    }else
    {
      ind       <- sort.ind
    }
  }else # no order (use the default X ordering)
  {
    ind         <- c(1:nrow(X))
  }
  
  
  # order features
  bm.sparse.ind <- bm.sparse[ind]
  # normalize
  bm.sparse.ind <- bm.sparse.ind/max(c(abs(bm.sparse.ind),1)) # for glmnet or SVM or other
  
  res <- NULL
  # use the different languages plot
  switch(mod$language,
         bin=,
         bininter={
           # case 'bininter' here...
           res <- abs(bm.sparse.ind)
         },
         ter = ,
         terinter = ,
         svm = ,
         logreg = ,
         glmnet = {
           # case 'terinter' here...
           res <- bm.sparse.ind
         },
         ratio = {
           # case 'ratio' here...
           if(is.na(mod$intercept_)) # for the case when it is NA as for the regression
           {
             bm.sparse.ind.ratio <- bm.sparse.ind
             bm.sparse.ind.ratio[bm.sparse.ind.ratio>0] <- bm.sparse.ind.ratio[bm.sparse.ind.ratio>0] / 1
           }else
           {
             bm.sparse.ind.ratio <- bm.sparse.ind
             bm.sparse.ind.ratio[bm.sparse.ind.ratio>0] <- bm.sparse.ind.ratio[bm.sparse.ind.ratio>0]/mod$intercept_  
           }
           bm.sparse.ind.ratio <- bm.sparse.ind.ratio/max(c(abs(bm.sparse.ind.ratio),1)) # normalize
           res <- bm.sparse.ind.ratio
         },
         {
           warning("normModelCoeffs: this language is not implemented")
         }
  )
  return(res)
}



# TODO convert in ggplot

#' Analyzes the score construction and model
#'
#' @description Analyzes the score construction and model
#' @param mod: a model object where the score will be computed
#' @param X: the data matrix with variables in the rows and observations in the columns
#' @param y: the class vector
#' @param clf: an object containing the different parameters of the classifier
#' @param plot: plot graphical interpretation of if TRUE, (default:TRUE)
#' @return an object containing statistics on a given model
#' @export
disectModel <- function(mod, X, y, clf, plot = TRUE)
{
  if(!isModel(obj = mod))
  {
    stop("disectModel: please provide a valid model!")    
  }
  
  if(isModelSota(obj = mod))
  {
    stop("disectModel: please provide a valid ternary model! This method is not adapted to SOTA models.")    
  }
  
  # discet result
  dres <- list()
  dres$mod <- mod
  dres$y <- y
  
  # Compute the positive and negative score
  pos <- which(mod$coeffs_ > 0)
  neg <- which(mod$coeffs_ < 0)
  pos <- list(indices_ = mod$indices_[pos], coeffs_ = mod$coeffs_[pos])
  neg <- list(indices_ = mod$indices_[neg], coeffs_ = mod$coeffs_[neg])
  
  if(length(pos$indices_) == 0) # send an zero filled vector (another idea to send out a NULL object)
  {
    pos_score <- rep(0, ncol(X)); names(pos_score) <- colnames(X)
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
  
  if(length(neg$indices_) == 0) # send an zero filled vector (another idea to send out a NULL object)
  {
    neg_score <- rep(0, ncol(X)); names(neg_score) <- colnames(X)
    neg_score <- t(as.matrix(neg_score))
  } else if(length(neg$indices_) == 1) # if the positive part is of length one
  {
    neg_data <- as.matrix(X[neg$indices_,])
    neg_score <- neg_data
  } else # If greater than one
  {
    neg_data <- X[neg$indices_,]
    neg_score <- (-neg$coeffs_) %*% as.matrix(neg_data) # compute score ^y
  }
  
  dres$scores <- t(data.frame(t(pos_score), t(neg_score)*-1, mod$score_))
  rownames(dres$scores) <- c("pos","neg","score")  
  
  # normalization function from momr. No need to importa all the package for this
  normFreqTC <- function (dat) 
  {
    if (!is.matrix(dat)) 
    {
      dat <- as.matrix(dat)
    }
    res <- dat
    for (i in 1:ncol(res)) {
      res[, i] <- res[, i]/sum(res[, i])
    }
    return(res)
  }
  
  # normalize by score
  dres$scores.norm <- t(normFreqTC(dres$scores))
  
  
  if(plot)
  {
    pmar <- par()$mar
    par(mfrow=c(2,1))
    par(mar = c(5,5,4.1,0))
    plotScoreBarcode(dscore = dres$scores, y = y, nb.col.levels = 30, main="score composition")
    
    # classification of samples as follwing the score
    par(mar = c(5,5,0,0))
    ord.score <- order(dres$mod$score_)
    plot(dres$mod$score_[ord.score], pch='', xlab="observatons", ylab="score_")
    points((1:length(y))[y[ord.score]=="-1"],dres$mod$score_[ord.score][y[ord.score]=="-1"], pch='+', col="orange")
    points((1:length(y))[y[ord.score]=="1"],dres$mod$score_[ord.score][y[ord.score]=="1"], pch='+', col="darkgreen")
    abline(h=mod$intercept_, col="darkgray", lty=2)
    abline(v=max(which(dres$mod$score_[ord.score]<mod$intercept_)), col="darkgray", lty=2)
    legend("topleft",legend = c(paste("accuracy:",signif(mod$accuracy_,3)), 
                                paste("intercept:", signif(mod$intercept_,3))), 
           border = "white", bty="n")
    par(mar = pmar)
  }
  
  return(dres)
}



#' Plots the barcode of the total score as well as positive and negative components
#'
#' @description Plots the barcode of the total score as well as positive and negative components
#' @importFrom viridis viridis_pal
#' @param dscore: an object containing different statistics on a model
#' @param y: the class vector
#' @param clf: an object containing the different parameters of the classifier
#' @param nb.col.levels: number of distinct colors from the viridis palette (default:30)
#' @param main: a title for the graphic
#' @return nothing
plotScoreBarcode <- function(dscore, y, nb.col.levels = 30, main="")
{
  ord.y <- order(y)
  cat.y <- table(y)
  
  # reorder by class
  dscore <- dscore[,ord.y]
  
  # build the color space and breaks
  cols <- list()
  cols$val <- c(seq(min(dscore),max(dscore), length.out = nb.col.levels+1))
  cols$col <- viridis_pal()(nb.col.levels)
  
  # build the image
  x.coord <- (1:nrow(dscore))
  y.coord <- (1:ncol(dscore))
  image(y.coord, x.coord, t(dscore[rev(x.coord),]), breaks = cols$val, col = cols$col, axes=FALSE, xlab="", ylab="", main=main)
  
  # add axis  
  axis(2, at = 1:nrow(dscore), labels=rownames(dscore)[rev(x.coord)],  las = 1, tick=FALSE)
  class.x <- c(cat.y[1]/2, cat.y[1]+cat.y[2]/2)
  axis(1, at = class.x, labels=c("-1","1"), tick=FALSE)
  
  # add grid
  grid(ny = 3, nx=1, col="white", lwd=1, lty=1)
  abline(v=cat.y[1], lwd=1, lty=1, col="white")
  box(col = "black", cex = 0.5, which = "plot")
  
}



#' Analyze the results from a given classifier
#'
#' @description Analyze the results from a given classifier.
#' @param score: this is the y^ of a given model
#' @param y: the class to be predted
#' @param main: title of the graph
#' @param ci: the point shape for the graph
#' @param percent: color for the graph
#' @return a roc object
#' @export
plotAUC <- function(score, y, main="", ci = TRUE, percent = TRUE)
{
  require(pROC)
  rocobj <-  roc(response = y, predictor = score, percent = percent, ci = ci, of = "se", sp = seq(0, 100, 5))
  plot(rocobj, ci.type="shape", ci.col="grey80", main=main)
  # compute information on the threshold
  rocobj2 <- roc(response = y, predictor = score, percent = percent, ci = TRUE, of = "auc")
  resa  = coords(rocobj2, x = "best", input = "threshold", best.method = "youden")
  abline(v=resa[2], col="red", lty=2); abline(h=resa[3], col="red", lty=2)
  legend("bottomright",legend = c(paste("auc:",signif(rocobj$auc,3)),
                                  paste("ci:",signif(rocobj2$ci,3)[1],"-",signif(rocobj2$ci,3)[3]),
                                  paste("threshold:",signif(resa[1],3))))
  return(rocobj)
}


# plotAUCg <- function(mod = NULL, verbose = FALSE, show.intercept = TRUE)
# {
#   if(!isModel(mod))
#   {
#     if(verbose) warning("plotAUCg: please provide a valid model.")
#     return(NULL)
#   }
#   
#   # make the roc object
#   robj <- roc(y ~ mod$score_, plot=FALSE, ci = TRUE)
#   # get the coordonates of the best
#   resa  = coords(robj, x = "best", input = "threshold", best.method = "youden")
#   # build the legend
#   legend = c(paste("AUC:",signif(robj$auc,3)),
#              paste("CI:",signif(robj$ci,3)[1],"-",signif(robj$ci,3)[3]),
#              paste("Threshold:",signif(resa[1],3)))
#   # make the plot
#   g <- ggroc(robj) +
#     annotate(geom = "text", x = 0.25, y = 0.1,  label = paste(legend, collapse = "\n"), hjust=0) +
#     theme_bw()
#   
#   # add intercept
#   if(show.intercept)
#   {
#     g <- g + annotate(geom = "text", x = resa["specificity"], y = resa["sensitivity"],  label = "+", colour = "red", size = 10)
#   }
#   
#   return(g)
# }
# 
# 

#' Plot the AUC of a given classifier
#'
#' @description Analyze the results from a given classifier.
#' @import pROC
#' @param mod: a predomics model object (default = NULL)
#' @param score: this is the y^ of a given model
#' @param y: the class to be predicted
#' @param main: title of the graph
#' @param ci: the point shape for the graph
#' @param show.intercept: plot or not the intercept on the graph (default:TRUE)
#' @return a ggplot object
#' @export
plotAUCg <- function(mod = NULL, score, y, main = "", ci = TRUE, show.intercept = TRUE)
{
  
  if(isModel(mod))
  {
    print("plotAUCg: using the model object to get the score")
    score <- mod$score_
  }
  
  if(length(y) != length(score))
  {
    warning("plotAUCg: the score should be the same length as y. Returning empty-handed !")
    return(NULL)
  }
  
  # for the ratio models we may have infinite values which will make the CI estimation to fail
  # here is a workaround
  
  ind <- is.infinite(score) | is.na(score) | is.nan(score)
  
  # compute the roc object
  rocobj <- roc(y[!ind] ~ score[!ind], plot=FALSE, ci = ci)
  tpr <- mod$confusionMatrix_[1,1]/ sum(mod$confusionMatrix_[,1])
  fpr <- mod$confusionMatrix_[1,2]/ sum(mod$confusionMatrix_[,2])
  res.ci <- ci(rocobj)
  legend <- paste0("AUC = ",signif(res.ci[2],2), 
                   "\nCI = [", signif(res.ci[1],2),";", signif(res.ci[3],2),"]")
  
  g <- ggroc(rocobj, colour = "black", linetype = 1, size = 1, legacy.axes = TRUE) +   
    theme_bw() + 
    annotate("text",
             x = 0.75,
             y = 0.1,
             label = legend, hjust = 0) +
    guides(fill = "none") +
    ggtitle(main) +
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      strip.background = element_rect(colour = "white", fill = "white", size = 0.3),
      strip.text.x = element_text(colour = "darkred", size = 10))
  
  resa  = coords(rocobj, x = "best", input = "threshold", best.method = "youden")
  
  # add intercept
  if(show.intercept)
  {
    if(sum(ind)==0) # everything fine
    {
      g <- g + annotate("point", 
                        x = fpr, 
                        y = tpr, 
                        colour = "red", 
                        size = 10, 
                        pch = "+")
    }else # potential no numbers
    {
      g <- g + annotate("point", 
                        x = 1-resa["specificity"], 
                        y = resa["sensitivity"], 
                        colour = "red", 
                        size = 10, 
                        pch = "+")
    }
   
  }
  return(g)
}



# # With a roc object:
# library(pROC)
# rocobj <- roc(response = y, predictor = score)
# rocobj.ci <- ci(rocobj, of="thresholds", thresholds="all", boot.n=100)
# plot(x = 1-rocobj$specificities, y = rocobj$sensitivities, type="l")
# points(x = 1-rocobj.ci$specificity[,"2.5%"], y = rocobj.ci$sensitivity[,"2.5%"], type="l", col="red")
# points(x = 1-rocobj.ci$specificity[,"97.5%"], y = rocobj.ci$sensitivity[,"97.5%"], type="l", col="red")
# points(x = 1-rocobj.ci$specificity[,"50%"], y = rocobj.ci$sensitivity[,"50%"], type="l", col="red")





# ggroc <- function(roc, showAUC = TRUE, interval = 0.2, breaks = seq(0, 1, interval)){
#   require(pROC)
#   if(class(rocobj) != "roc")
#     simpleError("Please provide roc object from pROC package")
#   plotx <- rev(rocobj$specificities)
#   ploty <- rev(rocobj$sensitivities)
#   
#   ggplot(NULL, aes(x = plotx, y = ploty)) +
#     geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5) + 
#     geom_step() +
#     scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks, expand = c(0.001,0.001)) + 
#     scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001)) +
#     theme_bw() + 
#     theme(axis.ticks = element_line(color = "grey80")) +
#     coord_equal() + 
#     annotate("text", x = interval/2, y = interval/2, vjust = 0, label = paste("AUC =",sprintf("%.3f",rocobj$auc)))
# }



#' # plot a horizontal barplot
#' #' @export
#' plotBarplot <- function(v, rev=TRUE, xlim=range(v), main=""){
#'   if(rev) v <- rev(v)
#'   barplot(v, las=2, horiz=TRUE, col="black", main=main, xlim=xlim)
#' }




################################################################
# PRINTING DIFFERENT, OBJECTS
################################################################

#' Prints a model object as text.
#'
#' @description Prints a model object as text
#' @param mod: a model to plot
#' @param method: an object containing the different parameters of the classifier
#' @param score: which score to show in the fit (default:fit_)
#' @export
printModel <- function(mod, method = "short", score = "fit_")
{
  if(!isModel(obj = mod))
  {
    print("printModel: please provide a valid predomics model object.")
    return(NULL)
  }
  
  if(!score %in% names(mod))
  {
    print("printModel: please provide a valid score that is found as an attribute in the model object.")
    return(NULL)
  }
  
  # if(isModelSota(mod))
  # {
  #   method <- "sota"
  # }
  
  switch(method, 
         short={
           if(!isModelSota(mod))
           {
             if(all(myAssertNotNullNorNa(mod$coeffs_)))
             {
               # Bin(inter) or Ter(inter)
               ind.pos <- sign(mod$coeffs_) == 1
               ind.neg <- sign(mod$coeffs_) == -1
               
               # BTR
               if(mod$language == "ratio")
               {
                 signs <- rep("+",length(mod$coeffs_))
                 if(all(!ind.pos))
                 {
                   term.pos <- "0"
                 }else
                 {
                   term.pos <- paste("(",paste(signs[ind.pos], mod$indices_[ind.pos], sep = "", collapse = ""),")",sep="")  
                 }
                 
                 if(all(!ind.neg))
                 {
                   term.neg <- "0"
                 }else
                 {
                   term.neg <- paste("(",paste(signs[ind.neg], mod$indices_[ind.neg], sep = "", collapse = ""),")",sep="")
                 }
                 res <- paste(term.pos, "/", term.neg)
                 mod.inter <- paste(mod$sign_, signif(mod$intercept_,2))
                 decision <- paste(" then class =", colnames(mod$confusionMatrix_)[2])
                 mod.fit <- mod[[score]]
                 if(!is.na(mod.fit)) res <- paste("|",res," ", mod.inter, decision, "| (F=",signif(mod.fit,4),sep="")
                 res <- paste(res,"|K=",mod$eval.sparsity,"|Le=",  mod$learner,"|La=", mod$language,")",sep="")
               }else
               {
                 # Bin(inter) or Ter(inter)  
                 signs <- rep("+",length(mod$coeffs_))
                 if(all(!ind.pos))
                 {
                   term.pos <- "0"
                 }else
                 {
                   term.pos <- paste("(",paste(signs[ind.pos], mod$indices_[ind.pos], sep = "", collapse = ""),")",sep="")  
                 }
                 
                 if(all(!ind.neg))
                 {
                   term.neg <- "0"
                 }else
                 {
                   term.neg <- paste("(",paste(signs[ind.neg], mod$indices_[ind.neg], sep = "", collapse = ""),")",sep="")
                 }
                 res <- paste(term.pos, " - ", term.neg)
                 mod.inter <- paste(mod$sign_, signif(mod$intercept_,2))
                 decision <- paste(" then class =",colnames(mod$confusionMatrix_)[2])
                 mod.fit <- mod[[score]]
                 if(!is.na(mod.fit)) res <- paste("|",res," ", mod.inter, decision, "| (F=",signif(mod.fit,4),sep="")
                 res <- paste(res,"|K=",mod$eval.sparsity,"|Le=",  mod$learner,"|La=", mod$language,")",sep="")
               }
             } # end exist coeffs
           }else
           {
             # SOTA
             if(!is.null(mod$coeffs_))
             {
               coeffs <- signif(as.numeric(mod$coeffs_),2)
             }else{
               coeffs <- ""
             }
             res <- mod$indices_
             res <- paste(coeffs, res, sep = " * ")
             res <- paste(res, collapse = " ")
             res <- paste(res, mod$learner, sep="|")
             mod.fit <- mod[[score]]
             if(!is.na(mod.fit)) res <- paste("|",res," ","| (F=",signif(mod.fit,4),sep="")
             res <- paste(res,"|K=",mod$eval.sparsity,"|Le=",  mod$learner,"|La=", mod$language,")",sep="")
           }
         },
         long={
           if(!isModelSota(mod))
           {
             if(all(myAssertNotNullNorNa(mod$coeffs_)))
             {
               # Bin(inter) or Ter(inter)
               ind.pos <- sign(mod$coeffs_) == 1
               ind.neg <- sign(mod$coeffs_) == -1
               
               # BTR
               if(mod$language == "ratio")
               {
                 signs <- rep("+ ",length(mod$coeffs_))
                 if(all(!ind.pos))
                 {
                   term.pos <- "0"
                 }else
                 {
                   term.pos <- paste("(",paste(signs[ind.pos], mod$names_[ind.pos], sep = "", collapse = ""),")",sep="")  
                 }
                 
                 if(all(!ind.neg))
                 {
                   term.neg <- "0"
                 }else
                 {
                   term.neg <- paste("(",paste(signs[ind.neg], mod$names_[ind.neg], sep = "", collapse = ""),")",sep="")
                 }
                 res <- paste(term.pos, "/", term.neg)
                 mod.inter <- paste(mod$sign_, signif(mod$intercept_,2))
                 decision <- paste(" then class =", colnames(mod$confusionMatrix_)[2])
                 mod.fit <- mod[[score]]
                 if(!is.na(mod.fit)) res <- paste("|",res," ", mod.inter, decision, "| (F=",signif(mod.fit,4),sep="")
                 res <- paste(res,"|K=",mod$eval.sparsity,"|Le=",  mod$learner,"|La=", mod$language,")",sep="")
               }else
               {
                 # Bin(inter) or Ter(inter)  
                 signs <- rep("+ ",length(mod$coeffs_))
                 if(all(!ind.pos))
                 {
                   term.pos <- "0"
                 }else
                 {
                   term.pos <- paste("(",paste(signs[ind.pos], mod$names_[ind.pos], sep = "", collapse = ""),")",sep="")  
                 }
                 
                 if(all(!ind.neg))
                 {
                   term.neg <- "0"
                 }else
                 {
                   term.neg <- paste("(",paste(signs[ind.neg], mod$names_[ind.neg], sep = "", collapse = ""),")",sep="")
                 }
                 res <- paste(term.pos, " - ", term.neg)
                 mod.inter <- paste(mod$sign_, signif(mod$intercept_,2))
                 decision <- paste(" then class =",colnames(mod$confusionMatrix_)[2])
                 mod.fit <- mod[[score]]
                 if(!is.na(mod.fit)) res <- paste("|",res," ", mod.inter, decision, "| (F=",signif(mod.fit,4),sep="")
                 res <- paste(res,"|K=",mod$eval.sparsity,"|L=",  mod$learner,"|La=", mod$language,")",sep="")
               }
             } # end exist coeffs
           }else
           {
             # SOTA
             if(!is.null(mod$coeffs_))
             {
               coeffs <- signif(as.numeric(mod$coeffs_),2)
             }else{
               coeffs <- ""
             }
             res <- mod$names_
             res <- paste(coeffs, res, sep = " * ")
             res <- paste(res, collapse = " ")
             res <- paste(res, mod$learner, sep="|")
             mod.fit <- mod[[score]]
             res <- paste("|",res," ","| (F=",signif(mod.fit,4),sep="")
             res <- paste(res,"|K=",mod$eval.sparsity,"|L=",  mod$learner,"|La=", mod$language,")",sep="")
           }
         },
         str={
           res <- str(mod)
         },
         {
           warning('This method does not exist! Try one of these: short, long or str')
         }
  )
  return(res)
}


#' Prints a population of model objects as text.
#'
#' @description Prints a population of model objects as text
#' @param obj: a population of models to plot
#' @param method: if "digested" a short sumary (one line) will be printed, otherwise the method will contain the 
#' specific way to print a model through the printModel() routine
#' @param score: which score to show in the fit (default:fit_)
#' @param indent: a string (default:'tab---') that will precede each element of the object.
#' @export
printPopulation <- function(obj, method = "short", score = "fit_", indent="")
{
  # sanity check
  if(!isPopulation(obj))
  {
    return(NULL)
    warning("printPopulation: the object to print is not a valid Population of models")
  }
  
  switch(method, 
         digested={
           spar <- populationGet_X("eval.sparsity")(obj)
           pop.name <- unique(spar)
           if(length(pop.name)==1)
           {
             attribute.name <- paste0("k_",pop.name)
           }else
           {
             attribute.name <- "k_mixed"
           }
           ptdf <- populationToDataFrame(obj) # convert model population to a dataframe for easy access
           #attribute.name <- ""
           attribute.value <- paste(length(obj), "models ...",
                                    paste(paste(ptdf$learner, ptdf$language, signif(ptdf$fit_,2), ptdf$eval.sparsity, sep="_")[1:min(5,length(obj))],collapse = "; "))
           cat(paste(indent,attribute.name,": ",attribute.value, "\n",sep=""))
         },
         short =, 
         str =,
         long ={
           for(i in 1:length(obj))
           {
             mod <- obj[[i]]
             print(paste(i, printModel(mod = mod, method = method, score = score), sep=":"))
           }
         },
         {
           print('printPopulation: please provide a valid method (digested/short/long/str)')
         }
  )
}


#' Prints as text the detail on a given Classifier object
#'
#' @description This function prints a summary of a Classifier object.
#' @param obj: a Classifier object
#' @param indent: a string (default:'tab---') that will precede each element of the object.
#' @return NULL if the object is not a valid Classifier
#' @export
printClassifier <- function(obj, indent="\t--- ")
{
  # sanity check
  if(!isClf(obj))
  {
    return(NULL)
    warning("printClassifier: the object to print is not a valid Classifier")
  }
  
  # Global experiment information
  if(!is.null(obj$experiment))
  {
    cat("========== Experiment ==========\n")
    for(i in 1:length(obj$experiment)){
      attribute.name <- names(obj$experiment)[i]
      attribute.value <- obj$experiment[[i]]
      cat(paste(indent,attribute.name,": ",attribute.value, "\n",sep=""))
    }
  }
  
  # Detailed learner options
  if(!is.null(obj$params))
  {
    cat("========== Learner ==========\n")
    cat(paste(indent,"learner: ", obj$learner, "\n",sep=""))
    for(i in 1:length(obj$params))
    { 
      attribute.name <- names(obj$params)[i]
      attribute.value <- obj$params[[i]]
      if(length(attribute.value)>1) # for sparsity for instance
      {
        if(!is.list(attribute.value))
        {
          cat(paste(indent, attribute.name,": ",paste(attribute.value, collapse = ","), "\n",sep=""))
        }
      }else
      {
        if(is.function(attribute.value)) # In case of functions
        {
          cat(paste(indent, attribute.name,": function","\n",sep=""))
        }else
        {
          cat(paste(indent, attribute.name,": ",attribute.value, "\n",sep=""))    
        }
      }
    }
    if(obj$learner=="metal")
    {
      nclf <- length(obj$params$list.clfs)-1
      for(i in 1:nclf)
      {
        cat(paste(indent, obj$params$list.clfs[[i]]$experiment$id,"\n",sep = ""))
      }
    }
  }
  
  # Detailed learner options
  if(isModelCollection(obj$models))
  {
    cat("========== Model Collection ==========\n")
    printModelCollection(obj$models)
  }
  
}

#' Prints as text the detail on a given Experiment object
#'
#' @description This function prints a summary of an Experiment object.
#' @param obj: an Experiment object
#' @param indent: a string (default:'tab---') that will precede each element of the object.
#' @return NULL if the object is not a valid Experiment
#' @export
printExperiment <- function(obj, indent = "\t--- ")
{
  if(!isExperiment(obj))
  {
    return(NULL)
    warning("printExperiment: the object to print is not a valid experiment.")
  }
  
  cat("========== Experiment ==========\n")
  for(i in 1:length(obj$classifier$experiment))
  {
    attribute.name <- names(obj$classifier$experiment)[i]
    attribute.value <- obj$classifier$experiment[[i]]
    cat(paste(indent,attribute.name,": ",attribute.value, "\n",sep=""))
  }
  
  if(!is.null(obj$crossVal))
  {
    cat("========== Cross validation ==========\n")
    cat(paste(indent, "total folds",": ",ncol(obj$crossVal$scores$empirical.auc), "\n",sep=""))
    cat(paste(indent, "folds",": ",ncol(obj$crossVal$scores$empirical.auc)/length(obj$classifier$params$seed), "\n",sep=""))
    cat(paste(indent, "times",": ",length(obj$classifier$params$seed), "\n",sep=""))
    cat(paste(indent, "seeds",": ",paste(obj$classifier$params$seed, collapse = ", "), "\n",sep=""))
  }else
  {
    cat("========== Cross validation ==========\n")
    cat(paste(indent, "total folds",": ",0, "\n",sep=""))
    cat(paste(indent, "folds",": ",0, "\n",sep=""))
    cat(paste(indent, "times",": ",length(obj$classifier$params$seed), "\n",sep=""))
    cat(paste(indent, "seeds",": ",paste(obj$classifier$params$seed, collapse = ", "), "\n",sep=""))
  }
  
  # Detailed learner options
  if(!is.null(obj$classifier$params))
  {
    cat("========== Learner ==========\n")
    cat(paste(indent,"learner: ", obj$classifier$learner, "\n",sep=""))
    for(i in 1:length(obj$classifier$params))
    { 
      attribute.name <- names(obj$classifier$params)[i]
      attribute.value <- obj$classifier$params[[i]]
      if(length(attribute.value) > 1) # for sparsity for instance
      {
        if(!is.list(attribute.value))
        {
          cat(paste(indent, attribute.name,": ",paste(attribute.value, collapse = ","), "\n",sep=""))
        }
      }else
      {
        if(is.function(attribute.value)) # In case of functions
        {
          cat(paste(indent, attribute.name,": function","\n",sep=""))
        }else
        {
          cat(paste(indent, attribute.name,": ",attribute.value, "\n",sep=""))    
        }
      }
    }
    if(obj$classifier$learner == "metal")
    {
      nclf <- length(obj$classifier$params$list.clfs)-1
      for(i in 1:nclf)
      {
        cat(paste(indent, obj$classifier$params$list.clfs[[i]]$experiment$id,"\n",sep = ""))
      }
    }
  }
  
  # Detailed learner options
  if(isModelCollection(obj$classifier$models))
  {
    cat("========== Model Collection ==========\n")
    printModelCollection(obj$classifier$models)
  }
}


#' Prints as text the detail on a given ModelCollection object
#'
#' @description This function prints a ModelCollection object. For each k_sparsity it will show some detail of 
#' the maximum first models
#' @param obj: a ModelCollection object
#' @param indent: a string (default:'tab---') that will precede each element of the object for the "long" method.
#' @param method: the output method (default:long) will print for each k_sparsity a short information of the population of models, 
#' while the short method will output the number of models for each k_sparsity
#' @return NULL if the object is not a valid ModelCollection.
#' @export
printModelCollection <- function(obj, indent = "\t--- ", method = "long")
{
  if(!isModelCollection(obj))
  {
    return(NULL)
    warning("printModelCollection: the object to print is not a valid experiment.")
  }
  
  switch(method, 
         short={
           paste(names(obj), unlist(lapply(obj,length)), sep=": ")
         },
         long={
           # for each k-sparsity show the number of models and some information on the first ones.
           for(i in 1:length(obj))
           {
             printPopulation(obj = obj[[i]], method = "digested", indent = indent)
           }
         },
         {
           print('printModelCollection: please provide a valid method (short/long)')
         }
  )
}


#' Prints as text the detail on a given object from the predomics package. 
#'
#' @description This function will summarize any of the predomics package objects such as can be an Experiment, 
#' a Model, a Population of models or a ModelCollection
#' @param obj: an object from the predomics object
#' @return NULL
#' @export
printy <- function(obj)
{
  type = NA
  if(isModel(obj))
  {
    type <- "model"
  }
  if(isPopulation(obj))
  {
    type <- "population"
  }
  if(isClf(obj))
  {
    type <- "classifier"
  }
  if(isExperiment(obj))
  {
    type <- "experiment"
  }
  if(isModelCollection(obj))
  {
    type <- "model.collection"
  }
  
  switch(type, 
         model={
           print(paste("Summary of Model object"))
           printModel(mod = obj, method = "long")
         },
         population={
           print(paste("Summary of a population of models with",length(obj),"models"))
           printPopulation(obj = obj[1:min(5,length(obj))], method = "long")
           if(length(obj) > 5) print("...")
         },
         model.collection={
           print(paste("Summary of a ModelCollection object with",length(obj),"populations of models"))
           printModelCollection(obj = obj, method = "long")
         },
         experiment={
           print(paste("Summary of Experiment object"))
           printExperiment(obj)
         },
         classifier={
           print(paste("Summary of Classifier object"))
           printClassifier(obj)
         },
         {
           print('printy: please provide valid predomics model')
         }
  )
  
}



#' Prints as text the detail on a given experiment along with summarized results (if computed)
#'
#' @description This function takes a population of models and makes three plots, feature prevalence in population, 
#' feature abundance by class and feature prevalence by class
#' @param pop: a population of models
#' @param X: the X dataset where to compute the abundance and prevalence
#' @param y: the target class
#' @param res_clf: the results of the classifier as well as the config object
#' @param makeplot: make a pdf file with the resulting plots (default:TRUE)
#' @param name: the suffix of the pdf file (default:"")
#' @param ord.feat: which ordering approch to use for the features (default:importance) in the models, anything 
#' else will compute automatic hierarchical ordering based on the manhattan distance
#' @param make.network: build a network and print it out in the pdf
#' @param network.layout: the network layout by default is circular (layout_in_circle) and will be a weighted Fruchterman-Reingold otherwise
#' @param network.alpha: threshold of significance for the network (default:1e-4)
#' @param verbose: print out informaiton
#' @param pdf.dims: dimensions of the pdf object (default: c(w = 25, h = 20))
#' @param filter.perc: filter by prevalence percentage in the population between 0 and 1 (default:0.05)
#' @param k_penalty: the sparsity penalty needed to select the best models of the population (default:0.75/100).
#' @param k_max: select the best population below a given threshold. If (default:0) no selection is performed.
#' @return plots if makeplot is FALSE
#' @export
analyzePopulationFeatures <- function(pop, X, y, res_clf, makeplot = TRUE, name = "", ord.feat = "importance", 
                                      make.network = TRUE, network.layout = "circular", network.alpha = 0.0001, 
                                      verbose = TRUE, pdf.dims = c(width = 25, height = 20), filter.perc = 0.05, 
                                      k_penalty = 0.75/100, 
                                      k_max = 0)
{
  pop <- selectBestPopulation(pop, p = 0.05, k_penalty = k_penalty, k_max = k_max)
  if(verbose) print(paste("There are",length(pop), "models in this population"))
  
  if(length(pop) == 1)
  {
    print("analyzePopulationFeatures: only one model after filtering. Plot can not be built... returing empty handed.")
    return(NULL)
  }
  
  pop.df <- populationToDataFrame(pop)
  pop.noz <- listOfModelsToDenseCoefMatrix(clf = res_clf$classifier, X = X, y = y, list.models = pop)
  if(verbose) print(paste("Pop noz object is created with", nrow(pop.noz), "features and", ncol(pop.noz), "models"))
  
  # clean not conform bin/bininter models
  tocheck <- pop.df$language == "bin" | pop.df$language == "bininter"
  todelete <- apply(pop.noz[,tocheck] < 0, 2, any)
  
  if(any(todelete))
  {
    if(verbose) print(paste(sum(todelete)," bin/bininter models contain negative coefficients ... deleting them"))
  }
  
  if(length(todelete) != 0)
  {
    # clean population and recompute things
    pop <- pop[!todelete]  
  }
  
  # make the feature annots
  fa <- makeFeatureAnnot(pop = pop, X = X, y = y, clf = res_clf$classifier)
  
  # filter features that are very rare in the models
  pop.noz <- filterFeaturesByPrevalence(X = fa$pop.noz, perc.prevalence = filter.perc * 100)
  
  if(!is.null(dim(pop.noz)))
  {
    if(nrow(pop.noz)==0)
    {
      pop.noz <- filterFeaturesByPrevalence(fa$pop.noz, perc.prevalence = 0)
    }  
  }
  if(verbose) print(paste("After filtering pop noz object is created with", nrow(pop.noz), "features and", ncol(pop.noz), "models"))
  
  if(verbose) print(paste("Ordering features by", ord.feat))
  switch(ord.feat,
         prevalence=
         {
           # Here we order based on the prevalence of features in the models
           prev <- getFeaturePrevalence(features = rownames(pop.noz), X = X, y = y, prop=TRUE)
           ind.feat <- order(prev$all, prev$`1`, prev$`-1`, decreasing = TRUE)
           features <- rownames(pop.noz)[ind.feat]
         },
         importance=
         {
           # Here we order based on a the crossval feature importance
           if(isExperiment(res_clf))
           {
             lr <- list(res_clf)
             names(lr) <- paste(res_clf$classifier$learner, res_clf$classifier$params$language, sep=".")
             
             feat.import <- mergeMeltImportanceCV(list.results = lr, 
                                                  filter.cv.prev = 0, 
                                                  min.kfold.nb = FALSE, 
                                                  learner.grep.pattern = "*", 
                                                  nb.top.features = NULL,
                                                  feature.selection = rownames(pop.noz),
                                                  scaled.importance = FALSE,
                                                  make.plot = TRUE)
           }
           
           if(is.null(feat.import))
           {
             print("analyzeImportanceFeatures: no feature importance data found... returning empty handed.")
             return(NULL)
           }
           
           # the most important features along with the order
           features <- rev(levels(feat.import$summary$feature))
           
         },
         hierarchical=
         {
           # Here we order based on a hierarchial clustering on the feature to model coefficient matrix
           hc.feat <- hclust(d = dist((pop.noz), method = "manhattan"), method = "ward.D")
           ind.feat <- hc.feat$order
           features <- rownames(pop.noz)[ind.feat]
         },
         {
           warning('This method does not exist !')
         }
  )
  
  if(ord.feat == "importance")
  {
    g6 <- feat.import$g
    
    g7 <- plotPrevalence(features = features, X, y)
    g7 <- g7 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    #plot abundance
    g8 <- plotAbundanceByClass(features = features, X, y)
    g8 <- g8 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    #tmp <- pop.noz; colnames(tmp) <- gsub("metal_","",colnames(pop.dense.noz))
    g9 <- plotFeatureModelCoeffs(feat.model.coeffs = pop.noz[features, ], 
                                 vertical.label = FALSE, col = c("deepskyblue1", "white", "firebrick1"))
    g9 <- g9 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    
    if(makeplot)
    {
      if(verbose) print(paste("Making plots in a dedicated pdf"))
      
      pdf(file=paste("population features",name,".pdf", sep=""), w=pdf.dims[1], h=pdf.dims[2])
      
      grid.arrange(g6, g9, g8, g7, ncol=4, widths = c(2,1,1,1))
      
      if(make.network)
      {
        if(verbose) print(paste("Making the network of co-occurance of features in the population of models"))
        #fa <- makeFeatureAnnot(pop = pop, X = X, y = y, clf = clf)
        
        try(makeFeatureModelPrevalenceNetworkCooccur(pop.noz = pop.noz, 
                                                     feature.annot = fa$feature.df[rownames(pop.noz),], 
                                                     alpha = network.alpha, 
                                                     verbose = verbose, 
                                                     layout = network.layout),
            silent = TRUE)
      }
      dev.off()
      
    }else
    {
      return(grid.arrange(g6, g9, g8, g7, ncol=4, widths = c(2,1,1,1)))
    }
  }else
  {
    g7 <- plotPrevalence(features = features, X, y)
    g7 <- g7 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    #plot abundance
    g8 <- plotAbundanceByClass(features = features, X, y)
    g8 <- g8 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    #tmp <- pop.noz; colnames(tmp) <- gsub("metal_","",colnames(pop.dense.noz))
    g9 <- plotFeatureModelCoeffs(feat.model.coeffs = pop.noz[features, ], 
                                 vertical.label = FALSE, col = c("deepskyblue1", "white", "firebrick1"))
    
    if(makeplot)
    {
      if(verbose) print(paste("Making plots in a dedicated pdf"))
      
      pdf(file=paste("population features",name,".pdf", sep=""), width = pdf.dims[1], height = pdf.dims[2])
      
      grid.arrange(g9, g8, g7, ncol=3, widths = c(2,1,1))
      
      if(make.network)
      {
        if(verbose) print(paste("Making the network of co-occurance of features in the population of models"))
        #fa <- makeFeatureAnnot(pop = pop, X = X, y = y, clf = clf)
        
        try(makeFeatureModelPrevalenceNetworkCooccur(pop.noz = pop.noz, 
                                                     feature.annot = fa$feature.df[rownames(pop.noz),], 
                                                     alpha = network.alpha, 
                                                     verbose = verbose, 
                                                     layout = network.layout),
            silent = TRUE)
      }
      dev.off()
      
    }else
    {
      return(grid.arrange(g9, g8, g7, ncol=3, widths = c(2,1,1)))
    }
  }
}


#' Prints as text the detail on a given experiment along with summarized results (if computed)
#'
#' @description This function takes a population of models and makes three plots, feature prevalence in population, 
#' feature abundance by class and feature prevalence by class
#' @param clf_res: the result of an experiment or multiple exmeriments (list of experimenets)
#' @param X: the X dataset where to compute the abundance and prevalence
#' @param y: the target class
#' @param makeplot: make a pdf file with the resulting plots (default:TRUE)
#' @param name: the suffix of the pdf file (default:"")
#' @param verbose: print out informaiton
#' @param pdf.dims: dimensions of the pdf object (default: c(w = 25, h = 20))
#' @param filter.perc: filter by prevalence percentage in the population between 0 and 1 (default:0.05)
#' @param filter.cv.prev: keep only features found in at least (default: 0.25, i.e 25 percent) of the cross validation experiments 
#' @param nb.top.features: the maximum number (default: 100) of most important features to be shown. If this value is NULL 
#' or NA, all features be returned
#' @param scaled.importance: the scaled importance is the importance multipied by the prevalence in the folds. If (default = TRUE) this will be used, the mean mda 
#' will be scaled by the prevalence of the feature in the folds and ordered subsequently 
#' @param k_penalty: the sparsity penalty needed to select the best models of the population (default:0.75/100).
#' @param k_max: select the best population below a given threshold. If (default:0) no selection is performed.
#' @return plots if makeplot is FALSE
#' @export
analyzeImportanceFeatures <- function(clf_res, X, y, makeplot = TRUE, name = "", 
                                      verbose = TRUE, pdf.dims = c(width = 25, height = 20), 
                                      filter.perc = 0.05, filter.cv.prev = 0.25, 
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
    
    # clean not conform bin/bininter models
    tocheck <- pop.df$language == "bin" | pop.df$language == "bininter"
    todelete <- apply(pop.noz[,tocheck] < 0, 2, any)
    
    if(any(todelete))
    {
      if(verbose) print(paste(sum(todelete)," bin/bininter models contain negative coefficients ... deleting them"))
    }
    
    if(length(todelete) != 0)
    {
      # clean population and recompute things
      pop <- pop[!todelete]  
    }
    
    # make the feature annots
    fa <- makeFeatureAnnot(pop = pop, X = X, y = y, clf = clf_res$classifier)
    
    # filter features that are very rare in the models
    pop.noz <- filterFeaturesByPrevalence(X = fa$pop.noz, perc.prevalence = filter.perc *100)
    
    if(nrow(pop.noz)==0)
    {
      pop.noz <- filterFeaturesByPrevalence(fa$pop.noz, perc.prevalence = 0)
    }
    if(verbose) print(paste("After filtering pop noz object is created with", nrow(pop.noz), "features and", ncol(pop.noz), "models"))
    # get the feature importance information if it exists
    lr <- list(clf_res)
    names(lr) <- paste(clf_res$classifier$learner, clf_res$classifier$params$language, sep=".")
    
    feat.import <- mergeMeltImportanceCV(list.results = lr, 
                                         filter.cv.prev = filter.cv.prev, 
                                         min.kfold.nb = FALSE, 
                                         learner.grep.pattern = "*", 
                                         nb.top.features = nb.top.features, 
                                         feature.selection = NULL,
                                         scaled.importance = scaled.importance,
                                         make.plot = TRUE)
    if(is.null(feat.import))
    {
      print("analyzeImportanceFeatures: no feature importance data found... returning empty handed.")
      return(NULL)
    }
    
    # the most important features along with the order
    features.import <- rev(levels(feat.import$summary$feature))
    
    # the most importance features found in generalization are not necessarely the same as those found in the best models, and we need to merge them to have a unified view
    pop.noz.import <- as.matrix(matrix(data = 0, nrow = length(features.import), ncol = ncol(pop.noz)))
    colnames(pop.noz.import) <- colnames(pop.noz)
    rownames(pop.noz.import) <- features.import
    features.pop.noz <- intersect(features.import, rownames(pop.noz))
    pop.noz.import[features.pop.noz, 1:ncol(pop.noz.import)] <- as.matrix(as.data.frame(pop.noz)[features.pop.noz, 1:ncol(pop.noz)])
    
    # ordering of the pop.noz.import
    ord.avail <- order(rowSums(abs(pop.noz.import)), decreasing = TRUE)
    # order on the prevalence on best models
    if(mode != "regression")
    {
      prev <- getFeaturePrevalence(features = features.import, X = X, y = y, prop=TRUE)
      ord.prev <- order(prev$all, prev$`1`, prev$`-1`, decreasing = TRUE)  
    }
    
    hc.mod <- hclust(d = dist(t(pop.noz.import), method = "manhattan"), method = "ward.D")
    
    # get the importance graphic
    g6 <- feat.import$g
    g6 <- g6 + ggtitle("")
    
    # get the prevalence graphic
    if(mode == "regression")
    {
      g7 <- plotPrevalence(features = features.import, X, y = NULL)
    }else
    {
      g7 <- plotPrevalence(features = features.import, X, y)
    }
    g7 <- g7 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    
    # get the abundance graphic
    
    g8 <- plotAbundanceByClass(features = features.import, X, y)  
    g8 <- g8 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    
    # get the feature to model graphic
    g9 <- plotFeatureModelCoeffs(feat.model.coeffs = pop.noz.import, 
                                 vertical.label = FALSE, 
                                 col = c("deepskyblue1", "white", "firebrick1"))
    
    if(makeplot)
    {
      if(verbose) print(paste("Making plots in a dedicated pdf"))
      pdf(file=paste("population features",name,".pdf", sep=""), w=pdf.dims[1], h=pdf.dims[2])
      grid.arrange(g6, 
                   #g9, 
                   g8, 
                   g7, 
                   ncol=3, 
                   widths=c(2,1,1))
      dev.off()
      
    }else
    {
      grid.arrange(g6, 
                   #g9, 
                   g8, 
                   g7, 
                   ncol=3, 
                   widths=c(2,1,1))
    }
    
  }else # multiple experiments
  {
    lr <- clf_res
    feat.import <- mergeMeltImportanceCV(list.results = lr, 
                                         filter.cv.prev = filter.cv.prev, 
                                         min.kfold.nb = FALSE, 
                                         learner.grep.pattern = "*", 
                                         nb.top.features = nb.top.features, 
                                         make.plot = TRUE)
    
    if(is.null(feat.import))
    {
      warning("analyzeImportantFeatures: These learners have no feature importance as implemented in the BTR ones. Sending an empty plot.")
      return(ggplot() + theme_void())
    }
    
    # the most important features along with the order
    features.import <- rev(levels(feat.import$summary$feature))
    
    # get the importance graphic
    g6 <- feat.import$g
    #g6 <- g6 + ggtitle("")
    # get the prevalence graphic
    if(mode == "regression")
    {
      g7 <- plotPrevalence(features = features.import, X, y = NULL)
    }else
    {
      g7 <- plotPrevalence(features = features.import, X, y)
    }
    g7 <- g7 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    # get the abundance graphic
    g8 <- plotAbundanceByClass(features = features.import, X, y)
    g8 <- g8 + theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())
    
    if(makeplot)
    {
      if(verbose) print(paste("Making plots in a dedicated pdf"))
      pdf(file=paste("population features",name,".pdf", sep=""), w=pdf.dims[1], h=pdf.dims[2])
      grid.arrange(g6, g8, g7, ncol=3, widths=c(2,1,1))
      dev.off()
      
    }else
    {
      grid.arrange(g6, g8, g7, ncol=3, widths=c(3,1,1))
    }
  } # end multiple experiments
}



#' Prints as text the detail on a given experiment along with summarized results (if computed)
#'
#' @description This function takes a population of models and creates a table with annotation on the features, 
#' such as prevalence in the models and dataset as well as different statistics
#' @param pop: a population of models
#' @param X: the X dataset where to compute the abundance and prevalence
#' @param y: the target class
#' @param clf: an object containing the different parameters of the classifier
#' @return a list with two data.frames one containing the coefficients per each model and the other a data.frame on the features
#' @export
makeFeatureAnnot <- function(pop, X, y, clf)
{
  if(!isPopulation(pop))
  {
    warning("makeFeatureAnnot: unvalid population, returning NULL")
    return(NULL)
  }
  check.X_y_w(X = X, y = y)
  
  pop.noz <- listOfModelsToDenseCoefMatrix(clf = clf, X = X, y = y, list.models = pop)
  if(any(is.na(rownames(pop.noz))))
  {
    print("makeFeatureAnnot: some features are NA in pop.noz ... omitting them")
    pop.noz <- pop.noz[!is.na(rownames(pop.noz)),]
  }
  
  # order on the prevalence on best models
  if(clf$params$objective == "cor")
  {
    prev <- getFeaturePrevalence(features = rownames(pop.noz), X = X, y = NULL, prop = TRUE)  
  }else
  {
    prev <- getFeaturePrevalence(features = rownames(pop.noz), X = X, y = y, prop = TRUE)
  }
  
  # a) prevalence in each class and in X
  feat.preval.X <- as.data.frame(prev); colnames(feat.preval.X) <- paste("prev.prop",names(prev))
  # b) enrichment analyses based on prevalence
  if(clf$params$objective == "cor")
  {
    prev.enrichment <- plotPrevalence(features = rownames(pop.noz), X, y = NULL, plot=FALSE)
  }else
  {
    prev.enrichment <- plotPrevalence(features = rownames(pop.noz), X, y, plot=FALSE)
  }
  
  feat.preval.chisq <- as.data.frame(prev.enrichment[c("chisq.p", "chisq.q")])
  # c) prevalence in models
  mod.prev <- rowSums(abs(pop.noz))
  # d) mean coefficient in real models
  pop.noz.na <- pop.noz # in order not to take into account the zeros we transform them to NA
  pop.noz.na[pop.noz == 0] <- NA
  coeff.mean <- rowMeans(pop.noz.na, na.rm = TRUE)
  # f) differential abundance analyses
  suppressWarnings(feat.abund.wilcox <- plotAbundanceByClass(features = rownames(pop.noz), X, y, plot=FALSE)[,c("p","q","status")])
  colnames(feat.abund.wilcox) <- c("wilcox.p", "wilcox.q", "wilcox.class")
  if(sum(dim(feat.preval.chisq)) == 0)
  {
    feat.preval.chisq <- rep(NA, nrow(feat.abund.wilcox))
  }
  # put everything together
  feature.df <- data.frame(feat.preval.X, 
                           feat.preval.chisq, 
                           mod.prev, 
                           feat.abund.wilcox,
                           coeff.mean,
                           check.names = FALSE)
  return(list(pop.noz = pop.noz,
              feature.df = feature.df))
}


#' Prints as text the detail on a given experiment along with summarized results (if computed)
#'
#' @description This function will use the coocur package to compute the co-occurance of features in a population of models
#' @param pop.noz: a data.frame of in features in the rows and models in the columns. 
#' This table contains the feature coefficients in the models and is obtained by makeFeatureAnnot()
#' @param feature.annot: a data frame with annotation on features obtained by makeFeatureAnnot()
#' @param alpha: the significane p-value of the co-occurance.
#' @param verbose: print out information during run
#' @param layout: the network layout by default is circular (layout_in_circle) and will be a weighted Fruchterman-Reingold otherwise
#' @return plots a graph
#' @export
makeFeatureModelPrevalenceNetworkCooccur <- function(pop.noz, feature.annot, alpha = 0.05, verbose = TRUE, layout = "circlular")
{
  require(igraph)
  require(cooccur)
  
  if(verbose) print("Creating the cooccur object")
  print(system.time(cooccur.pop.noz <- cooccur(mat = as.data.frame(pop.noz !=0), 
                                               type = "spp_site", 
                                               thresh = "TRUE", 
                                               spp_names = TRUE
  )))
  
  #prob.table(cooccur.pop.noz)
  
  edges <- data.frame(from = cooccur.pop.noz$results$sp1_name,
                      to = cooccur.pop.noz$results$sp2_name,
                      cooccur.pop.noz$results
  )
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  
  #-----
  if(verbose) print(paste(nrow(edges),"edges are found"))
  rownames(edges) <- paste(edges$from, edges$to, sep=" => ")
  
  
  edges$sign <- rep(0, nrow(edges))  
  edges$sign[edges$p_lt < alpha] <- -1
  edges$sign[edges$p_gt < alpha] <- 1
  
  # delete the non significant edges
  edges <- edges[edges$sign != 0,]
  
  if(verbose) print(paste(nrow(edges),"edges are kept after filtering by p-value threshold", alpha))
  
  
  # BUILD NETWORK
  #-------------------------------------------------------------------------------------------
  # ANNOTATION of the edges
  allnodes <- unique(c(as.character(edges$from), as.character(edges$to)))
  edges.annot <- feature.annot[allnodes,] 
  edges.annot <- data.frame(rownames(edges.annot),edges.annot);
  colnames(edges.annot)[match("name",colnames(edges.annot))] <- "name_long"; 
  colnames(edges.annot)[1] <- "name"
  
  # create the igraph object
  gD <- graph.data.frame(d = edges, directed = TRUE, 
                         vertices = edges.annot)
  
  if(verbose) print(paste("The network is built"))
  
  # Calculate degree for all nodes
  degAll <- igraph::degree(gD, v = V(gD), mode = "all")
  gD <- set.vertex.attribute(gD, "degree", index = V(gD)$name, value = degAll)
  # Calculate betweenness for all nodes
  betAll <- igraph::betweenness(gD, v = V(gD), directed = FALSE) / (((vcount(gD) - 1) * (vcount(gD)-2)) / 2)
  betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll)); rm(betAll)
  # Add new node/edge attributes based on the calculated node properties/similarities
  gD <- set.vertex.attribute(gD, "betweenness", index = V(gD)$name, value = betAll.norm)
  
  # Calculate edge properties and add to the network
  E(gD)$color <- c("#DC143C","#A6A6A6")[as.factor(factor(E(gD)$sign, levels=c('-1','1')))]
  
  #Calculate Dice similarities between all pairs of nodes
  dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
  # The following function will transform a square matrix to an edge driven one and add values to each edge
  F1 <- function(x) {data.frame(dice = dsAll[which(V(gD)$name == as.character(x$from)), which(V(gD)$name == as.character(x$to))])}
  library(plyr)
  edges.ext <- ddply(edges, .variables=c("from", "to"), function(x) data.frame(F1(x))); dim(edges.ext)
  
  gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
  E(gD)[as.character(edges.ext$from) %--% as.character(edges.ext$to)]$similarity <- as.numeric(edges.ext$dice)
  
  # Check the attributes
  summary(gD)
  #l <- layout.fruchterman.reingold(gD)
  
  if(layout == "circular")
  {
    l <- layout_in_circle(gD, order = V(gD))  
  }else
  {
    l <- layout_with_fr(gD, weights=E(gD)$weight)
  }
  
  
  plot(gD,
       vertex.label = V(gD)$name_long,
       vertex.color = c("deepskyblue1", "firebrick1")[factor(V(gD)$wilcox.class)],
       vertex.size=V(gD)$mod.prev/max(V(gD)$mod.prev)*10 + 1,
       edge.arrow.size=0,
       asp=TRUE,
       rescale=TRUE,
       layout=l,
       #edge.arrow.mode = E(gD)$infOrt,
       vertex.label.cex = 0.7,
       vertex.label.dist=0,
       edge.width = E(gD)$prob_cooccur*10
  )
  
  #return(gD)
}


#' Prints as text the detail on a given experiment along with summarized results (if computed)
#'
#' @description This function will use the miic package to compute the co-occurance of features in a population of models
#' @param pop.noz: a data.frame of in features in the rows and models in the columns. 
#' This table contains the feature coefficients in the models and is obtained by makeFeatureAnnot()
#' @param feature.annot: a data frame with annotation on features obtained by makeFeatureAnnot()
#' @param cor.th: a threshold abtained on the partial correlation value
#' @param verbose: print out information during run
#' @param layout: the network layout by default is circular (layout_in_circle) and will be a weighted Fruchterman-Reingold otherwise
#' @return plots a graph
#' @export
makeFeatureModelPrevalenceNetworkMiic <- function(pop.noz, feature.annot, cor.th = 0.3, verbose = TRUE, layout = "circlular")
{
  require(igraph)
  require(miic)
  
  if(verbose) print("Creating the miic object")
  mobj <- miic(inputData = as.data.frame(t(pop.noz)))
  #miic.plot(mobj, igraphLayout = igraph::layout.fruchterman.reingold)
  
  #-----
  # load the edge information for spectral3off2 network
  edges <- mobj$retained.edges.summary
  if(verbose) print(paste(nrow(edges),"edges are found"))
  #dim(edges) # 350 edges
  colnames(edges)[1:2] <- c("from","to")
  rownames(edges) <- paste(edges$from, edges$to, sep=" => ")
  
  # clean those edges that have no link
  edges <- edges[abs(edges$partial_correlation) > cor.th,]
  if(verbose) print(paste(nrow(edges),"edges are kept after filtering by absolute correlation threshold", cor.th))
  
  
  # BUILD NETWORK
  #-------------------------------------------------------------------------------------------
  # ANNOTATION of the edges
  allnodes <- unique(c(edges$from,edges$to))
  edges.annot <- feature.annot[allnodes,] 
  edges.annot <- data.frame(rownames(edges.annot),edges.annot);
  colnames(edges.annot)[match("name",colnames(edges.annot))] <- "name_long"; 
  colnames(edges.annot)[1] <- "name"
  
  # create the igraph object
  gD <- graph.data.frame(d = edges, directed = TRUE, 
                         vertices = edges.annot)
  
  if(verbose) print(paste("The network is built"))
  
  # Calculate degree for all nodes
  degAll <- igraph::degree(gD, v = V(gD), mode = "all")
  gD <- set.vertex.attribute(gD, "degree", index = V(gD)$name, value = degAll)
  # Calculate betweenness for all nodes
  betAll <- igraph::betweenness(gD, v = V(gD), directed = FALSE) / (((vcount(gD) - 1) * (vcount(gD)-2)) / 2)
  betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll)); rm(betAll)
  # Add new node/edge attributes based on the calculated node properties/similarities
  gD <- set.vertex.attribute(gD, "betweenness", index = V(gD)$name, value = betAll.norm)
  
  # Calculate edge properties and add to the network
  E(gD)$color <- c("#DC143C","#A6A6A6")[as.factor(factor(sign(E(gD)$infOrt), levels=c('-1','1')))]
  E(gD)$infOrt[E(gD)$infOrt==  1] <- 0; 
  E(gD)$infOrt[E(gD)$infOrt== -2] <- 1
  
  #Calculate Dice similarities between all pairs of nodes
  dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
  # The following function will transform a square matrix to an edge driven one and add values to each edge
  F1 <- function(x) {data.frame(dice = dsAll[which(V(gD)$name == as.character(x$from)), which(V(gD)$name == as.character(x$to))])}
  library(plyr)
  edges.ext <- ddply(edges, .variables=c("from", "to"), function(x) data.frame(F1(x))); dim(edges.ext)
  
  gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
  E(gD)[as.character(edges.ext$from) %--% as.character(edges.ext$to)]$similarity <- as.numeric(edges.ext$dice)
  
  # Check the attributes
  summary(gD)
  
  if(layout == "circular")
  {
    l <- layout_in_circle(gD, order = V(gD))  
  }else
  {
    l <- layout_with_fr(gD, weights=E(gD)$weight)
  }
  
  plot(gD,
       vertex.label = V(gD)$name_long,
       vertex.color = c("deepskyblue1", "firebrick1")[factor(V(gD)$wilcox.class)],
       vertex.size=V(gD)$mod.prev/max(V(gD)$mod.prev)*10 + 3,
       edge.arrow.size=.2,
       asp=TRUE,
       rescale=TRUE,
       layout=l,
       edge.arrow.mode = E(gD)$infOrt,
       vertex.label.cex = 0.7,
       vertex.label.dist=0,
       edge.width = log(E(gD)$log_confidence)
  )
  
  #return(gD)
}

