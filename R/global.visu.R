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
# @date: October 2023
# @date: November 2024


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

#' Plot Comparative Empirical Score for Multiple Methods
#'
#' This function generates a plot comparing the empirical scores (such as AUC or
#' accuracy) across multiple methods for different values of the sparse
#' parameter (k sparse). The plot includes lines representing each method and
#' points indicating the empirical score at each k sparse value. Horizontal
#' lines indicate major thresholds (e.g., AUC = 0.5 or accuracy = majority
#' class).
#'
#' @param digested.results A list containing the empirical results of the
#'   models, including performance scores for various methods.
#' @param ylim A numeric vector of length 2 specifying the limits for the
#'   y-axis. Default is `c(0.5, 1)`.
#' @param score A string specifying which score to visualize, e.g., "auc_" or
#'   "accuracy_". Default is `"auc_"`.
#' @param main A string specifying the title of the plot. Default is an empty
#'   string.
#'
#' @details The function plots empirical scores (such as AUC or accuracy) for
#' different methods across various values of k sparse. It handles multiple
#' methods and adds horizontal lines to indicate important thresholds, such as
#' AUC = 0.5 or the majority class in classification tasks.
#'
#' The plot is created using \code{ggplot2}, and different methods can be
#' assigned different colors and point shapes. If no empirical data is provided,
#' a blank plot with axis labels and horizontal lines for thresholds is
#' displayed.
#'
#' @return A ggplot object visualizing the comparative empirical scores across
#'   multiple methods.
#'
#' @examples
#' # Assuming digested.results contains the performance scores for methods
#' plotComparativeEmpiricalScore(digested.results, ylim = c(0.5, 1), score = "auc_", main = "Comparison of AUC across Methods")
#'
#' @author Edi Prifti (IRD)
#'
#' @import ggplot2
#' @import RColorBrewer
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


#' Plot Comparative Cross-Validation (CV) Performance for Multiple Methods
#'
#' This function generates a plot comparing the cross-validation (CV)
#' performance (such as AUC, accuracy, or other scores) across multiple methods
#' for different values of the sparse parameter (k sparse). The plot includes
#' lines representing each method and points indicating the empirical or
#' generalization performance score at each k sparse value. Optional error bars
#' (confidence intervals) can be included in the plot.
#'
#' @param digested.results A list containing the CV results of the models,
#'   including performance scores for various methods.
#' @param ylim A numeric vector of length 2 specifying the limits for the
#'   y-axis. Default is `c(0.5, 1)`.
#' @param generalization A logical value (`TRUE` or `FALSE`). If `TRUE`, the
#'   plot shows the generalization performance (cross-validation performance
#'   across folds). If `FALSE`, it shows empirical performance. Default is
#'   `TRUE`.
#' @param score A string specifying which score to visualize, e.g., "auc_",
#'   "accuracy_", "recall_", etc. Default is `"auc_"`.
#' @param ci A logical value (`TRUE` or `FALSE`). If `TRUE`, confidence
#'   intervals (error bars) are shown in the plot. Default is `TRUE`.
#' @param main A string specifying the title of the plot. Default is an empty
#'   string.
#'
#' @details The function plots cross-validation (CV) scores (such as AUC,
#' accuracy, etc.) for different methods across various values of k sparse. It
#' handles both generalization (cross-validation) and empirical tasks, and
#' includes optional error bars representing confidence intervals for each
#' score.
#'
#' The plot is created using \code{ggplot2}, and different methods can be
#' assigned different colors and point shapes. Horizontal lines indicate
#' important thresholds, such as AUC = 0.5 or the majority class in
#' classification tasks.
#'
#' @return A ggplot object visualizing the comparative cross-validation (CV)
#'   scores across multiple methods.
#'
#' @examples
#' # Assuming digested.results contains the cross-validation performance scores for methods
#' plotComparativeCV(digested.results, ylim = c(0.5, 1), score = "auc_", ci = TRUE, main = "Comparison of AUC across Methods")
#'
#' # You can customize the plot by adjusting the score, error bars (ci), and other parameters.
#'
#' @author Edi Prifti (IRD)
#'
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom stats sd
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


#' Plot Comparative Best Cross-Validation (CV) Performance for Multiple Methods
#'
#' This function generates a plot comparing the best cross-validation (CV)
#' performance (such as AUC, accuracy, or other scores) across multiple methods
#' for different values of the sparse parameter (k sparse). The plot includes
#' lines representing each method and points indicating the best CV performance
#' score at each k sparse value. Optional error bars (confidence intervals) can
#' be included in the plot. The best method for each k sparse is also indicated
#' in the legend.
#'
#' @param digested.results A list containing the best CV results of the models,
#'   including performance scores for various methods.
#' @param ylim A numeric vector of length 2 specifying the limits for the
#'   y-axis. Default is `c(0.5, 1)`.
#' @param generalization A logical value (`TRUE` or `FALSE`). If `TRUE`, the
#'   plot shows the generalization performance (cross-validation performance
#'   across folds). If `FALSE`, it shows empirical performance. Default is
#'   `TRUE`.
#' @param score A string specifying which score to visualize, e.g., "auc_",
#'   "accuracy_", "recall_", etc. Default is `"auc_"`.
#' @param ci A logical value (`TRUE` or `FALSE`). If `TRUE`, confidence
#'   intervals (error bars) are shown in the plot. Default is `TRUE`.
#' @param main A string specifying the title of the plot. Default is an empty
#'   string.
#'
#' @details The function plots the best cross-validation (CV) scores (such as
#' AUC, accuracy, etc.) for different methods across various values of k sparse.
#' It handles both generalization (cross-validation) and empirical tasks, and
#' includes optional error bars representing confidence intervals for each
#' score.
#'
#' The plot is created using \code{ggplot2}, and different methods can be
#' assigned different colors and point shapes. Horizontal lines indicate
#' important thresholds, such as AUC = 0.5 or the majority class in
#' classification tasks.
#'
#' The plot also includes a legend with the best method for each k sparse value.
#'
#' @return A ggplot object visualizing the best cross-validation (CV) scores
#'   across multiple methods.
#'
#' @examples
#' # Assuming digested.results contains the best cross-validation performance scores for methods
#' plotComparativeBestCV(digested.results, ylim = c(0.5, 1), score = "auc_", ci = TRUE, main = "Comparison of Best AUC across Methods")
#'
#' # You can customize the plot by adjusting the score, error bars (ci), and other parameters.
#'
#' @author Edi Prifti (IRD)
#'
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom stats sd
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


#' Plot Comparative Results for Multiple Methods and Cross-Validation Scores
#'
#' This function generates plots comparing multiple performance metrics (such as
#' AUC, accuracy, recall, precision, F1-score, etc.) across different methods.
#' It can display results from both empirical and cross-validation (CV)
#' evaluations, with options to show the best results across methods or by
#' k-sparsity values. The function can handle both classification and regression
#' tasks and supports visualizing both the empirical and generalization
#' performance.
#'
#' @param digested.results A list containing the performance results, including
#'   both empirical and cross-validation (CV) scores for various methods.
#' @param plot A logical value (`TRUE` or `FALSE`). If `TRUE`, the function will
#'   generate and display the plots. Default is `TRUE`.
#' @param ylim A numeric vector of length 2 specifying the limits for the
#'   y-axis. Default is `c(0.5, 1)`.
#' @param best A logical value (`TRUE` or `FALSE`). If `TRUE`, the function will
#'   plot the best results across methods, regardless of the k-sparsity. Default
#'   is `FALSE`.
#' @param ci A logical value (`TRUE` or `FALSE`). If `TRUE`, confidence
#'   intervals (error bars) will be shown in the plots. Default is `FALSE`.
#' @param main A string specifying the title of the plots. Default is an empty
#'   string.
#' @param mode A string specifying the type of model being analyzed. Options are
#'   `"classification"` or `"regression"`. Default is `"classification"`.
#'
#' @details The function generates multiple plots comparing performance metrics
#'   such as AUC, accuracy, recall, precision, F1-score, and correlation, across
#'   multiple methods. The plots can show:
#' - Empirical performance for each method.
#' - Cross-validation performance (generalization) for each method.
#' - The best results across methods, either by k-sparsity or regardless of k-sparsity.
#'
#'   The plots are generated using \code{ggplot2} and arranged in a grid using
#'   the \code{multiplot} function. The user can choose to visualize the results
#'   for classification or regression models.
#'
#' @return If `plot = TRUE`, the function displays the plots. If `plot = FALSE`,
#'   the function returns a list of ggplot objects for further manipulation.
#'
#' @examples
#' # Assuming digested.results contains the performance scores for methods
#' plotComparativeResults(digested.results, plot = TRUE, ylim = c(0.5, 1),
#' best = TRUE, ci = TRUE, main = "Comparison of Results")
#'
#' @author Edi Prifti (IRD)
#'
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom gridExtra grid.arrange
#' @importFrom stats sd
#' @export
plotComparativeResults <- function(digested.results, 
                                   plot = TRUE, 
                                   ylim = c(0.5, 1), 
                                   best = FALSE, 
                                   ci = FALSE, 
                                   main = "", 
                                   mode = "classification")
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


#' Plot Comparative Results for Best Performance Across Multiple Methods
#'
#' This function generates plots comparing the best performance metrics (such as
#' AUC, accuracy, recall, precision, F1-score, etc.) across multiple methods. It
#' visualizes both empirical performance and cross-validation (CV) scores for
#' different methods, with the option to plot results for different scores. The
#' function can handle both classification and regression tasks.
#'
#' @param digested.results A list containing the performance results, including
#'   both empirical and cross-validation (CV) scores for various methods.
#' @param plot A logical value (`TRUE` or `FALSE`). If `TRUE`, the function will
#'   generate and display the plots. Default is `TRUE`.
#' @param ylim A numeric vector of length 2 specifying the limits for the
#'   y-axis. Default is `c(0.5, 1)`.
#'
#' @details The function generates multiple plots comparing the best performance
#' metrics such as AUC, accuracy, recall, precision, F1-score, and correlation,
#' across multiple methods. The plots include:
#' - Empirical performance for each method.
#' - Cross-validation performance (generalization) for each method.
#'
#' The function can visualize both classification and regression models, and
#' displays the best performance across different methods. The plots are
#' arranged using the `multiplot` function.
#'
#' @return If `plot = TRUE`, the function displays the plots. If `plot = FALSE`,
#' the function returns a list of ggplot objects for further manipulation.
#'
#' @examples
#' # Assuming digested.results contains the performance scores for methods
#' plotComparativeResultsBest(digested.results, plot = TRUE, ylim = c(0.5, 1))
#'
#' # You can customize the plot by adjusting the score, error bars (ci), and other parameters.
#'
#' @author Edi Prifti (IRD)
#'
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom gridExtra grid.arrange
#' @importFrom stats sd
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


#' Plot MGS Quality for Genes
#'
#' This function visualizes the quality of Multi-Genome Scoring (MGS) for a
#' dataset. It generates multiple plots, including barcode plots and individual
#' signal metric plots for the first 50 genes (or fewer, depending on the
#' dataset size). The function also computes and visualizes various signal
#' metrics, with results displayed in color-coded plots based on the values of
#' the computed metrics.
#'
#' @param dat A data frame or matrix of gene data. Rows represent genes and
#'   columns represent individuals.
#' @param main A string for the title of the plots. Default is `"mgs"`.
#' @param return.scores A logical value (`TRUE` or `FALSE`). If `TRUE`, the
#'   function returns the computed signal scores for the genes. Default is
#'   `TRUE`.
#'
#' @details The function generates a series of plots to assess the quality of
#' MGS for a dataset. The process includes:
#' - Creating a barcode plot for the first 50 genes.
#' - Generating individual plots for each of the computed signal metrics, with the colors of the points indicating the magnitude of the signal (from black to red).
#'
#' The function calls `plotBarcode` for visualizing the dataset and
#' `computeSignalMetrics` for calculating the signal metrics.
#'
#' If the dataset contains fewer than 50 genes, all available genes are used
#' instead. The signal metrics are calculated and displayed in individual plots,
#' with color gradients indicating the score value.
#'
#' @return If `return.scores = TRUE`, the function returns a data frame
#' containing the computed signal metrics for the dataset.
#'
#' @examples
#' # Assuming `dat` is a data frame of gene data
#' plotMGSQuality(dat)
#'
#' # To get the computed signal metrics
#' scores <- plotMGSQuality(dat, return.scores = TRUE)
#'
#' @author Edi Prifti (IRD)
#'
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
#' Create Multiple Plots on One Page
#'
#' This function arranges multiple plots on a single page using grid layout. It
#' allows for custom arrangement of the plots and can save the plots to a file
#' if needed. It supports both direct plotting and generating a multi-plot
#' layout from a list of plots.
#'
#' @param ... One or more ggplot objects to be displayed.
#' @param plotlist A list of ggplot objects to be displayed. This is an
#'   alternative to passing the plots as `...`.
#' @param file Optional; a character string specifying the file path to save the
#'   plot as a file (e.g., PNG, PDF). Default is `NULL`, which means the plot is
#'   shown in the R graphics window.
#' @param cols The number of columns to arrange the plots in. Default is `1`.
#'   This is used to calculate the layout if `layout` is not provided.
#' @param layout A matrix specifying the layout of the plots on the page. If
#'   `NULL`, the layout is automatically calculated based on the number of plots
#'   and the `cols` parameter.
#'
#' @details The function arranges multiple ggplot objects in a grid layout, with
#' the number of columns determined by the `cols` argument. The function will
#' automatically adjust the number of rows to fit all the plots. If `layout` is
#' provided, it will override the `cols` argument to control the layout.
#'
#' If `file` is provided, the function will save the multi-plot layout to the
#' specified file. The supported formats depend on the device used (e.g., PNG,
#' PDF).
#'
#' @return No return value. The function displays the plots or saves them to a
#' file if `file` is specified.
#'
#' @examples
#' library(ggplot2)
#' p1 <- ggplot(mtcars, aes(mpg, disp)) + geom_point()
#' p2 <- ggplot(mtcars, aes(mpg, hp)) + geom_point()
#' p3 <- ggplot(mtcars, aes(disp, hp)) + geom_point()
#'
#' # Display the plots in a 2x2 grid
#' multiplot(p1, p2, p3, cols=2)
#'
#' # Save the plots to a file
#' multiplot(p1, p2, p3, file="my_plots.pdf", cols=2)
#'
#' @author Edi Prifti (IRD)
#'
#' @import grid
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


#' Plot Feature Prevalence and Enrichment
#'
#' This function visualizes the prevalence of features across different groups
#' in a dataset and computes feature enrichment, providing an optional
#' statistical test for enrichment. The plot displays the prevalence of features
#' in each group, and if enrichment data is available, it adds significance
#' markers.
#'
#' @param features A character vector of feature names to be plotted.
#' @param X A data matrix or data frame where each row is an observation and
#'   each column is a feature.
#' @param y A vector of class labels (e.g., 1 and -1 for binary classification)
#'   corresponding to the rows in `X`.
#' @param topdown Logical; whether to arrange the features in a top-down order
#'   (default is `TRUE`).
#' @param main A string for the title of the plot.
#' @param plot Logical; if `TRUE`, the function will display the plot. If
#'   `FALSE`, the function will return the enrichment statistics.
#' @param col.pt Colors for points in the plot (default is `c("deepskyblue4",
#'   "firebrick4")`).
#' @param col.bg Colors for bars in the plot (default is `c("deepskyblue1",
#'   "firebrick1")`).
#' @param zero.value The value to replace in `y` for missing values (default is
#'   `0`).
#'
#' @details The function computes and visualizes the prevalence of the specified
#' features across different groups in the dataset, showing the percentage of
#' occurrences in each class. If statistical enrichment tests are available, the
#' function will display these using asterisks for significant features. The
#' enrichment is computed using chi-square tests.
#'
#' @return If `plot = TRUE`, the function returns a ggplot object displaying the
#' feature prevalence and enrichment. If `plot = FALSE`, the function returns
#' the enrichment results.
#'
#' @examples
#' # Example usage
#' features <- c("feature1", "feature2", "feature3")
#' X <- data.frame(feature1 = rnorm(100), feature2 = rnorm(100), feature3 = rnorm(100))
#' y <- sample(c(1, -1), 100, replace = TRUE)
#'
#' # Plot feature prevalence
#' plotPrevalence(features, X, y, main = "Feature Prevalence Plot")
#'
#' # Get enrichment statistics without plotting
#' plotPrevalence(features, X, y, plot = FALSE)
#'
#' @author Edi Prifti (IRD)
#'
#' @import ggplot2
#' @importFrom reshape2 melt
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



#' Plot Feature Abundance by Class
#'
#' This function visualizes the abundance of features across different classes
#' in a dataset. It creates boxplots to show the distribution of feature
#' abundances in each class, along with statistical tests for the significance
#' of differences between the classes. The function supports both classification
#' and regression tasks.
#'
#' @param features A character vector of feature names to be plotted.
#' @param X A data matrix or data frame where each row represents an observation
#'   and each column represents a feature.
#' @param y A vector of class labels (e.g., 1 and -1 for binary classification,
#'   or continuous for regression) corresponding to the rows in `X`.
#' @param topdown Logical; whether to arrange the features in a top-down order
#'   (default is `TRUE`).
#' @param main A string for the title of the plot.
#' @param plot Logical; if `TRUE`, the function will display the plot. If
#'   `FALSE`, the function will return the statistical results of the test.
#' @param col.pt Colors for points in the plot (default is `c("deepskyblue4",
#'   "firebrick4")`).
#' @param col.bg Colors for boxplots in the plot (default is `c("deepskyblue1",
#'   "firebrick1")`).
#'
#' @details This function computes and visualizes the abundance of features in
#' each class (group) of the dataset. It creates a boxplot for each feature and
#' computes a non-parametric test for differences in abundance between classes.
#' The function supports both classification (e.g., binary or multi-class) and
#' regression tasks (using continuous values in `y`).
#'
#' In classification mode, the plot compares the two classes (e.g., 1 vs -1 for
#' binary classification) and adds significance markers (e.g., asterisks) for
#' features with significant differences in abundance between classes. In
#' regression mode, it compares feature abundance across all observations and
#' computes correlations with the response variable.
#'
#' @return If `plot = TRUE`, the function returns a ggplot object displaying the
#' feature abundance by class. If `plot = FALSE`, it returns the statistical
#' results of the test (p-values and q-values).
#'
#' @examples
#' # Example usage for classification task
#' features <- c("feature1", "feature2", "feature3")
#' X <- data.frame(feature1 = rnorm(100), feature2 = rnorm(100), feature3 = rnorm(100))
#' y <- sample(c(1, -1), 100, replace = TRUE)
#'
#' # Plot feature abundance
#' plotAbundanceByClass(features, X, y, main = "Feature Abundance Plot")
#'
#' # Get statistical results without plotting
#' plotAbundanceByClass(features, X, y, plot = FALSE)
#'
#' @author Edi Prifti (IRD)
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
plotAbundanceByClass <- function(features, X, y, topdown = TRUE, 
                                 main = "", plot = TRUE, 
                                 col.pt = c("deepskyblue4", "firebrick4"), 
                                 col.bg = c("deepskyblue1", "firebrick1"))
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



#' Plot Feature Model Coefficients
#'
#' This function visualizes the coefficients of features across different models
#' in a heatmap-style plot. It allows for customized ordering of the features
#' and models, and it supports vertical or horizontal axis labels.
#'
#' @param feat.model.coeffs A matrix or data frame where rows represent
#'   features, and columns represent models. The values in the matrix correspond
#'   to the coefficients of the features in each model.
#' @param topdown Logical; if `TRUE`, the features will be arranged from top to
#'   bottom in the plot. If `FALSE`, the original order is preserved.
#' @param main A string for the title of the plot.
#' @param col A vector of colors for the heatmap. Default is `c("deepskyblue1",
#'   "white", "firebrick1")`.
#' @param vertical.label Logical; if `TRUE`, the feature labels on the y-axis
#'   will be rotated vertically.
#'
#' @details The function generates a heatmap-style plot of feature coefficients
#' across multiple models. The color intensity represents the coefficient
#' values, with positive values shown in one color, negative values in another,
#' and zeros in a neutral color. The function supports customization of the plot
#' layout and label orientation.
#'
#' The `topdown` argument controls the ordering of the features in the plot.
#' When set to `TRUE`, the features are arranged in reverse order, and if set to
#' `FALSE`, the original order is maintained. The function also allows the user
#' to rotate the feature labels vertically if needed.
#'
#' @return A `ggplot` object displaying the heatmap of feature coefficients
#' across models.
#'
#' @examples
#' # Example usage
#' features <- c("feature1", "feature2", "feature3")
#' models <- c("model1", "model2", "model3")
#' coeffs <- matrix(runif(9), nrow = 3, ncol = 3)
#' rownames(coeffs) <- features
#' colnames(coeffs) <- models
#'
#' # Plot feature model coefficients
#' plotFeatureModelCoeffs(coeffs, main = "Feature Coefficients Heatmap")
#'
#' @author Edi Prifti (IRD)
#'
#' @import ggplot2
#' @importFrom reshape2 melt
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



#' Plot Population of Models
#'
#' This function visualizes a population of models by plotting their feature
#' importance or coefficients. It supports sorting of features by discriminance
#' or importance, and allows for customization of colors and the number of
#' columns in the layout.
#'
#' @param pop A population of models (either a list of models or a single
#'   model).
#' @param X A data matrix or data frame where rows represent observations and
#'   columns represent features.
#' @param y A vector of class labels (e.g., 1 and -1 for binary classification)
#'   corresponding to the rows in `X`.
#' @param sort.features Logical; if `TRUE`, the features will be sorted by their
#'   discriminance with respect to `y`. Default is `FALSE`.
#' @param sort.ind A vector of indices for sorting the features. If `NULL`, the
#'   function will determine the order based on discriminance. Default is
#'   `NULL`.
#' @param col.sign A vector of two colors (default is `c("deepskyblue1",
#'   "firebrick1")`) used to represent positive and negative coefficients,
#'   respectively.
#' @param ncol The number of columns to arrange the plots in (default is `10`).
#' @param slim Logical; if `TRUE`, the plots will be simplified. Default is
#'   `FALSE`.
#' @param importance Logical; if `TRUE`, the feature importance will be plotted.
#'   Default is `FALSE`.
#'
#' @details This function generates a series of plots for a population of
#' models, displaying the feature coefficients or importances. If the population
#' consists of multiple models, each model's coefficients or importances are
#' plotted in a separate subplot. The function supports feature sorting based on
#' their discriminance or importance, and customizes the layout of the plots.
#'
#' @return If the population consists of a single model, a plot of the model's
#' feature coefficients or importance is returned. If the population consists of
#' multiple models, a grid of plots is displayed using `grid.arrange`.
#'
#' @examples
#' # Example usage for a population of models
#' pop <- list(model1, model2, model3)  # Assume these are pre-defined models
#' X <- data.frame(feature1 = rnorm(100), feature2 = rnorm(100))
#' y <- sample(c(1, -1), 100, replace = TRUE)
#'
#' # Plot the population
#' plotPopulation(pop, X, y, sort.features = TRUE, ncol = 3)
#'
#' @author Edi Prifti (IRD)
#'
#' @import gridExtra
#' @import ggplot2
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

#' Plot Model Coefficients and Importance
#'
#' This function visualizes the coefficients of a model, with the option to plot
#' feature importance. It supports various model types, including SOTA models
#' and Random Forests, and can display the results in a variety of layouts,
#' including plots of feature coefficients, feature importance, or decision
#' trees.
#'
#' @param mod A model object (e.g., a fitted machine learning model such as SVM,
#'   logistic regression, or random forest).
#' @param X A data matrix or data frame where each row is an observation and
#'   each column is a feature.
#' @param y A vector of class labels (e.g., 1 and -1 for binary classification,
#'   or continuous for regression) corresponding to the rows in `X`.
#' @param sort.features Logical; if `TRUE`, the features will be sorted by their
#'   discriminance with respect to `y`. Default is `FALSE`.
#' @param sort.ind A vector of indices for sorting the features. If `NULL`, the
#'   function will determine the order based on discriminance. Default is
#'   `NULL`.
#' @param feature.name Logical; if `TRUE`, the feature names will be displayed
#'   on the plot.
#' @param col.sign A vector of two colors for positive and negative coefficients
#'   (default is `c("deepskyblue1", "firebrick1")`).
#' @param main A string for the title of the plot.
#' @param slim Logical; if `TRUE`, the plot will be simplified (e.g., removing
#'   axis labels and ticks). Default is `FALSE`.
#' @param importance Logical; if `TRUE`, the feature importance will be plotted.
#'   Default is `FALSE`.
#' @param res_clf Optional; a result object containing cross-validation or
#'   feature importance data. This is used when `importance` is `TRUE` for
#'   models that have a cross-validation feature.
#'
#' @details This function generates a plot of model coefficients or feature
#' importance, depending on the model type and the `importance` parameter. The
#' function supports SOTA models (like SVM or logistic regression), where
#' coefficients are visualized, and Random Forest models, where decision trees
#' can be plotted. The plot can also include feature importance, if available.
#'
#' If the model is a `Random Forest`, the function will plot a decision tree.
#' For other models, it will plot the feature coefficients with color-coded bars
#' for positive and negative coefficients.
#'
#' @return If the model is a `Random Forest`, the function returns a decision
#' tree plot. For other models, it returns a `ggplot` object showing the feature
#' coefficients or importance.
#'
#' @examples
#' # Example usage for a classification model
#' model <- train(logistic_regression_model)  # Assume this is a pre-trained model
#' X <- data.frame(feature1 = rnorm(100), feature2 = rnorm(100))
#' y <- sample(c(1, -1), 100, replace = TRUE)
#'
#' # Plot the model's coefficients
#' plotModel(model, X, y, main = "Model Coefficients Plot")
#'
#' # Plot the model's importance (if available)
#' plotModel(model, X, y, importance = TRUE)
#'
#' @author Edi Prifti (IRD)
#'
#' @import ggplot2
#' @importFrom reshape2 melt
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


#' Plot Model Performance Scores
#'
#' This function visualizes the performance of a model by plotting the predicted
#' scores (`yhat`) against the true class labels (`y`). It supports both
#' classification (AUC, accuracy) and regression (R2, correlation) tasks. For
#' classification tasks, it displays a boxplot with jittered points, while for
#' regression tasks, it shows a scatter plot with a linear regression line.
#'
#' @param mod A model object containing the attribute `score_` (predicted
#'   scores) and other model evaluation metrics such as accuracy, AUC, R2, etc.
#' @param y A vector of true class labels or continuous values, corresponding to
#'   the predicted scores in `mod$score_`.
#' @param col.sign A vector of two colors for positive and negative class labels
#'   (default is `c("deepskyblue1", "firebrick1")`).
#' @param main A string for the title of the plot (default is an empty string).
#'
#' @details This function checks the validity of the model and the input data,
#' then creates a plot based on the model's prediction performance. For
#' classification tasks, it uses a boxplot to show the distribution of predicted
#' scores, while for regression tasks, it uses a scatter plot with a linear
#' regression line.
#'
#' The function also displays performance metrics in the plot title, such as
#' accuracy and AUC for classification tasks, or correlation coefficient (Rho),
#' R-squared (R2), and standard error of regression (SER) for regression tasks.
#'
#' If the model is of type `SOTA` or lacks the required attributes (`score_`,
#' `y`), the function will return `NULL` and display a corresponding error
#' message.
#'
#' @return A `ggplot` object displaying the model's performance plot, either a
#' boxplot for classification tasks or a scatter plot with a regression line for
#' regression tasks.
#'
#' @examples
#' # Example usage for a classification model
#' model <- train(logistic_regression_model)  # Assume this is a pre-trained model
#' X <- data.frame(feature1 = rnorm(100), feature2 = rnorm(100))
#' y <- sample(c(1, -1), 100, replace = TRUE)
#'
#' # Plot the model's performance score
#' plotModelScore(model, y, main = "Classification Model Performance")
#'
#' # Example usage for a regression model
#' model <- train(regression_model)  # Assume this is a pre-trained model
#' y <- rnorm(100)  # Continuous response variable
#'
#' # Plot the model's performance score
#' plotModelScore(model, y, main = "Regression Model Performance")
#'
#' @author Edi Prifti (IRD)
#'
#' @import ggplot2
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



#' Normalize Model Coefficients
#'
#' This function normalizes the model coefficients to a common scale. It is
#' designed to handle different types of models and normalizes their
#' coefficients for comparison purposes. The function also offers the ability to
#' sort the features based on their importance or discriminatory power relative
#' to the target variable.
#'
#' @param mod A model object that includes the attribute `coeffs_` (the model
#'   coefficients). The model should also include an associated `language`
#'   attribute specifying the type of model (e.g., `bin`, `glmnet`, `svm`).
#' @param X A matrix or data frame of feature values used in the model.
#' @param y A vector of true labels or target values corresponding to `X`.
#' @param sort.features A logical value indicating whether to sort the features
#'   based on their importance in relation to the target variable (default is
#'   `FALSE`).
#' @param sort.ind A vector of indices for sorting the features. If `NULL`, the
#'   function will determine the sorting based on feature discriminance (default
#'   is `NULL`).
#'
#' @details The function normalizes the coefficients of a model based on the
#' type of model. It handles models such as logistic regression (`logreg`), SVM
#' (`svm`), GLM (`glmnet`), and others. If `sort.features` is set to `TRUE`, the
#' features will be sorted based on their discriminatory power, as calculated by
#' a feature selection method (e.g., filtering based on statistical
#' significance).
#'
#' The normalized coefficients are scaled to lie between -1 and 1, which allows
#' for a fair comparison of feature importance across different models.
#'
#' @return A numeric vector containing the normalized coefficients of the model,
#' or `NULL` if the model does not have valid coefficients.
#'
#' @examples
#' # Example usage for a logistic regression model
#' model <- train(logistic_regression_model)  # Assume this is a pre-trained model
#' X <- data.frame(feature1 = rnorm(100), feature2 = rnorm(100))
#' y <- sample(c(1, -1), 100, replace = TRUE)
#'
#' # Normalize the model coefficients
#' normalized_coeffs <- normModelCoeffs(model, X, y)
#'
#' @author Edi Prifti (IRD)
#'
#' @import ggplot2
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

#' Dissects the model by separating positive and negative coefficients
#'
#' This function dissects a model by separating its positive and negative
#' coefficients, calculates the corresponding scores for each group (positive
#' and negative coefficients), and normalizes them. It also provides a plot
#' showing the composition of the score.
#'
#' @param mod A valid model object.
#' @param X The matrix of features (design matrix).
#' @param y The class labels (response variable).
#' @param clf The classifier used (not currently utilized in the function).
#' @param plot Logical, if `TRUE`, a plot will be generated showing the score
#'   composition and a classification of samples based on the score.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{mod}{The provided model.}
#'   \item{y}{The response variable.}
#'   \item{scores}{A matrix containing positive, negative, and raw scores.}
#'   \item{scores.norm}{Normalized scores.}
#' }
#'
#' @details The function works by first identifying the positive and negative
#'   coefficients from the model. It then calculates the corresponding scores
#'   for both the positive and negative coefficients. The scores are normalized
#'   by dividing each score by the total sum of the scores. Finally, the
#'   function provides an optional plot that visualizes the score composition.
#'
#'   The plot shows:
#' \itemize{
#'   \item A barcode plot of the score composition.
#'   \item A classification of the samples according to the model's score with the intercept line.
#' }
#'
#' @author Edi Prifti (IRD)
#'
#' @examples
#' \dontrun{
#' # Assuming `mod`, `X`, and `y` are already defined
#' dissectResult <- disectModel(mod = mod, X = X, y = y, plot = TRUE)
#' }
#'
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



#' Plot a barcode representation of model scores
#'
#' This function creates a barcode-style heatmap to visualize the scores of a
#' model, ordered by class labels. It uses color gradients to represent the
#' score values and displays the relationship between the scores and the actual
#' classes.
#'
#' @param dscore A matrix of model scores with rows representing different
#'   features (or samples) and columns representing the model's prediction
#'   scores for each sample.
#' @param y A vector of class labels corresponding to the samples.
#' @param nb.col.levels An integer specifying the number of color levels to
#'   represent the scores. Default is 30.
#' @param main A title for the plot. Default is an empty string.
#'
#' @return This function generates a barcode-style heatmap plot.
#'
#' @details The function visualizes the model scores by reordering them
#' according to the class labels (`-1` and `1`). It uses a color gradient to
#' represent the range of scores and adds a grid for better visual distinction.
#' The plot also includes axes to indicate the feature names and the class
#' labels.
#'
#' The color palette is generated using `viridis` for better visibility of
#' scores across different value ranges. The breaks are set to cover the entire
#' range of the scores.
#'
#' @author Edi Prifti (IRD)
#'
#' @examples
#' \dontrun{
#' # Assuming `dscore` is a matrix of scores and `y` is a vector of class labels
#' plotScoreBarcode(dscore, y, nb.col.levels = 30, main = "Model Score Barcode")
#' }
#'
#' @export
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



#' Plot AUC and ROC Curve with Confidence Intervals
#'
#' This function generates a Receiver Operating Characteristic (ROC) curve and
#' computes the Area Under the Curve (AUC) along with the corresponding
#' confidence intervals. It also highlights the best threshold using Youden's
#' index.
#'
#' @param score A numeric vector containing the predicted scores from the model.
#' @param y A numeric or factor vector containing the true class labels. The
#'   labels should be binary, with two levels (e.g., 1 and -1, or 0 and 1).
#' @param main A string representing the title of the plot. Default is an empty
#'   string.
#' @param ci A logical value indicating whether to compute and display the
#'   confidence intervals for the AUC. Default is `TRUE`.
#' @param percent A logical value indicating whether to express the ROC curve in
#'   percentage scale. Default is `TRUE`.
#'
#' @return A `roc` object from the `pROC` package, containing the ROC curve and
#'   AUC information.
#'
#' @details The function uses the `pROC` package to compute and plot the ROC
#' curve and AUC. The best threshold is determined using Youden’s index, and it
#' is displayed on the plot with vertical and horizontal lines. The plot
#' includes the AUC value and its confidence intervals, as well as the best
#' threshold on the curve.
#'
#' @author Edi Prifti (IRD)
#'
#' @examples
#' \dontrun{
#' # Assuming `score` is a vector of predicted scores and `y` is the true labels
#' plotAUC(score, y, main = "ROC Curve with AUC", ci = TRUE)
#' }
#'
#' @export
plotAUC <- function(score, y, main="", ci = TRUE, percent = TRUE)
{
  require(pROC)
  rocobj <-  pROC::roc(response = y, predictor = score, percent = percent, ci = ci, of = "se", sp = seq(0, 100, 5))
  plot(rocobj, ci.type="shape", ci.col="grey80", main=main)
  # compute information on the threshold
  rocobj2 <- pROC::roc(response = y, predictor = score, percent = percent, ci = TRUE, of = "auc")
  resa  = coords(rocobj2, x = "best", input = "threshold", best.method = "youden")
  abline(v=resa[2], col="red", lty=2); abline(h=resa[3], col="red", lty=2)
  legend("bottomright",legend = c(paste("auc:",signif(rocobj$auc,3)),
                                  paste("ci:",signif(rocobj2$ci,3)[1],"-",signif(rocobj2$ci,3)[3]),
                                  paste("threshold:",signif(resa[1],3))))
  return(rocobj)
}



#' Plot AUC with ROC Curve and Confidence Intervals
#'
#' This function generates a ROC (Receiver Operating Characteristic) curve for a
#' given model or score, along with the corresponding AUC (Area Under the Curve)
#' value and its confidence intervals. Optionally, it can also display the
#' intercept point on the curve.
#'
#' @param mod An optional model object. If provided, the function will use
#'   `mod$score_` as the predicted scores. If not provided, the `score` argument
#'   must be supplied.
#' @param score A numeric vector containing the predicted scores (either
#'   provided directly or obtained from `mod`).
#' @param y A numeric or factor vector containing the true class labels. The
#'   labels should be binary (e.g., 1 and -1).
#' @param main A string representing the title of the plot. Default is an empty
#'   string.
#' @param ci A logical value indicating whether to compute and display the
#'   confidence intervals for the AUC. Default is `TRUE`.
#' @param show.intercept A logical value indicating whether to display the
#'   intercept point on the ROC curve. Default is `TRUE`.
#'
#' @return A `ggplot` object representing the ROC curve with AUC and its
#'   confidence intervals.
#'
#' @details The function computes the ROC curve and the AUC using the `pROC`
#' package. If the `mod` object is provided, the function will use `mod$score_`
#' as the predicted score. The plot includes the ROC curve, AUC, confidence
#' intervals, and optionally the intercept point. The intercept is represented
#' as a red `+` symbol on the plot.
#'
#' @author Edi Prifti (IRD)
#'
#' @examples
#' \dontrun{
#' # Assuming `mod` is a trained model and `y` is the true labels
#' plotAUCg(mod, y, main = "ROC Curve with AUC", ci = TRUE)
#' }
#'
#' @import pROC
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



# # plot a horizontal barplot
# #' @export
# plotBarplot <- function(v, rev=TRUE, xlim=range(v), main=""){
#   if(rev) v <- rev(v)
#   barplot(v, las=2, horiz=TRUE, col="black", main=main, xlim=xlim)
# }




################################################################
# PRINTING DIFFERENT, OBJECTS
################################################################

#' Print Model Information
#'
#' This function prints information about a given model, either in a short,
#' long, or structured format. The function provides a summary of the model,
#' including the model coefficients, intercept, evaluation score, learner, and
#' language, depending on the selected format.
#'
#' @param mod A model object. It can be any predomics model object.
#' @param method A string specifying the format in which the model summary will
#'   be printed. Possible values are "short" (default), "long", and "str".
#'   "short" gives a compact summary, "long" provides a detailed summary, and
#'   "str" prints the structure of the model.
#' @param score A string specifying the score attribute to be displayed. Default
#'   is "fit_".
#'
#' @return A string representing the model summary in the chosen format.
#'
#' @details
#' - The "short" method provides a brief overview of the model with information such as the coefficients,
#' intercept, decision boundary, and evaluation score (if available).
#' - The "long" method gives a more detailed version of the model summary, including the coefficients for
#' both positive and negative terms, along with other model attributes such as
#' learner type, language, and sparsity.
#' - The "str" method prints the structure of the model using the `str()` function.
#'
#' If a SOTA (state-of-the-art) model is provided, the function adjusts the
#' output accordingly, displaying the model's coefficients and attributes in a
#' simplified format.
#'
#' @examples
#' \dontrun{
#' # Assuming 'mod' is a trained model
#' printModel(mod, method = "short", score = "fit_")
#' }
#'
#' @author Edi Prifti (IRD)
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


#' Print Information about a Population of Models
#'
#' This function prints detailed information about a population of models. It
#' supports multiple methods for displaying the model summaries, such as
#' providing a "digested" view, a "short" version, or a more detailed "long"
#' view. It can also print the structure of each model within the population.
#'
#' @param obj A population of models. This should be a valid object returned by
#'   a model training procedure.
#' @param method A string specifying the format in which the population summary
#'   will be printed. Possible values are "digested", "short", "long", and
#'   "str". "digested" provides a summarized view of the population's
#'   properties, "short" gives a brief summary of each model, "long" provides a
#'   more detailed view, and "str" prints the structure of each model in the
#'   population.
#' @param score A string specifying the score attribute to be used when printing
#'   models. Default is "fit_".
#' @param indent A string used for indentation when printing information,
#'   helpful when displaying hierarchical data.
#'
#' @return None. The function prints the information directly.
#'
#' @details
#' - The "digested" method provides an overview of the population, summarizing key attributes such as the
#' sparsity and learner type.
#' - The "short" method gives a brief summary of each model, including its learner, language, and evaluation
#' score.
#' - The "long" method offers a detailed description of each model, including all relevant information about
#' coefficients and evaluation metrics.
#' - The "str" method prints the structure of each model using `str()`.
#'
#' @examples
#' \dontrun{
#' # Assuming 'population' is a valid population of models
#' printPopulation(population, method = "short")
#' }
#'
#' @author Edi Prifti (IRD)
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


#' Print Information about a Classifier Object
#'
#' This function prints detailed information about a classifier object,
#' including information about the experiment, the learner settings, and the
#' models within the classifier. It provides a structured view of the
#' classifier's parameters and any relevant attributes to facilitate
#' understanding and debugging.
#'
#' @param obj A classifier object, typically returned by a classifier training
#'   function. The object should contain details about the experiment, model
#'   parameters, and possibly a collection of models.
#' @param indent A string for indentation used when printing information. It
#'   allows for hierarchical display, making the output easier to read and
#'   understand. Default is "   --- ".
#'
#' @return None. The function prints the information directly to the console.
#'
#' @details
#' - If the classifier object contains an experiment attribute, the function prints details of the experiment.
#' - If the classifier has parameters (`params`), it prints the learner type, its parameters, and, if applicable,
#' any models used in the classifier.
#' - If the classifier includes a model collection, details of the models are printed as well.
#'
#' @examples
#' \dontrun{
#' # Assuming 'classifier' is a valid classifier object
#' printClassifier(classifier)
#' }
#'
#' @author Edi Prifti (IRD)
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

#' Print Information about an Experiment Object
#'
#' This function prints detailed information about an experiment object,
#' including details about the experiment, cross-validation, and the classifier
#' used. It helps in understanding and inspecting the components of an
#' experiment.
#'
#' @param obj An experiment object that contains information about the
#'   classifier, cross-validation, and the models used in the experiment.
#' @param indent A string for indentation, used to structure the printed output
#'   in a hierarchical manner. Default is "   --- ".
#'
#' @return None. The function prints the information directly to the console.
#'
#' @details
#' - Prints details about the experiment, including information about the classifier and the experiment settings.
#' - If cross-validation data exists, prints the number of folds, times, and seeds used in the cross-validation.
#' - Prints detailed learner options such as learner type, parameters, and the models involved in the experiment.
#'
#' @examples
#' \dontrun{
#' # Assuming 'experiment' is a valid experiment object
#' printExperiment(experiment)
#' }
#'
#' @author Edi Prifti (IRD)
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


#' Print Information about a Model Collection
#'
#' This function prints detailed information about a collection of models. It
#' allows for summarizing the model collection in either a short or long format,
#' providing insights into the models' sparsity, performance, and other relevant
#' details.
#'
#' @param obj A model collection object containing multiple models.
#' @param indent A string for indentation, used to structure the printed output
#'   in a hierarchical manner. Default is "   --- ".
#' @param method A string specifying the format for printing. Valid options are
#'   "short" and "long". Default is "long".
#'
#' @return None. The function prints the information directly to the console.
#'
#' @details
#' - In "short" mode, it prints the names of the models in the collection along with the number of models in each category.
#' - In "long" mode, it prints detailed information about each model, including the k-sparsity and a summary of the models' characteristics.
#' - The "long" mode will call the `printPopulation` function for each model in the collection to show its details.
#'
#' @examples
#' \dontrun{
#' # Assuming 'model_collection' is a valid model collection object
#' printModelCollection(model_collection, method = "short")
#' }
#'
#' @author Edi Prifti (IRD)
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


#' Print Summary of Predomics Object
#'
#' This function prints a summary of a given object, identifying its type
#' (model, population, classifier, experiment, or model collection) and calling
#' the appropriate print function to display relevant information about the
#' object.
#'
#' @param obj An object that can be of type model, population, classifier,
#'   experiment, or model collection.
#'
#' @return None. The function prints the summary of the object directly to the
#'   console.
#'
#' @details The function checks the type of the provided object using `isModel`,
#' `isPopulation`, `isClf`, `isExperiment`, and `isModelCollection` functions.
#' Based on the object type, it prints a summary:
#' - **Model**: Calls `printModel` with a detailed description of the model.
#' - **Population**: Calls `printPopulation`, showing a summary of the population of models.
#' - **Model Collection**: Calls `printModelCollection` to show a summary of a collection of models.
#' - **Experiment**: Calls `printExperiment` to display experiment details.
#' - **Classifier**: Calls `printClassifier` for classifier details.
#' If the object type is not recognized, an error message is printed.
#'
#' @examples
#' \dontrun{
#' # Assuming 'model', 'population', 'classifier', 'experiment', and 'model_collection' are valid objects
#' printy(model)
#' printy(population)
#' printy(classifier)
#' printy(experiment)
#' printy(model_collection)
#' }
#'
#' @author Edi Prifti (IRD)
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



#' Analyze Features in a Population of Models
#'
#' This function analyzes features in a population of models, allowing for the
#' visualization and examination of feature importance, prevalence, and model
#' coefficients. It can generate a variety of plots to understand the
#' distribution and importance of features in the given population.
#'
#' @param pop A population of models, typically obtained from
#'   `modelCollectionToPopulation` or similar functions.
#' @param X The data matrix containing features (rows represent features,
#'   columns represent samples).
#' @param y The response variable (class labels or continuous values depending
#'   on the model).
#' @param res_clf The classifier used for the analysis, typically a result from
#'   a classification experiment.
#' @param makeplot Logical. If `TRUE`, the function generates plots and saves
#'   them as a PDF. If `FALSE`, it returns the analysis results without
#'   plotting.
#' @param name A string representing the name of the analysis or output (used
#'   for saving files).
#' @param ord.feat A string indicating the ordering method for features. Options
#'   are:
#'   - "prevalence": Order by the prevalence of features across models.
#'   - "importance": Order by feature importance based on cross-validation.
#'   - "hierarchical": Order by hierarchical clustering of the feature-to-model coefficient matrix.
#' @param make.network Logical. If `TRUE`, generates a network of feature
#'   co-occurrence across the population of models.
#' @param network.layout A string indicating the layout of the network. Default
#'   is "circular". Other options may include "fr" for Fruchterman-Reingold
#'   layout.
#' @param network.alpha A numeric value controlling the alpha transparency of
#'   the network plot.
#' @param verbose Logical. If `TRUE`, prints additional information during
#'   execution.
#' @param pdf.dims A vector of two numbers specifying the width and height of
#'   the PDF output (in inches).
#' @param filter.perc A numeric value between 0 and 1 specifying the minimum
#'   prevalence of a feature to be included in the analysis.
#' @param k_penalty A penalty value for model selection in the population
#'   filtering.
#' @param k_max The maximum number of models to include in the final population
#'   after filtering.
#'
#' @return If `makeplot = TRUE`, returns a PDF with visualizations of feature
#'   importance, prevalence, and model coefficients. If `makeplot = FALSE`,
#'   returns a list of the analysis results including the normalized scores and
#'   feature importance.
#'
#' @details The function performs a variety of analyses on a population of
#' models:
#' - It filters models based on feature prevalence.
#' - It orders features by various metrics such as prevalence, importance, or hierarchical clustering.
#' - It generates plots of feature prevalence, model coefficients, and other characteristics.
#' - If requested, it also generates a network of feature co-occurrence across the models.
#'
#' @examples
#' \dontrun{
#' # Assuming 'pop' is a valid population of models, 'X' is the feature matrix, and 'y' is the response variable
#' analyzePopulationFeatures(pop = pop, X = X, y = y, res_clf = res_clf, makeplot = TRUE, name = "population_analysis")
#' }
#'
#' @author Edi Prifti (IRD)
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


#' Analyze Feature Importance for Machine Learning Models
#'
#' This function analyzes the importance of features in a set of machine
#' learning models. It computes various plots related to feature importance,
#' prevalence, and effect sizes. The function can handle both classification and
#' regression tasks. It can process a single experiment or multiple experiments
#' and generate corresponding visualizations in a PDF file.
#'
#' @param clf_res An object of class \code{experiment} or a list of experiments
#'   containing machine learning models to analyze.
#' @param X A data frame or matrix containing the feature data used in the
#'   model.
#' @param y A vector containing the target variable (binary or continuous values
#'   depending on the task).
#' @param makeplot Logical, if `TRUE`, plots will be generated and saved as a
#'   PDF. Default is `TRUE`.
#' @param name A string to specify the name used in output files (e.g., for
#'   saving the PDF).
#' @param verbose Logical, if `TRUE`, the function will print progress messages.
#'   Default is `TRUE`.
#' @param pdf.dims Numeric vector specifying the dimensions of the output PDF
#'   (width and height). Default is `c(width = 25, height = 20)`.
#' @param filter.perc Numeric, percentage threshold used to filter out features
#'   that appear in less than `filter.perc` of the models. Default is `0.05`
#'   (5\%).
#' @param filter.cv.prev Numeric, cross-validation threshold used to filter the
#'   importance of features based on their performance. Default is `0.25`.
#' @param nb.top.features Numeric, the number of top features to select based on
#'   importance. Default is `100`.
#' @param scaled.importance Logical, if `TRUE`, scales the feature importance
#'   scores. Default is `FALSE`.
#' @param k_penalty Numeric, penalty factor for selecting top features in
#'   models. Default is `0.75/100`.
#' @param k_max Numeric, the maximum number of features to consider. Default is
#'   `0` (no limit).
#'
#' @details This function analyzes feature importance and creates visualizations
#' of features that contribute most to the model predictions. It can handle
#' classification and regression tasks. The function computes several types of
#' graphics:
#' \itemize{
#'   \item Feature Importance: Plots the importance of features across models.
#'   \item Prevalence of Features: Shows the prevalence of features across different groups (e.g., class 1 and class -1 in classification tasks).
#'   \item Abundance of Features: Shows how frequently features appear across the dataset.
#'   \item Feature Model Coefficients: Visualizes the coefficients of features in the models.
#' }
#' The results are saved as a PDF document and also plotted directly within R.
#'
#' @return The function returns `NULL` if no models are found or after the plot
#' has been saved. It generates a PDF containing multiple plots: feature
#' importance, prevalence, abundance, and model coefficients.
#'
#' @examples
#' # Assume clf_res is a list of experiment results, and X and y are your data
#' result <- analyzeImportanceFeatures(clf_res, X, y, makeplot = TRUE, name = "Feature_Analysis", verbose = TRUE)
#'
#' # You can access the plots via result if you choose not to save them as PDFs
#'
#' @author Edi Prifti (IRD)
#'
#' @seealso \code{\link{modelCollectionToPopulation}},
#' \code{\link{plotPrevalence}}, \code{\link{plotAbundanceByClass}},
#' \code{\link{plotFeatureModelCoeffs}}
#'
#' @import ggplot2
#' @import gridExtra
#' @importFrom stats dist hclust
#' @importFrom reshape2 melt dcast
#' @importFrom cowplot ggsave2
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
makeFeatureModelPrevalenceNetworkCooccur <- function(pop.noz, 
                                                     feature.annot, 
                                                     alpha = 0.05, 
                                                     verbose = TRUE, 
                                                     layout = "circlular")
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
  # library(plyr) => took this out to force changing to tidyr
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
makeFeatureModelPrevalenceNetworkMiic <- function(pop.noz, 
                                                  feature.annot, 
                                                  cor.th = 0.3, 
                                                  verbose = TRUE, 
                                                  layout = "circlular")
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
  # library(plyr) => took this out to force changing to tidyr
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

