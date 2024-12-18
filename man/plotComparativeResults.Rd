% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global.visu.R
\name{plotComparativeResults}
\alias{plotComparativeResults}
\title{Plot Comparative Results for Multiple Methods and Cross-Validation Scores}
\usage{
plotComparativeResults(
  digested.results,
  plot = TRUE,
  ylim = c(0.5, 1),
  best = FALSE,
  ci = FALSE,
  main = "",
  mode = "classification"
)
}
\arguments{
\item{digested.results}{A list containing the performance results, including
both empirical and cross-validation (CV) scores for various methods.}

\item{plot}{A logical value (`TRUE` or `FALSE`). If `TRUE`, the function will
generate and display the plots. Default is `TRUE`.}

\item{ylim}{A numeric vector of length 2 specifying the limits for the
y-axis. Default is `c(0.5, 1)`.}

\item{best}{A logical value (`TRUE` or `FALSE`). If `TRUE`, the function will
plot the best results across methods, regardless of the k-sparsity. Default
is `FALSE`.}

\item{ci}{A logical value (`TRUE` or `FALSE`). If `TRUE`, confidence
intervals (error bars) will be shown in the plots. Default is `FALSE`.}

\item{main}{A string specifying the title of the plots. Default is an empty
string.}

\item{mode}{A string specifying the type of model being analyzed. Options are
`"classification"` or `"regression"`. Default is `"classification"`.}
}
\value{
If `plot = TRUE`, the function displays the plots. If `plot = FALSE`,
  the function returns a list of ggplot objects for further manipulation.
}
\description{
This function generates plots comparing multiple performance metrics (such as
AUC, accuracy, recall, precision, F1-score, etc.) across different methods.
It can display results from both empirical and cross-validation (CV)
evaluations, with options to show the best results across methods or by
k-sparsity values. The function can handle both classification and regression
tasks and supports visualizing both the empirical and generalization
performance.
}
\details{
The function generates multiple plots comparing performance metrics
  such as AUC, accuracy, recall, precision, F1-score, and correlation, across
  multiple methods. The plots can show:
- Empirical performance for each method.
- Cross-validation performance (generalization) for each method.
- The best results across methods, either by k-sparsity or regardless of k-sparsity.

  The plots are generated using \code{ggplot2} and arranged in a grid using
  the \code{multiplot} function. The user can choose to visualize the results
  for classification or regression models.
}
\examples{
# Assuming digested.results contains the performance scores for methods
plotComparativeResults(digested.results, plot = TRUE, ylim = c(0.5, 1),
best = TRUE, ci = TRUE, main = "Comparison of Results")

}
\author{
Edi Prifti (IRD)
}
