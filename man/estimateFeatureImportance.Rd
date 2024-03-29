% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global.lib.R
\name{estimateFeatureImportance}
\alias{estimateFeatureImportance}
\title{Estimates the importance of each feature in the model object}
\usage{
estimateFeatureImportance(
  mod,
  X,
  y,
  clf,
  attribute = "unpenalized_fit_",
  plot.importance = FALSE
)
}
\arguments{
\item{mod:}{a model object}

\item{X:}{the data matrix with variables in the rows and observations in the columns}

\item{y:}{the response vector}

\item{clf:}{the classifier parameter object}

\item{attribute:}{which attribute should be used to compute the importance (default:unpenalized_fit_)}

\item{plot.importance:}{should the function plot the improtance of the features (default:FALSE)}
}
\value{
a model object with the importance of each feature computed. Negative importance of a feature means that the feature is not beneficial.
}
\description{
Estimates the importance of each feature in the model object
}
