% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predomics.R
\name{runCrossval}
\alias{runCrossval}
\title{Compute the cross-validation emprirical and generalization scores}
\usage{
runCrossval(X, y, clf, lfolds = NULL, nfolds = 10, return.all = FALSE)
}
\arguments{
\item{X:}{the data matrix with variables in the rows and observations in the columns}

\item{y:}{the response vector}

\item{clf:}{the classifier parameter object}

\item{nfolds:}{the number of folds for the cross-validation}

\item{return.all:}{return all results from the crossvalidation for feature stability testing}
}
\value{
a list containing empirical, generalisation scores for each fold as well as a matrix with the mean values.
}
\description{
Compute the cross-validation emprirical and generalization scores.
}
