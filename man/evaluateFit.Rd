% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global.lib.R
\name{evaluateFit}
\alias{evaluateFit}
\title{Evaluates the fitting score of a model object}
\usage{
evaluateFit(mod, X, y, clf, force.re.evaluation = FALSE, mode = "train")
}
\arguments{
\item{mod}{: a model object}

\item{X:}{the data matrix with variables in the rows and observations in the columns}

\item{y:}{the response vector}

\item{clf:}{the classifier parameter object}

\item{force.re.evaluation:}{re-evaluate all the scores even if they exist (default:FALSE)}

\item{mode:}{A choice from c("train", "test") indicates wether we wish to learn the threthold 
of the model (default:"train") or not "test" for the c("terinter","bininter","ratio") languages}
}
\value{
a model object with the fitting score
}
\description{
Evaluates the fitting score of a model object.
}
