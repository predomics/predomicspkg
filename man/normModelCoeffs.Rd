% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global.visu.R
\name{normModelCoeffs}
\alias{normModelCoeffs}
\title{Normalize the model coefficients needed for the plot}
\usage{
normModelCoeffs(mod, X, y, sort.features = FALSE, sort.ind = NULL)
}
\arguments{
\item{mod:}{a model to plot}

\item{X:}{the data matrix with variables in the rows and observations in the columns}

\item{y:}{the class vector}

\item{sort.features:}{wether the features need to be sorted by correlation with 'y' or not (default:FALSE)}

\item{sort.ind:}{computing sorting can take time if computed for every model and can be computed outside the function and passed as a parameter}
}
\value{
the normalized coefficients
}
\description{
Normalize the model coefficients needed for the plot
}
