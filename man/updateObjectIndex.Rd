% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global.lib.R
\name{updateObjectIndex}
\alias{updateObjectIndex}
\title{updateObjectIndex}
\usage{
updateObjectIndex(obj, features = NULL)
}
\arguments{
\item{obj:}{the object can be a model, population, or modelCollection}

\item{features:}{the list of features which overrides the clf$data$features if this exists.}
}
\value{
the same object type as input, but updated
}
\description{
Update the index of a model, population, or modelCollection.
}
