% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global.visu.R
\name{plotAUC}
\alias{plotAUC}
\title{Analyze the results from a given classifier}
\usage{
plotAUC(score, y, main = "", ci = TRUE, percent = TRUE)
}
\arguments{
\item{score:}{this is the y^ of a given model}

\item{y:}{the class to be predted}

\item{main:}{title of the graph}

\item{ci:}{the point shape for the graph}

\item{percent:}{color for the graph}
}
\value{
a roc object
}
\description{
Analyze the results from a given classifier.
}
