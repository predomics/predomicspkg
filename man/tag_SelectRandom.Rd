% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/terga2.lib.R
\name{tag_SelectRandom}
\alias{tag_SelectRandom}
\title{Randomly tag selected individuals parenting}
\usage{
tag_SelectRandom(clf, pop, nbToSelect)
}
\arguments{
\item{clf:}{The classifier object}

\item{pop:}{The population on which the selection process will be performed.}

\item{nbToSelect:}{the number of individuals we are going to select in the population}
}
\value{
the population given as an input with `nbToSelect` individuals with `selected = TRUE`
}
\description{
This function turns the selected switch on when an individual is 
selected to survive the generation and be among the pool of parents for the 
next generation.
}
