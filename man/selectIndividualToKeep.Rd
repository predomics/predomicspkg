% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/terga2.lib.R
\name{selectIndividualToKeep}
\alias{selectIndividualToKeep}
\title{Select the Best Individual by Sparsity Level from a Population}
\usage{
selectIndividualToKeep(clf, pop)
}
\arguments{
\item{clf}{A classifier object containing parameters for selection, including
`sparsity`, which defines the target sparsity levels.}

\item{pop}{A list representing the population, where each individual has
attributes like `eval.sparsity` (sparsity level) and `fit_` (fitness
score).}
}
\value{
A list of the best individuals for each sparsity level specified in
  `clf$params$sparsity`, selected based on fitness.
}
\description{
This function selects the best individual from a population for each
specified sparsity level, based on fitness. It groups individuals by their
sparsity, then identifies the highest-ranking individual in each group.
}
\details{
The function first groups individuals in `pop` by their sparsity
levels, as specified in `clf$params$sparsity`. For each group, it selects the
individual with the highest fitness score (`fit_`). Only sparsity levels with
at least one individual in the population are considered.

This approach allows for maintaining a diverse set of individuals across
different sparsity levels, focusing on the best-performing individual within
each sparsity category.
}
\examples{
\dontrun{
pop <- list(
  list(eval.sparsity = 3, fit_ = 0.9),
  list(eval.sparsity = 4, fit_ = 0.8),
  list(eval.sparsity = 3, fit_ = 0.85)
)
clf <- list(params = list(sparsity = c(3, 4, 5)))
best_individuals <- selectIndividualToKeep(clf, pop)
print(best_individuals)
}

}
