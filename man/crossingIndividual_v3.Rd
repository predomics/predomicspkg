% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/terga2.lib.R
\name{crossingIndividual_v3}
\alias{crossingIndividual_v3}
\title{Create a New Individual by Sampling Genes from Two Parents with Controlled
Sparsity}
\usage{
crossingIndividual_v3(X, y, clf, parent1, parent2)
}
\arguments{
\item{X}{A matrix or data frame of feature values, where each row represents
a feature and each column represents a sample.}

\item{y}{A response vector or target variable for supervised learning. The
length of \code{y} should match the number of columns in \code{X}.}

\item{clf}{A classifier object containing parameters for generating a new
individual.}

\item{parent1}{A list representing the first parent, containing attributes
`indices_` (gene indices) and `eval.sparsity` (sparsity level).}

\item{parent2}{A list representing the second parent, also containing
`indices_` and `eval.sparsity`.}
}
\value{
A list representing the new individual, with a subset of genes
  sampled from both parents. If insufficient unique genes are available, the
  function raises an error.
}
\description{
This function generates a new individual by sampling genes from the combined
gene pool of two parents, with controlled sparsity based on each parent's
sparsity level. It ensures sufficient unique features in the gene pool to
meet the child's sparsity requirements.
}
\details{
The function starts by creating a combined feature reservoir
(`feat_reserv`) from the gene indices of both parents. It checks for adequate
unique features to satisfy the maximum sparsity of either parent. If either
parent lacks enough unique genes to meet its declared sparsity, the function
raises an error.

The child’s sparsity level is randomly chosen from the sparsity of one of the
parents. If the feature reservoir has duplicates, probabilities are assigned
based on feature frequency. The function then samples a subset of genes from
the reservoir for the child, ensuring the child’s sparsity matches the chosen
level.
}
\examples{
\dontrun{
X <- matrix(rnorm(100), nrow = 10) # Random features
y <- sample(c(0, 1), 10, replace = TRUE) # Random binary response
clf <- list() # Placeholder for classifier settings
parent1 <- list(indices_ = c(1, 3, 5, 7), eval.sparsity = 4)
parent2 <- list(indices_ = c(2, 4, 6, 8), eval.sparsity = 4)
child <- crossingIndividual_v3(X, y, clf, parent1, parent2)
print(child)
}

}
