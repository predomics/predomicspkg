% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global.lib.R
\name{mergeMeltScoreEmpirical}
\alias{mergeMeltScoreEmpirical}
\title{mergeMeltScoreEmpirical}
\usage{
mergeMeltScoreEmpirical(list.results.digest, k_catalogue, score = "fit_")
}
\arguments{
\item{list.results.digest:}{a list of digest objects one for each learner used. For example, list(res.terda.digest, res.terga.digest, res.terbeam.digest)}

\item{k_catalogue:}{the k_catalogue that will serve to build the result matrix}

\item{score:}{which score is to be used for value (default: fit_)}
}
\value{
a data.frame
}
\description{
mergeMeltScoreEmpirical returns a data frames that contain the performance of each digest in the list with their sparsity.
}
\details{
Merge a list of empirical scores form digest results
}
