% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minitools.R
\name{isclose}
\alias{isclose}
\title{tests weather two values are close}
\usage{
isclose(x, y, e = 1e-10)
}
\arguments{
\item{x:}{condition to be tested}

\item{y:}{message to be printed}
}
\value{
TRUE when the distance of two numbers is smaller than a given value
}
\description{
Asserts wether two vectors of the same length are close in value
  below a given threshold
}
