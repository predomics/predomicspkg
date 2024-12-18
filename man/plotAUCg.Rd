% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/global.visu.R
\name{plotAUCg}
\alias{plotAUCg}
\title{Plot AUC with ROC Curve and Confidence Intervals}
\usage{
plotAUCg(mod = NULL, score, y, main = "", ci = TRUE, show.intercept = TRUE)
}
\arguments{
\item{mod}{An optional model object. If provided, the function will use
`mod$score_` as the predicted scores. If not provided, the `score` argument
must be supplied.}

\item{score}{A numeric vector containing the predicted scores (either
provided directly or obtained from `mod`).}

\item{y}{A numeric or factor vector containing the true class labels. The
labels should be binary (e.g., 1 and -1).}

\item{main}{A string representing the title of the plot. Default is an empty
string.}

\item{ci}{A logical value indicating whether to compute and display the
confidence intervals for the AUC. Default is `TRUE`.}

\item{show.intercept}{A logical value indicating whether to display the
intercept point on the ROC curve. Default is `TRUE`.}
}
\value{
A `ggplot` object representing the ROC curve with AUC and its
  confidence intervals.
}
\description{
This function generates a ROC (Receiver Operating Characteristic) curve for a
given model or score, along with the corresponding AUC (Area Under the Curve)
value and its confidence intervals. Optionally, it can also display the
intercept point on the curve.
}
\details{
The function computes the ROC curve and the AUC using the `pROC`
package. If the `mod` object is provided, the function will use `mod$score_`
as the predicted score. The plot includes the ROC curve, AUC, confidence
intervals, and optionally the intercept point. The intercept is represented
as a red `+` symbol on the plot.
}
\examples{
\dontrun{
# Assuming `mod` is a trained model and `y` is the true labels
plotAUCg(mod, y, main = "ROC Curve with AUC", ci = TRUE)
}

}
\author{
Edi Prifti (IRD)
}
