#' sota.glmnet: sota.glmnet classifier parameter function
#'
#' @title sota.glmnet
#' @description sota.glmnet herits from terda and does not use the randomized rounding, using thus only the glmnet component
#' @return an object containing a list of parameters for this classifier
#' @export
sota.glmnet <- function(...)
{
  # check params
  params = names(list(...))
  stopifnot (!(('language' %in% params) || ('nRR' %in% params)))
  
  # return appropriate terda object  
  clf                   <- terda(language = "logreg", nRR=0,...)
  #clf$params$intercept  <-  0
  return(clf)
}
