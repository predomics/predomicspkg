

#' Asserts a condition and prints a message or stops the block
#'
#' @description Asserts a condition and prints a message or stops the block
#' @param condition: condition to be tested
#' @param message: message to be printed
#' @param stop: if TRUE stop the block
myAssert <- function(condition, message, stop=TRUE)
{
  if(!condition)
  {
    if(stop)
    {
      stop(message)
    }else
    {
      print(message)
    }
  }
  return(TRUE)
}

#' Asserts the existance of an object and prints a message or stops the block
#'
#' @description Asserts the existance of an object and prints a message or stops the block
#' @param obj: condition to be tested
#' @param message: message to be printed
#' @param stop: if TRUE stop the block
#' @export
myAssertNotNullNorNa <- function(obj, message="", stop=FALSE)
{
  assert <- TRUE
  if(is.null(obj))
  {
    assert <- FALSE
    if(!assert & stop)
    {
      stop(message)
    }
    return(assert)
  }
  
  if(all(is.na(obj)))
  {
    assert <- FALSE
    if(!assert & stop)
    {
      stop(message)
    }
    return(assert)
  }
  
  if(all(is.infinite(obj)))
  {
    assert <- FALSE
    if(!assert & stop)
    {
      stop(message)
    }
    return(assert)
  }
  
  return(assert)
}


#' tests weather two values are close
#'
#' @description Asserts wether two vectors of the same length are close in value below a given threshold
#' @param x: condition to be tested
#' @param y: message to be printed
#' @return TRUE when the distance of two numbers is smaller than a given value
isclose <- function(x, y, e=1e-10)
{
  if(length(x)!=length(y)){
    stop("x and y should have the same length")
  }
  
  if(length(x)>1)
  {
    return(all(abs(x-y) < e))
  }else
  {
    return(abs(x-y) < e)  
  }
}
