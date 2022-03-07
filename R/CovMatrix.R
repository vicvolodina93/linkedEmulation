CovMatrix <- function(Design, cls){
  #' Constructs a squared exponential correlation matrix
  #'
  #' @param Design a matrix of inputs for design set
  #' @param cls a vector of correlation length parameters
  #'
  #' @return A matrix of correlation calculated for design points
  #'
  n = dim(Design)[1] # number of training points
  if(is.null(n))
    n <- length(Design)
  p = dim(Design)[2] # number of parameters
  if(is.null(p))
    p <- 1
  stopifnot(p == length(cls))
  if(!is.numeric(as.matrix(Design))){
    for (i in 1:p){
      Design[,i] = as.numeric(as.character(Design[,i]))
    }
  }
  D = as.matrix(dist(scale(Design,center=FALSE,scale=cls),method="euclidean",diag=TRUE,upper=TRUE))
  return(exp(-(D^2)))
}
