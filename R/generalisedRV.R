# http://wwwf.imperial.ac.uk/~gmontana/software/grv/GRV_test.R

trace.matrix <- function(x)
{
  return(sum(diag(x)))
}

norm.matrix <- function(mat)
{
  return(sqrt(trace.matrix(mat%*%t(mat))))
}


create.G <- function(n,D)
{
  D <- as.matrix(D)
  A <- -(1/2)*D^2
  ones <- matrix(data=1,ncol=n,nrow=n)
  centmat <- (diag(n)-ones/n)
  G <- centmat%*%A%*%centmat
  return(G)
}

#' calculates the generalised RV coefficient of two matrices
#' @param GX matrix
#' @param GY matrix
#' @return RV coefficient
GRV <- function(GX,GY)
{
  return(trace.matrix(GX%*%GY)/(norm.matrix(GX)*norm.matrix(GY)))
}
