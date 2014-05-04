factorX <- function(X){
  r <- dim(X)[1]
  c <- dim(X)[2]
  Y <- matrix(rep(1, r*c), nrow = r, ncol= c)  
  return(X - Y)
}