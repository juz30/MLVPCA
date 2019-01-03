EmpiricalS=function(X){
  # X need to be centered
  n=nrow(X)
  xcov=t(X)%*%X/n
  return(xcov)
}